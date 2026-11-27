#include "world.hpp"

#include "../solver/particle.hpp"
#include "../solver/distance_constraints.hpp"
#include "../solver/angle_constraints.hpp"
#include "../solver/volume_constraints.hpp"
#include "../solver/actuators.hpp"
#include "../solver/collision.hpp"
#include "cluster.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <print>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace xpbd
{
World::World() : deltaTick(std::chrono::duration<double>(timeScale / tickRate))
{
    updateThread = std::thread(&World::update_loop, this);
}

World::~World()
{
    stop_update();
    if (updateThread.joinable())
        updateThread.join();
}

void World::pause_update()
{
    std::lock_guard<std::mutex> lock(updateMutex);
    updateState = UpdateState::Paused;
    updateCv.notify_all();
}

void World::resume_update()
{
    std::lock_guard<std::mutex> lock(updateMutex);
    updateState = UpdateState::Running;
    nextUpdateTimePoint = std::chrono::steady_clock::now();
    updateCv.notify_all();
}

void World::step_once_update()
{
    std::lock_guard<std::mutex> lock(updateMutex);
    updateState = UpdateState::StepOnce;
    nextUpdateTimePoint = std::chrono::steady_clock::now();
    updateCv.notify_all();
}

void World::stop_update()
{
    {
        std::lock_guard<std::mutex> lock(updateMutex);
        updateState = UpdateState::Stopped;
    }
    updateCv.notify_all();
}

void World::update_loop()
{
    std::unique_lock<std::mutex> lock(updateMutex);
    nextUpdateTimePoint = std::chrono::steady_clock::now();

    while (updateState != UpdateState::Stopped)
    {
        updateCv.wait(lock, [this]
                      { return updateState == UpdateState::Running ||
                               updateState == UpdateState::StepOnce ||
                               updateState == UpdateState::Stopped; });
        if (updateState == UpdateState::Stopped)
            break;

        lock.unlock();
        const auto start = std::chrono::steady_clock::now();
        step();
        lastStepTime = std::chrono::steady_clock::now() - start;

        {
            std::lock_guard<std::mutex> sceneLock(sceneMutex);
            auto snapshot = make_draw_snapshot();
            std::lock_guard<std::mutex> drawLock(clusterDrawMutex);
            clusterDraw = std::move(snapshot);
            newFrame = true;
        }
        execute_enqueued_commands();

        lock.lock();
        if (updateState == UpdateState::StepOnce)
        {
            updateState = UpdateState::Paused;
            continue;
        }

        nextUpdateTimePoint += deltaTick;
        updateCv.wait_until(lock, nextUpdateTimePoint,
                            [this]
                            { return updateState != UpdateState::Running; });
    }
}

void World::step()
{
    std::lock_guard<std::mutex> lock(sceneMutex);
    if (substeps == 0 || tickRate <= 0.0f)
        return;

    const float dt = 1.0f / tickRate * timeScale / substeps;

    for (size_t sub = 0; sub < substeps; ++sub)
    {
        iterate(scene.particle_storage, dt, gravity);
        reset_distance_constraints_lambdas(scene.distance_constraint_storage);
        reset_angle_constraints_lambdas(scene.angle_constraint_storage);
        reset_volume_constraints_lambdas(scene.volume_constraint_storage);

        solve_distance_constraints(scene.particle_storage, scene.distance_constraint_storage, dt);
        solve_angle_constraints(scene.particle_storage, scene.angle_constraint_storage, dt);
        solve_volume_constraints(scene.particle_storage, scene.volume_constraint_storage, dt);

        generate_collider_points_aabbs(polygonColliders_aabbs, scene.particle_storage, scene.polygon_collider_storage);
        generate_collider_points_aabbs(pointsColliders_aabbs, scene.particle_storage, scene.points_collider_storage);

        const size_t polygonCount = scene.polygon_collider_storage.size();
        const size_t pointsCount = scene.points_collider_storage.size();
        spatialHashAABB.clear();
        for (size_t i = 0; i < polygonCount; ++i)
        {
            spatialHashAABB.add_aabb(polygonColliders_aabbs[i], i);
        }

        std::vector<PolygonCache> caches(polygonCount);
        for (size_t i = 0; i < polygonCount; ++i)
        {
            const auto *polygon = scene.polygon_collider_storage.at(i);
            auto &cache = caches[i];
            const size_t n = polygon->indices.size();
            cache.positions_size = n;
            cache.positions.resize(n);
            cache.edges.resize(n);
            cache.edgesLen2.resize(n);
            for (size_t j = 0; j < n; ++j)
            {
                const auto *particle = scene.particle_storage.get(polygon->indices[j]);
                cache.positions[j] = particle ? particle->pos : glm::vec2(0.0f);
            }
            for (size_t j = 0, previous = n ? n - 1 : 0; j < n; previous = j++)
            {
                const glm::vec2 edge = cache.positions[j] - cache.positions[previous];
                cache.edges[previous] = edge;
                cache.edgesLen2[previous] = glm::dot(edge, edge);
            }
        }

        pecc.clear();
        for (size_t i = 0; i < polygonCount + pointsCount; ++i)
        {
            const bool sourceIsPolygon = i < polygonCount;
            const PointsCollider *sourceCollider = sourceIsPolygon
                                                       ? scene.polygon_collider_storage.at(i)
                                                       : scene.points_collider_storage.at(i - polygonCount);
            const AABB &sourceBox = sourceIsPolygon
                                        ? polygonColliders_aabbs[i]
                                        : pointsColliders_aabbs[i - polygonCount];

            for (size_t target = 0; target < polygonCount; ++target)
            {
                if (sourceIsPolygon && i == target)
                    continue;
                if (!sourceBox.collide(polygonColliders_aabbs[target]))
                    continue;

                const PointsCollider &points = *sourceCollider;
                const PointsCollider &polygon = *scene.polygon_collider_storage.at(target);
                const AABB &intersection = sourceBox.get_intersection(polygonColliders_aabbs[target]);

                PointsPolygonCollision collision;
                collision.points = points.indices;
                collision.polygon = polygon.indices;
                collision.static_friction = 0.5f * (points.static_friction + polygon.static_friction);
                collision.kinetic_friction = 0.5f * (points.kinetic_friction + polygon.kinetic_friction);
                collision.compliance = 0.5f * (points.compliance + polygon.compliance);
                collision.box = intersection;
                collision.cache = caches[target];
                add_point_edge_collision_constraints_from_points_polygon_collision(pecc, scene.particle_storage, collision);
            }
        }

        solve_point_edge_collision_constraints(scene.particle_storage, pecc, dt);
        update_velocities(scene.particle_storage, dt);
        apply_point_edge_collision_constraints_kinetic_friction(scene.particle_storage, pecc, dt);

        for (size_t i = 0; i < scene.angular_motor_storage.size(); ++i)
            solve_angularMotor(scene.particle_storage, *scene.angular_motor_storage.at(i), dt);
        for (size_t i = 0; i < scene.angular_servo_storage.size(); ++i)
            solve_angularServo(scene.particle_storage, *scene.angular_servo_storage.at(i), dt);
    }
}

void World::init()
{
    std::lock_guard<std::mutex> lock(sceneMutex);
    scene.clear();
    polygonColliders_aabbs.clear();
    pointsColliders_aabbs.clear();
    pecc.clear();
    spatialHashAABB.clear();
}

void World::set_tickRate(float rate)
{
    if (rate <= 0.0f)
        return;
    tickRate = rate;
    deltaTick = std::chrono::duration<double>(timeScale / rate);
}

void World::spawn_from_json(const std::string &name, const glm::vec2 &position)
{
    std::lock_guard<std::mutex> lock(sceneMutex);
    scene.spawn_entity(name, position);
}

void World::spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
{
    std::lock_guard<std::mutex> lock(sceneMutex);
    if (segments < 3)
        segments = 3;
    std::vector<ParticleId> ids;
    ids.reserve(segments);
    for (size_t i = 0; i < segments; ++i)
    {
        const float angle = 2.0f * 3.14159265358979323846f * float(i) / float(segments);
        ids.push_back(scene.add_particle(pos + radius * glm::vec2(std::cos(angle), std::sin(angle)), mass / float(segments)));
    }
    for (size_t i = 0; i < segments; ++i)
        scene.add_distance_constraint(ids[i], ids[(i + 1) % segments], compliance);
    scene.add_volume_constraint(ids, compliance);
    scene.add_polygon_collider(ids, 0.4f, 0.3f, 0.0f, 1, 1);
    // scene.add_points_collider(ids, 0.4f, 0.3f, 0.0f, 1, 1);
}

void World::enqueue_command(Command command)
{
    std::lock_guard<std::mutex> lock(commandQueueMutex);
    commandQueue.push(std::move(command));
}

void World::execute_enqueued_commands()
{
    std::queue<Command> local;
    {
        std::lock_guard<std::mutex> lock(commandQueueMutex);
        std::swap(local, commandQueue);
    }
    while (!local.empty())
    {
        Command command = std::move(local.front());
        local.pop();
        std::visit([this](const auto &arg)
                   {
                       using T = std::decay_t<decltype(arg)>;
                       if constexpr (std::is_same_v<T, init_command>)
                           init();
                       else if constexpr (std::is_same_v<T, set_tickRate_command>)
                           set_tickRate(arg.tickRate);
                       else if constexpr (std::is_same_v<T, spawn_from_json_command>)
                           spawn_from_json(arg.name, arg.position);
                       else if constexpr (std::is_same_v<T, spawn_polygon_command>)
                           spawn_polygon(arg.pos, arg.radius, arg.segments, arg.mass, arg.compliance);
                       else if constexpr (std::is_same_v<T, remove_particle_command>)
                           scene.remove_particle(arg.id);
                       else if constexpr (std::is_same_v<T, remove_entity_command>)
                           scene.remove_entity(arg.id, arg.mode);
                       else if constexpr (std::is_same_v<T, update_angular_motor_command>)
                           scene.update_angular_motor(arg.id, arg.target_rpm, arg.max_torque);
                       //
                   },
                   command);
    }
}

cluster::Cluster World::make_draw_snapshot() const
{
    auto handle_key = [](ParticleId id) -> uint64_t
    {
        return (uint64_t(id.generation) << 32) | id.index;
    };

    cluster::Cluster result;
    std::unordered_map<uint64_t, size_t> particleIndices;
    result.particles.reserve(scene.particle_storage.size());
    for (size_t i = 0; i < scene.particle_storage.size(); ++i)
    {
        auto id = scene.particle_storage.handle_at(i);
        const auto *p = scene.particle_storage.at(i);
        if (!id || !p) continue;
        particleIndices.emplace(handle_key(*id), result.particles.size());
        result.particles.push_back({p->pos, p->prev, p->vel, p->w});
    }
    auto indexOf = [&](ParticleId id) -> size_t
    {
        auto it = particleIndices.find(handle_key(id));
        return it == particleIndices.end() ? SIZE_MAX : it->second;
    };
    for (size_t i = 0; i < scene.distance_constraint_storage.size(); ++i)
    {
        const auto &c = *scene.distance_constraint_storage.at(i);
        result.distance_constraints.push_back({indexOf(c.id1), indexOf(c.id2), c.min_dist, c.max_dist, c.compressionDamping, c.extensionDamping, c.compliance, c.lambda});
    }
    for (size_t i = 0; i < scene.polygon_collider_storage.size(); ++i)
    {
        cluster::PointsCollider c{};
        const auto &src = *scene.polygon_collider_storage.at(i);
        for (auto id : src.indices) c.indices.push_back(indexOf(id));
        c.static_friction = src.static_friction;
        c.kinetic_friction = src.kinetic_friction;
        c.compliance = src.compliance;
        c.layer = src.layer;
        c.mask = src.mask;
        result.polygonColliders.push_back(std::move(c));
    }
    for (size_t i = 0; i < scene.points_collider_storage.size(); ++i)
    {
        cluster::PointsCollider c{};
        const auto &src = *scene.points_collider_storage.at(i);
        for (auto id : src.indices) c.indices.push_back(indexOf(id));
        c.static_friction = src.static_friction;
        c.kinetic_friction = src.kinetic_friction;
        c.compliance = src.compliance;
        c.layer = src.layer;
        c.mask = src.mask;
        result.pointsColliders.push_back(std::move(c));
    }
    return result;
}
} // namespace xpbd

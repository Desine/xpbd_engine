#include "world.hpp"
#include "types.hpp"
#include "particle.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"
#include "cluster.hpp"
#include "commands.hpp"

#include "xpbd/solver/particles.hpp"
#include "xpbd/solver/distance_constraints.hpp"
#include "xpbd/solver/angle_constraints.hpp"
#include "xpbd/solver/volume_constraints.hpp"
#include "xpbd/solver/actuators.hpp"
#include "xpbd/solver/collision.hpp"
#include "xpbd/solver/spatial_hashing.hpp"

#include "../json_body_loader.hpp"
#include "../xpbd.hpp"

#include "glm/glm.hpp"
#include <vector>
#include <condition_variable>
#include <chrono>
#include <queue>
#include <thread>
#include "tracy/Tracy.hpp"

namespace xpbd
{
    World::World()
        : deltaTick(std::chrono::duration<double>(1.0 / tickRate))
    {
        updateThread = std::thread(&World::update_loop, this);
    }
    World::~World()
    {
        {
            std::lock_guard<std::mutex> lk(updateMutex);
            updateState = UpdateState::Stopped;
        }
        updateCv.notify_all();
        if (updateThread.joinable())
            updateThread.join();
    }

    void World::pause_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::Paused;
        updateCv.notify_all();
    }

    void World::resume_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::Running;
        nextUpdateTimePoint = std::chrono::steady_clock::now();
        updateCv.notify_all();
    }

    void World::step_once_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::StepOnce;
        nextUpdateTimePoint = std::chrono::steady_clock::now();
        updateCv.notify_all();
    }

    void World::stop_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::Stopped;
        updateCv.notify_all();
    }

    void World::update_loop()
    {
        std::unique_lock<std::mutex> update_lk(updateMutex);
        nextUpdateTimePoint = std::chrono::steady_clock::now();

        while (updateState != UpdateState::Stopped)
        {
            updateCv.wait(update_lk, [this]
                          { return updateState == UpdateState::Running ||
                                   updateState == UpdateState::StepOnce ||
                                   updateState == UpdateState::Stopped; });
            if (updateState == UpdateState::Stopped)
                break;

            update_lk.unlock();
            auto start = std::chrono::steady_clock::now();
            step();
            lastStepTime = std::chrono::steady_clock::now() - start;
            {
                // submit draw
                std::lock_guard<std::mutex> draw_lk(clusterDrawMutex);
                clusterDraw = cluster;
                newFrame = true;
            }
            execute_enqueued_commands();
            update_lk.lock();
            if (updateState == UpdateState::StepOnce)
            {
                updateState = UpdateState::Paused;
                continue;
            }

            nextUpdateTimePoint += deltaTick;
            updateCv.wait_until(update_lk, nextUpdateTimePoint, [this]
                                { return updateState != UpdateState::Running; });
        }
    }

    void World::step()
    {
        ZoneScoped;

        size_t polygons_size = cluster.polygonColliders.size();
        size_t points_size = cluster.pointsColliders.size();

        float substep_time = 1.0f / tickRate * timeScale / substeps;

        bool reversePass_switch = 0;

        for (size_t sub = 0; sub < substeps; ++sub)
        {
            iterate(cluster.particles, substep_time, gravity);

            reset_distance_constraints_lambdas(cluster.distanceConstraints);
            reset_angle_constraints_lambdas(cluster.angleConstraints);
            reset_volume_constraints_lambdas(cluster.volumeConstraints);

            // for (size_t iter = 0; iter < iterations; ++iter) ITERATION
            //{                                                 ITERATION
            if (reversePass)
                reversePass_switch = !reversePass_switch;

            solve_distance_constraints(cluster.particles, cluster.distanceConstraints, substep_time, reversePass_switch);
            solve_angle_constraints(cluster.particles, cluster.angleConstraints, substep_time);
            solve_volume_constraints(cluster.particles, cluster.volumeConstraints, substep_time);

            polygonColliders_aabbs.clear();
            generate_collider_points_aabbs(polygonColliders_aabbs, cluster.particles, cluster.polygonColliders);
            pointsColliders_aabbs.clear();
            generate_collider_points_aabbs(pointsColliders_aabbs, cluster.particles, cluster.pointsColliders);

            spatialHashAABB.clear();
            for (size_t i = 0; i < polygons_size; ++i)
                spatialHashAABB.add_aabb(polygonColliders_aabbs[i], i);

            std::vector<PolygonCache> polygonsCache(polygons_size);
            for (size_t polygonId = 0; polygonId < polygons_size; ++polygonId)
            {
                PolygonCache &cache = polygonsCache[polygonId];
                std::vector<size_t> &polygonIndices = cluster.polygonColliders[polygonId].indices;
                size_t n = polygonIndices.size();
                cache.positions_size = n;

                cache.positions.resize(polygonIndices.size());
                for (size_t i = 0; i < n; ++i)
                    cache.positions[i] = cluster.particles[polygonIndices[i]].pos;

                cache.edges.resize(n);
                cache.edgesLen2.resize(n);
                for (size_t i = 0, j = n - 1; i < n; j = i++)
                {
                    glm::vec2 edge = cache.positions[i] - cache.positions[j];
                    cache.edges[j] = edge;
                    cache.edgesLen2[j] = glm::dot(edge, edge);
                }
            }

            // until polygons_size: polygon/polygon // after: point/polygon
            combinedPointsColliders.clear();
            combinedPointsColliders.reserve(polygons_size + points_size);
            combinedPointsColliders.insert(combinedPointsColliders.end(), cluster.polygonColliders.begin(), cluster.polygonColliders.end());
            combinedPointsColliders.insert(combinedPointsColliders.end(), cluster.pointsColliders.begin(), cluster.pointsColliders.end());

            combinedAABB.clear();
            combinedAABB.reserve(polygons_size + points_size);
            combinedAABB.insert(combinedAABB.end(), polygonColliders_aabbs.begin(), polygonColliders_aabbs.end());
            combinedAABB.insert(combinedAABB.end(), pointsColliders_aabbs.begin(), pointsColliders_aabbs.end());

            pecc.clear();
            for (size_t pid = 0; pid < combinedPointsColliders.size(); ++pid)
            {
                overlapping_aabb_ids.clear();
                spatialHashAABB.get_overlapping_aabb_ids_excludeId(overlapping_aabb_ids, combinedAABB[pid], pid);

                bool intersection_brute_force = false;
                if (intersection_brute_force)
                {
                    overlapping_aabb_ids.clear();
                    for (size_t i = 0; i < polygons_size; i++)
                        if (i != pid && combinedAABB[pid].intersects(polygonColliders_aabbs[i]))
                            overlapping_aabb_ids.emplace_back(i);
                }

                for (size_t i : overlapping_aabb_ids)
                {
                    float avgStaticFriction = (combinedPointsColliders[pid].staticFriction + cluster.polygonColliders[i].staticFriction) * 0.5f;
                    float avgKineticFriction = (combinedPointsColliders[pid].kineticFriction + cluster.polygonColliders[i].kineticFriction) * 0.5f;
                    float avgCompliance = (combinedPointsColliders[pid].compliance + cluster.polygonColliders[i].compliance) * 0.5f;

                    PointsPolygonCollision pointsPolygonCollision{
                        combinedPointsColliders[pid].indices,
                        cluster.polygonColliders[i].indices,
                        avgStaticFriction,
                        avgKineticFriction,
                        avgCompliance,
                        combinedAABB[pid].get_intersection(polygonColliders_aabbs[i]),
                        polygonsCache[i],
                    };
                    add_point_edge_collision_constraints_from_points_polygon_collision(pecc, cluster.particles, pointsPolygonCollision);
                }
            }

            solve_point_edge_collision_constraints(cluster.particles, pecc, substep_time);
            //} ITERATION

            update_velocities(cluster.particles, substep_time);
            apply_point_edge_collision_constraints_kinetic_friction(cluster.particles, pecc, substep_time);
            // apply_distance_constraints_damping(particles, distanceConstraints, substep_time);

            for (auto &motor : cluster.angularMotors)
            {
                solve_angularMotor(
                    cluster.particles,
                    motor,
                    substep_time);
            }
            for (auto &servo : cluster.angularServos)
            {
                solve_angularServo(
                    cluster.particles,
                    servo,
                    substep_time);
            }
        }
    }

    void World::init()
    {
        cluster.particles.clear();
        cluster.particles.reserve(1000);

        cluster.distanceConstraints.clear();
        cluster.distanceConstraints.reserve(1000);
        cluster.angleConstraints.clear();
        cluster.angleConstraints.reserve(1000);
        cluster.volumeConstraints.clear();
        cluster.volumeConstraints.reserve(100);

        cluster.polygonColliders.clear();
        cluster.polygonColliders.reserve(100);
        cluster.pointsColliders.clear();
        cluster.pointsColliders.reserve(100);

        cluster.angularMotors.clear();
        cluster.angularMotors.reserve(100);
        cluster.angularServos.clear();
        cluster.angularServos.reserve(100);

        polygonColliders_aabbs.reserve(1000);
        pointsColliders_aabbs.reserve(1000);
        combinedPointsColliders.reserve(1000);
        combinedAABB.reserve(1000);
        pecc.reserve(1000);
        overlapping_aabb_ids.reserve(1000);
    }

    void World::set_tickRate(size_t tickRate)
    {
        deltaTick = std::chrono::duration<double>(1.0 / tickRate);
    }
    void World::spawn_from_json(const std::string &name, const glm::vec2 &position)
    {
        std::lock_guard<std::mutex> clusterAdd_lk(clusterAddMutex);
        json_body_loader::load(clusterAdd, name, position);

        std::lock_guard<std::mutex> cluster_lk(clusterMutex);
        merge_clusters(cluster, clusterAdd);
    }
    void World::spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
    {
        std::lock_guard<std::mutex> clusterAdd_lk(clusterAddMutex);
        xpbd::spawn_polygon(clusterAdd, pos, radius, segments, mass, compliance);

        std::lock_guard<std::mutex> cluster_lk(clusterMutex);
        merge_clusters(cluster, clusterAdd);
    }

    void World::enqueue_command(Command variant)
    {
        std::lock_guard<std::mutex> lk(commandQueueMutex);
        commandQueue.push(variant);
    }
    void World::execute_enqueued_commands()
    {
        std::queue<Command> local;
        {
            std::lock_guard<std::mutex> lk(commandQueueMutex);
            std::swap(local, commandQueue);
        }

        while (!local.empty())
        {
            Command cmd = std::move(local.front());
            local.pop();

            std::visit([this](auto &&arg)
                       {
                using T = std::decay_t<decltype(arg)>;

                if constexpr (std::is_same_v<T, init_command>)
                    init();
                else if constexpr (std::is_same_v<T, set_tickRate_command>)
                    set_tickRate(arg.tickRate);
                else if constexpr (std::is_same_v<T, add_particle_command>)
                    add_particle(cluster, arg.pos, arg.mass, arg.vel);
                else if constexpr (std::is_same_v<T, add_distance_constraint_command>)
                    add_distance_constraint(cluster, arg.constraint.i1, arg.constraint.i2, arg.constraint.compliance, arg.constraint.minDist, arg.constraint.maxDist);
                else if constexpr (std::is_same_v<T, add_angle_constraint_command>)
                    add_angle_constraint(cluster, arg.constraint.i1, arg.constraint.i2, arg.constraint.i3, arg.constraint.i4, arg.constraint.compliance, arg.constraint.minAngle, arg.constraint.maxAngle);
                else if constexpr (std::is_same_v<T, add_volume_constraint_command>)
                    add_volume_constraint(cluster, arg.constraint.indices, arg.constraint.compliance, arg.constraint.restVolume);
                
                else if constexpr (std::is_same_v<T, add_polygon_collider_command>)
                    add_polygon_collider(cluster, arg.collider.indices, arg.collider.staticFriction, arg.collider.kineticFriction,
                                        arg.collider.compliance, arg.collider.layer, arg.collider.mask);
                else if constexpr (std::is_same_v<T, add_points_collider_command>)
                    add_points_collider(cluster, arg.collider.indices, arg.collider.staticFriction, arg.collider.kineticFriction,
                                        arg.collider.compliance, arg.collider.layer, arg.collider.mask);

                else if constexpr (std::is_same_v<T, add_angular_motor_command>)
                    add_angular_motor(cluster, arg.motor.i1, arg.motor.i2, arg.motor.i3, arg.motor.i4, arg.motor.targetRPM, arg.motor.maxTorque);
                else if constexpr (std::is_same_v<T, update_angular_motor_command>)
                    update_angular_motor(cluster, arg.index, arg.motor.targetRPM, arg.motor.maxTorque);
                
                else if constexpr (std::is_same_v<T, add_angular_servo_command>)
                    add_angular_servo(cluster, arg.servo.i1, arg.servo.i2, arg.servo.i3, arg.servo.i4, arg.servo.targetAngle, arg.servo.kp, arg.servo.kd, arg.servo.maxTorque);            
                else if constexpr (std::is_same_v<T, update_angular_servo_command>)
                    update_angular_servo(cluster, arg.index, arg.servo.targetAngle, arg.servo.kp, arg.servo.kd, arg.servo.maxTorque);
                    
                else if constexpr (std::is_same_v<T, spawn_from_json_command>)
                    spawn_from_json(arg.name, arg.position);
                else if constexpr (std::is_same_v<T, spawn_polygon_command>)
                    spawn_polygon(arg.pos, arg.radius, arg.segments, arg.mass, arg.compliance); }, cmd);
        }
    }
} // namespace xpbd
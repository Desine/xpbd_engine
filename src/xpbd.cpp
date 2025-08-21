#include "xpbd.hpp"
#include "utils.hpp"
#include "cstdio"
#include <algorithm>
#include <omp.h>
#include <thread>

#include "defines.hpp"
#include "renderer.hpp"

#include "profiler.hpp"

namespace xpbd
{
    void iterate(Particles &p, float dt, const glm::vec2 &gravity)
    {
        for (size_t i = 0; i < p.pos.size(); ++i)
        {
            if (p.w[i] == 0.0f)
                continue;

            p.vel[i] += gravity * dt;
            p.prev[i] = p.pos[i];
            p.pos[i] += p.vel[i] * dt;
        }
    }

    void update_velocities(Particles &p, float dt)
    {
        for (size_t i = 0; i < p.pos.size(); ++i)
            p.vel[i] = (p.pos[i] - p.prev[i]) / dt;
    }

    void reset_distance_constraints_lambdas(std::vector<DistanceConstraint> &dc)
    {
        for (auto &c : dc)
            c.lambda = 0;
    }
    void reset_volume_constraints_lambdas(std::vector<VolumeConstraint> &vc)
    {
        for (auto &c : vc)
            c.lambda = 0;
    }

    void solve_distance_constraint(Particles &p, DistanceConstraint &dc, float dt)
    {
        size_t i1 = dc.i1;
        size_t i2 = dc.i2;

        glm::vec2 &p1 = p.pos[i1];
        glm::vec2 &p2 = p.pos[i2];

        float w1 = p.w[i1];
        float w2 = p.w[i2];

        glm::vec2 delta = p1 - p2;
        float len = glm::length(delta);
        if (len < 1e-6f)
            return;

        float C = len - dc.restDist;
        glm::vec2 grad = delta / len;

        float alphaTilde = dc.compliance / (dt * dt);
        float denom = w1 + w2 + alphaTilde;
        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * dc.lambda) / denom;
        dc.lambda += deltaLambda;

        glm::vec2 correction = deltaLambda * grad;
        p1 += w1 * correction;
        p2 -= w2 * correction;
    }
    void solve_distance_constraints(Particles &p, std::vector<DistanceConstraint> &dc, float dt)
    {
        for (size_t i = 0; i < dc.size(); ++i)
        {
            solve_distance_constraint(p, dc[i], dt);
        }
    }
    void apply_distance_constraints_damping(Particles &p, std::vector<DistanceConstraint> &dc, float dt)
    {
        // does not work
        for (size_t i = 0; i < dc.size(); ++i)
        {
            size_t i1 = dc[i].i1;
            size_t i2 = dc[i].i2;
            if (i1 >= p.pos.size() || i2 >= p.pos.size())
                continue;

            glm::vec2 delta = p.pos[i1] - p.pos[i2];
            float len = glm::length(delta);
            if (len < 1e-6f)
                continue;

            glm::vec2 dir = delta / len;
            float rest = dc[i].restDist;
            float C = len - rest;

            float w1 = p.w[i];
            float w2 = p.w[i];
            float wSum = w1 + w2;
            if (wSum < 1e-6f)
                continue;

            glm::vec2 relVel = p.vel[i1] - p.vel[i2];
            float relVelAlongConstraint = glm::dot(relVel, dir);

            float damping = (C < 0.0f) ? dc[i].compressionDamping : dc[i].extensionDamping;

            float impulseMag = -damping * relVelAlongConstraint;
            glm::vec2 impulse = impulseMag * dir;

            p.vel[i1] += w1 / wSum * impulse;
            p.vel[i2] -= w2 / wSum * impulse;
        }
    }

    float compute_polygon_area(const std::vector<glm::vec2> &positions)
    {
        float area = 0.0f;
        size_t n = positions.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = positions[j];
            const glm::vec2 &p1 = positions[i];
            area += utils::cross_2d(p0, p1);
        }
        return 0.5f * area;
    }
    float compute_polygon_area(const Particles &p, const std::vector<size_t> &indices)
    {
        float area = 0.0f;
        size_t n = indices.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = p.pos[indices[j]];
            const glm::vec2 &p1 = p.pos[indices[i]];
            area += utils::cross_2d(p0, p1);
        }
        return 0.5f * area;
    }

    void solve_volume_constraint(Particles &p, VolumeConstraint &vc, float dt)
    {
        std::vector<size_t> &indices = vc.indices;
        size_t N = indices.size();

        float volume = compute_polygon_area(p, indices);
        float C = volume - vc.restVolume;

        float denom = 0.0f;
        std::vector<glm::vec2> grads(N);

        for (size_t i = 0; i < N; ++i)
        {
            const glm::vec2 &pi_prev = p.pos[indices[(i + N - 1) % N]];
            const glm::vec2 &pi_next = p.pos[indices[(i + 1) % N]];

            glm::vec2 grad = 0.5f * utils::perp_2d(pi_next, pi_prev);
            grads[i] = grad;

            float wi = p.w[indices[i]];
            denom += wi * glm::dot(grad, grad);
        }

        float alphaTilde = vc.compliance / (dt * dt);
        denom += alphaTilde;

        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * vc.lambda) / denom;
        vc.lambda += deltaLambda;

        for (size_t i = 0; i < N; ++i)
        {
            uint32_t idx = indices[i];
            float wi = p.w[idx];

            p.pos[idx] += wi * deltaLambda * grads[i];
        }
    }

    void solve_volume_constraints(Particles &p, std::vector<VolumeConstraint> &vc, float dt)
    {
        for (size_t i = 0; i < vc.size(); ++i)
        {
            solve_volume_constraint(p, vc[i], dt);
        }
    }

    void generate_collider_points_aabbs(
        std::vector<AABB> &out,
        const Particles &p,
        const std::vector<ColliderPoints> &colliderPoints)
    {
        out.clear();
        out.reserve(colliderPoints.size());

        for (auto &c : colliderPoints)
        {
            AABB current;
            current.l = std::numeric_limits<float>::max();
            current.r = -std::numeric_limits<float>::max();
            current.b = std::numeric_limits<float>::max();
            current.t = -std::numeric_limits<float>::max();

            for (auto i : c.indices)
            {
                glm::vec2 pos = p.pos[i];
                if (pos.x < current.l)
                    current.l = pos.x;
                if (pos.x > current.r)
                    current.r = pos.x;
                if (pos.y < current.b)
                    current.b = pos.y;
                if (pos.y > current.t)
                    current.t = pos.y;
            }

            out.emplace_back(current);
        }
    }

    bool point_in_polygon(const glm::vec2 &point, const std::vector<glm::vec2> &positions)
    {
        int windingNumber = 0;
        const size_t n = positions.size();

        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &v1 = positions[i];
            const glm::vec2 &v2 = positions[j];
            float cross = utils::cross_2d(v2 - v1, point - v1);

            if (v1.y <= point.y)
            {
                if (v2.y > point.y && cross > 0)
                    ++windingNumber;
            }
            else
            {
                if (v2.y <= point.y && cross < 0)
                    --windingNumber;
            }
        }
        return windingNumber != 0;
    }
    void get_PointEdgeCollisions_of_points_inside_polygon(
        std::vector<PointEdgeCollision> &out,
        const Particles &p,
        const std::vector<size_t> &pointIndices,
        const std::vector<size_t> &polygonIndices,
        const PolygonCache &cache)
    {
        out.reserve(pointIndices.size());

        const auto &polygonPositions = cache.positions;
        const auto &edges = cache.edges;
        const auto &edgesLen2 = cache.edgesLen2;
        size_t n = cache.positions_size;

        for (size_t pid : pointIndices)
        {
            const glm::vec2 &point = p.pos[pid];

            if (!point_in_polygon(point, polygonPositions))
                continue;

            float minDist = std::numeric_limits<float>::max();
            size_t nearestEdge = 0;

            for (size_t i = 0; i < n; ++i)
            {
                if (cache.edgesLen2[i] < 1e-12f)
                    continue;

                glm::vec2 toPoint = point - polygonPositions[i];
                float t = glm::dot(toPoint, edges[i]) / edgesLen2[i];
                t = glm::clamp(t, 0.0f, 1.0f);
                glm::vec2 closest = polygonPositions[i] + edges[i] * t;

                float dist = glm::dot(point - closest, point - closest);
                if (dist < minDist)
                {
                    minDist = dist;
                    nearestEdge = i;
                }
            }

            out.emplace_back(pid,
                             polygonIndices[nearestEdge],
                             polygonIndices[(nearestEdge + 1) % n]);
        }
    }
    void add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(
        std::vector<PointEdgeCollisionConstraints> &pecc,
        const Particles &particles,
        const PointsPolygonCollision &collision)
    {
        std::vector<size_t> filteredPointIndices;
        filteredPointIndices.reserve(collision.points.size());
        for (auto i : collision.points)
            if (collision.box.contains_point(particles.pos[i]))
                filteredPointIndices.emplace_back(i);

        if (filteredPointIndices.empty())
            return;

        std::vector<PointEdgeCollision> detections;
        get_PointEdgeCollisions_of_points_inside_polygon(detections, particles, filteredPointIndices, collision.polygon, collision.cache);

        pecc.reserve(pecc.size() + detections.size());
        for (const auto &c : detections)
            pecc.emplace_back(c.point, c.edge1, c.edge2, collision.staticFriction, collision.kineticFriction, collision.compliance, 0.0f);
    }

    void solve_point_edge_collision_constraints(
        Particles &p,
        std::vector<PointEdgeCollisionConstraints> &pecc,
        float dt)
    {
        for (size_t i = 0; i < pecc.size(); ++i)
        {
            size_t i_p = pecc[i].point;
            size_t i_e0 = pecc[i].edge1;
            size_t i_e1 = pecc[i].edge2;

            glm::vec2 &p_pos = p.pos[i_p];
            glm::vec2 &p_prev = p.prev[i_p];
            float w_p = p.w[i_p];

            glm::vec2 &e0_pos = p.pos[i_e0];
            glm::vec2 &e1_pos = p.pos[i_e1];
            glm::vec2 &e0_prev = p.prev[i_e0];
            glm::vec2 &e1_prev = p.prev[i_e1];
            float w_e0 = p.w[i_e0];
            float w_e1 = p.w[i_e1];

            glm::vec2 edge = e1_pos - e0_pos;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                return;

            float t = glm::clamp(glm::dot(p_pos - e0_pos, edge) / edgeLengthSq, 0.0f, 1.0f);
            glm::vec2 closest = e0_pos + t * edge;

            glm::vec2 n = p_pos - closest;
            float C = glm::length(n);
            if (C < 1e-6f)
                return;
            n /= C;

            glm::vec2 grad_p = n;
            glm::vec2 grad_e0 = -n * (1.0f - t);
            glm::vec2 grad_e1 = -n * t;

            float w_sum = w_p + w_e0 * (1.0f - t) * (1.0f - t) + w_e1 * t * t;
            if (w_sum < 1e-6f)
                return;

            float alpha = pecc[i].compliance / (dt * dt);
            float deltaLambda = (-C - alpha * pecc[i].lambda) / (w_sum + alpha);
            pecc[i].lambda += deltaLambda;

            p_pos += w_p * deltaLambda * grad_p;
            e0_pos += w_e0 * deltaLambda * grad_e0;
            e1_pos += w_e1 * deltaLambda * grad_e1;

            // --- Friction ---
            glm::vec2 tangent(-n.y, n.x);

            glm::vec2 p_disp = p_pos - p_prev;
            glm::vec2 e0_disp = e0_pos - e0_prev;
            glm::vec2 e1_disp = e1_pos - e1_prev;
            glm::vec2 edge_disp = e0_disp * (1.0f - t) + e1_disp * t;

            glm::vec2 rel_disp = p_disp - edge_disp;
            float tangential_disp = glm::dot(rel_disp, tangent);

            if (fabs(tangential_disp) < pecc[i].staticFriction * fabs(deltaLambda))
            {
                // static
                p_pos -= (w_p / w_sum) * tangential_disp * tangent;
                e0_pos += (w_e0 * (1.0f - t) / w_sum) * tangential_disp * tangent;
                e1_pos += (w_e1 * t / w_sum) * tangential_disp * tangent;
            }
            else
            {
                // kinetic
                float tangential_dir = (tangential_disp > 0.0f) ? 1.0f : -1.0f;
                float kineticImpulse = pecc[i].kineticFriction * fabs(deltaLambda);

                glm::vec2 frictionImpulse = -tangential_dir * kineticImpulse * tangent;

                p_pos += (w_p / w_sum) * frictionImpulse;
                e0_pos -= (w_e0 * (1.0f - t) / w_sum) * frictionImpulse;
                e1_pos -= (w_e1 * t / w_sum) * frictionImpulse;
            }
        }
    }
    void apply_point_edge_collision_constraints_kinetic_friction(
        Particles &p,
        const std::vector<PointEdgeCollisionConstraints> &pecc,
        float dt)
    {
        // does not work
        for (size_t i = 0; i < pecc.size(); ++i)
        {
            size_t p_i = pecc[i].point;
            size_t e0_i = pecc[i].edge1;
            size_t e1_i = pecc[i].edge2;

            glm::vec2 p_pos = p.pos[p_i];
            glm::vec2 &p_vel = p.vel[p_i];
            float w_p = p.w[p_i];

            glm::vec2 e0_pos = p.pos[e0_i];
            glm::vec2 e1_pos = p.pos[e1_i];
            glm::vec2 &e0_vel = p.vel[e0_i];
            glm::vec2 &e1_vel = p.vel[e1_i];
            float w_e0 = p.w[e0_i];
            float w_e1 = p.w[e1_i];

            glm::vec2 edge = e1_pos - e0_pos;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                continue;

            float t = glm::clamp(glm::dot(p_pos - e0_pos, edge) / edgeLengthSq, 0.0f, 1.0f);
            glm::vec2 closest = e0_pos + t * edge;

            glm::vec2 normal = p_pos - closest;
            float dist = glm::length(normal);
            if (dist < 1e-6f)
                continue;

            normal /= dist;
            glm::vec2 tangent(-normal.y, normal.x);

            glm::vec2 edge_vel = e0_vel * (1.0f - t) + e1_vel * t;
            glm::vec2 rel_vel = p_vel - edge_vel;
            float v_tangent = glm::dot(rel_vel, tangent);
            if (fabs(v_tangent) < 1e-6f)
                continue;

            float w_sum = w_p + w_e0 * (1.0f - t) * (1.0f - t) + w_e1 * t * t;
            if (w_sum < 1e-6f)
                continue;

            float maxFrictionImpulse = pecc[i].kineticFriction * pecc[i].lambda;

            float deltaImpulse = -v_tangent / (w_sum * dt);
            deltaImpulse = glm::clamp(deltaImpulse, -maxFrictionImpulse, maxFrictionImpulse);

            glm::vec2 frictionImpulse = deltaImpulse * tangent;

            p.vel[p_i] += w_p * frictionImpulse;
            p.vel[e0_i] -= w_e0 * frictionImpulse * (1.0f - t);
            p.vel[e1_i] -= w_e1 * frictionImpulse * t;
        }
    }

    // Class World
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
    }

    void World::resume_update()
    {
        {
            std::lock_guard<std::mutex> lk(updateMutex);
            updateState = UpdateState::Running;
            nextUpdateTimePoint = std::chrono::steady_clock::now();
        }
        updateCv.notify_all();
    }

    void World::step_once_update()
    {
        {
            std::lock_guard<std::mutex> lk(updateMutex);
            updateState = UpdateState::StepOnce;
            nextUpdateTimePoint = std::chrono::steady_clock::now();
        }
        updateCv.notify_all();
    }

    void World::update_loop()
    {
        std::unique_lock<std::mutex> update_lk(updateMutex);
        nextUpdateTimePoint = std::chrono::steady_clock::now();

        while (updateState != UpdateState::Stopped)
        {
            updateCv.wait(update_lk, [this]
                          { return updateState != UpdateState::Paused &&
                                   updateState != UpdateState::Stopped; });
            if (updateState == UpdateState::Stopped)
                break;

            update_lk.unlock();
            step();
            submit_draw();
            execute_enqueued_functions();
            update_lk.lock();

            if (updateState == UpdateState::StepOnce)
            {
                updateState = UpdateState::Paused;
                continue;
            }

            nextUpdateTimePoint += deltaTick;
            updateCv.wait_until(update_lk, nextUpdateTimePoint, [this]
                                { return updateState != UpdateState::Running &&
                                         updateState != UpdateState::Stopped; });
        }
    }

    void World::step()
    {
        ZoneScoped;

        size_t polygons_size = polygonColliders.size();
        size_t points_size = pointsColliders.size();

        float substep_time = 1.0f / tickRate * timeScale / substeps;

        for (size_t sub = 0; sub < substeps; ++sub)
        {
            iterate(particles, substep_time, gravity);

            reset_distance_constraints_lambdas(distanceConstraints);
            reset_volume_constraints_lambdas(volumeConstraints);

            for (size_t iter = 0; iter < iterations; ++iter)
            {
                solve_distance_constraints(particles, distanceConstraints, substep_time);
                solve_volume_constraints(particles, volumeConstraints, substep_time);

                std::vector<AABB> polygonColliders_aabbs(polygons_size);
                generate_collider_points_aabbs(polygonColliders_aabbs, particles, polygonColliders);
                std::vector<AABB> pointsColliders_aabbs(points_size);
                generate_collider_points_aabbs(pointsColliders_aabbs, particles, pointsColliders);

                spatialHashAABB.clear();
                for (size_t i = 0; i < polygons_size; ++i)
                    spatialHashAABB.add_aabb(polygonColliders_aabbs[i], i);

            
                std::vector<PolygonCache> polygonsCache(polygons_size);
                for (size_t polygonId = 0; polygonId < polygons_size; ++polygonId)
                {
                    PolygonCache &hash = polygonsCache[polygonId];
                    std::vector<size_t> &polygonIndices = polygonColliders[polygonId].indices;
                    size_t n = polygonIndices.size();
                    hash.positions_size = n;

                    hash.positions.resize(polygonIndices.size());
                    for (size_t i = 0; i < n; ++i)
                        hash.positions[i] = particles.pos[polygonIndices[i]];

                    hash.edges.resize(n);
                    hash.edgesLen2.resize(n);
                    for (size_t i = 0, j = n - 1; i < n; j = i++)
                    {
                        glm::vec2 edge = hash.positions[i] - hash.positions[j];
                        hash.edges[j] = edge;
                        hash.edgesLen2[j] = glm::dot(edge, edge);
                    }
                }

                // until polygons_size: polygon/polygon // after: point/polygon
                std::vector<ColliderPoints> combinedColliderPoints;
                combinedColliderPoints.reserve(polygons_size + points_size);
                combinedColliderPoints.insert(combinedColliderPoints.end(), polygonColliders.begin(), polygonColliders.end());
                combinedColliderPoints.insert(combinedColliderPoints.end(), pointsColliders.begin(), pointsColliders.end());

                std::vector<AABB> combinedAABB;
                combinedAABB.reserve(polygons_size + points_size);
                combinedAABB.insert(combinedAABB.end(), polygonColliders_aabbs.begin(), polygonColliders_aabbs.end());
                combinedAABB.insert(combinedAABB.end(), pointsColliders_aabbs.begin(), pointsColliders_aabbs.end());

                std::vector<std::vector<PointEdgeCollisionConstraints>> pecc_thread(omp_get_max_threads());
#pragma omp parallel for
                for (size_t pid = 0; pid < combinedColliderPoints.size(); ++pid)
                {
                    std::vector<size_t> overlapping_aabb_ids;
                    spatialHashAABB.get_overlapping_aabb_ids_excludeId(overlapping_aabb_ids, combinedAABB[pid], pid);

                    bool intersection_brute_force = true;
                    if (intersection_brute_force)
                    {
                        overlapping_aabb_ids.clear();
                        for (size_t i = 0; i < polygons_size; i++)
                            if (i != pid && combinedAABB[pid].intersects(polygonColliders_aabbs[i]))
                                overlapping_aabb_ids.emplace_back(i);
                    }

                    auto &local_pecc = pecc_thread[omp_get_thread_num()];
                    local_pecc.reserve(overlapping_aabb_ids.size() * 10);

                    for (size_t i : overlapping_aabb_ids)
                    {
                        float avgStaticFriction = (combinedColliderPoints[pid].staticFriction + polygonColliders[i].staticFriction) * 0.5f;
                        float avgKineticFriction = (combinedColliderPoints[pid].kineticFriction + polygonColliders[i].kineticFriction) * 0.5f;
                        float avgCompliance = (combinedColliderPoints[pid].compliance + polygonColliders[i].compliance) * 0.5f;

                        PointsPolygonCollision pointsPolygonCollision{
                            combinedColliderPoints[pid].indices,
                            polygonColliders[i].indices,
                            avgStaticFriction,
                            avgKineticFriction,
                            avgCompliance,
                            combinedAABB[pid].get_intersection(polygonColliders_aabbs[i]),
                            polygonsCache[i],
                        };
                        add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(local_pecc, particles, pointsPolygonCollision);
                    }
                }

                // combine pecc from pecc_thread
                size_t total_size = 0;
                for (auto &v : pecc_thread)
                    total_size += v.size();

                std::vector<PointEdgeCollisionConstraints> pecc;
                pecc.reserve(total_size);
                for (auto &v : pecc_thread)
                    pecc.insert(pecc.end(), v.begin(), v.end());
            

                solve_point_edge_collision_constraints(particles, pecc, substep_time);
            }

            update_velocities(particles, substep_time);
            // apply_point_edge_collision_constraints_kinetic_friction(particles, pecc, substep_time);
            // apply_distance_constraints_damping(particles, distanceConstraints, substep_time);
        }
    }

    void World::enqueue_function(std::function<void()> f)
    {
        std::lock_guard<std::mutex> lk(functionQueueMutex);
        queue.push(std::move(f));
    }
    void World::execute_enqueued_functions()
    {
        std::queue<std::function<void()>> local;
        {
            std::lock_guard<std::mutex> lk(functionQueueMutex);
            std::swap(local, queue);
        }
        while (!local.empty())
        {
            local.front()();
            local.pop();
        }
    }

    void World::init()
    {
        particles.pos.clear();
        particles.pos.reserve(5000);
        particles.prev.clear();
        particles.prev.reserve(5000);
        particles.vel.clear();
        particles.vel.reserve(5000);
        particles.w.clear();
        particles.w.reserve(5000);

        distanceConstraints.clear();
        distanceConstraints.reserve(5000);

        volumeConstraints.clear();
        volumeConstraints.reserve(1000);

        polygonColliders.clear();
        polygonColliders.reserve(1000);

        pointsColliders.clear();
        pointsColliders.reserve(1000);

        // polygonColliders_aabbs.reserve(1000);
        // pointsColliders_aabbs.reserve(1000);

        // polygonsHash.reserve(1000);

        // combinedColliderPoints.reserve(2000);
        // combinedAABB.reserve(2000);

        // pecc.reserve(5000);
    }
    void World::enqueue_init()
    {
        enqueue_function([this]
                         { init(); });
    }

    void World::submit_draw()
    {
        std::lock_guard<std::mutex> draw_lk(drawMutex);
        draw.particles = particles;
        draw.distanceConstraints = distanceConstraints;
        draw.volumeConstraints = volumeConstraints;
        draw.polygonColliders = polygonColliders;
        draw.pointsColliders = pointsColliders;
        newFrame = true;
    }

    void World::set_tickRate(size_t tickRate)
    {
        deltaTick = std::chrono::duration<double>(1.0 / tickRate);
    }
    void World::enqueue_set_tickRate(size_t tickRate)
    {
        enqueue_function([this, tickRate]
                         { set_tickRate(tickRate); });
    }

    void World::add_particle(const glm::vec2 &pos, float mass, const glm::vec2 &vel = {0.0f, 0.0f})
    {
        particles.pos.emplace_back(pos);
        particles.prev.emplace_back(pos);
        particles.vel.emplace_back(vel);
        float w = mass == 0 ? 0 : 1 / mass;
        particles.w.emplace_back(w);
    }
    void World::enqueue_add_particle(const glm::vec2 &pos, float mass, const glm::vec2 &vel = {0.0f, 0.0f})
    {
        enqueue_function([this, pos, mass, vel]
                         { add_particle(pos, mass, vel); });
    }

    void World::add_distance_constraint(size_t i1, size_t i2, float compliance, float restDist = FLOAT_DEFAULT)
    {
        if (restDist == FLOAT_DEFAULT)
            restDist = glm::distance(particles.pos[i1], particles.pos[i2]);
        float compressionDamping = 0.5f;
        float extensionDamping = 0.5f;
        float lambda = 0;
        distanceConstraints.emplace_back(i1, i2, restDist, compressionDamping, extensionDamping, compliance, lambda);
    }
    void World::enqueue_add_distance_constraint(size_t i1, size_t i2, float compliance, float restDist = FLOAT_DEFAULT)
    {
        enqueue_function([this, i1, i2, compliance, restDist]
                         { add_distance_constraint(i1, i2, compliance, restDist); });
    }

    void World::add_volume_constraint(const std::vector<size_t> &indices, float compliance, float restPressure = FLOAT_DEFAULT)
    {
        float restVolume = compute_polygon_area(particles, indices);
        if (restPressure != FLOAT_DEFAULT)
            restVolume *= restPressure;
        float lambda = 0;
        volumeConstraints.emplace_back(indices, restVolume, compliance, lambda);
    }
    void World::enqueue_add_volume_constraint(std::vector<size_t> indices, float compliance, float restPressure = FLOAT_DEFAULT)
    {
        enqueue_function([this, indices = std::move(indices), compliance, restPressure]
                         { add_volume_constraint(indices, compliance, restPressure); });
    }

    void World::add_polygon_collider(const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance)
    {
        polygonColliders.emplace_back(indices, staticFriction, kineticFriction, compliance);
    }
    void World::enqueue_add_polygon_collider(std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
    {
        enqueue_function([this, indices = std::move(indices), staticFriction, kineticFriction, compliance]
                         { add_polygon_collider(indices, staticFriction, kineticFriction, compliance); });
    }

    void World::add_points_collider(const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance)
    {
        pointsColliders.emplace_back(indices, staticFriction, kineticFriction, compliance);
    }
    void World::enqueue_add_points_collider(std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
    {
        enqueue_function([this, indices = std::move(indices), staticFriction, kineticFriction, compliance] mutable
                         { add_points_collider(indices, staticFriction, kineticFriction, compliance); });
    }

    void World::spawn_from_json(const std::string &name, const glm::vec2 &position)
    {
        json_body_loader::load(*this, name, position);
    }
    void World::enqueue_spawn_from_json(std::string name, glm::vec2 position)
    {
        enqueue_function([this, name = std::move(name), position]
                         { spawn_from_json(name, position); });
    }

    void World::spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
    {
        if (segments < 3)
            segments = 3;

        size_t start = particles.pos.size();
        size_t end = start + segments;
        std::vector<size_t> ids;

        float angleStep = 2.0f * float(M_PI) / segments;
        for (int i = 0; i < segments; ++i)
        {
            float angle = i * angleStep;
            glm::vec2 dir = glm::vec2(cosf(angle), sinf(angle));
            add_particle(dir * radius + pos, mass / segments);

            ids.emplace_back(start + i);
        }

        for (int i = 0; i < segments; ++i)
            add_distance_constraint(start + i, start + (i + 1) % segments, compliance);

        add_volume_constraint(ids, compliance);

        float staticFriction = 0.4f;
        float kineticFriction = 0.3f;
        add_polygon_collider(ids, staticFriction, kineticFriction, compliance);
    }
    void World::enqueue_spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
    {
        enqueue_function([this, pos, radius, segments, mass, compliance]
                         { spawn_polygon(pos, radius, segments, mass, compliance); });
    }

}
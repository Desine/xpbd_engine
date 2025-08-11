#include "xpbd.hpp"
#include "utils.hpp"
#include "cstdio"
#include <algorithm>
#include <omp.h>

#include "renderer.hpp"

#include "Remotery.h"

namespace xpbd
{
    bool should_tick(float &sec, const float &dt)
    {
        if (sec < dt)
            return false;
        sec -= dt;
        return true;
    }

    void iterate(Particles &p, const float &dt, const glm::vec2 &gravity)
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

    void update_velocities(Particles &p, const float &dt)
    {
        for (size_t i = 0; i < p.pos.size(); ++i)
            p.vel[i] = (p.pos[i] - p.prev[i]) / dt;
    }

    void add_particle(Particles &p, const glm::vec2 &pos, const float &mass)
    {
        p.pos.push_back(pos);
        p.prev.push_back(pos);
        p.vel.push_back({0, 0});
        float w = mass == 0 ? 0 : 1 / mass;
        p.w.push_back(w);
    }
    void add_particle(Particles &p, const glm::vec2 &pos, const float &mass, const glm::vec2 &vel)
    {
        p.pos.push_back(pos);
        p.prev.push_back(pos);
        p.vel.push_back(vel);
        float w = mass == 0 ? 0 : 1 / mass;
        p.w.push_back(w);
    }

    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist)
    {
        dc.i1.push_back(i1);
        dc.i2.push_back(i2);
        dc.restDist.push_back(restDist);
        dc.compliance.push_back(compliance);
        dc.lambda.push_back(0);

        dc.compressionDamping.push_back(0.1f);
        dc.extensionDamping.push_back(0.9f);
    }
    void add_distance_constraint_auto_restDist(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, Particles &p)
    {
        dc.i1.push_back(i1);
        dc.i2.push_back(i2);
        float restDist = glm::distance(p.pos[i1], p.pos[i2]);
        dc.restDist.push_back(restDist);
        dc.compliance.push_back(compliance);
        dc.lambda.push_back(0);

        dc.compressionDamping.push_back(0.5f);
        dc.extensionDamping.push_back(0.5f);
    }

    void solve_distance_constraints(Particles &p, DistanceConstraints &dc, float dt)
    {
        for (size_t i = 0; i < dc.i1.size(); ++i)
        {
            size_t i1 = dc.i1[i];
            size_t i2 = dc.i2[i];
            if (i1 >= p.pos.size() || i2 >= p.pos.size())
            {
                printf("Warning: Invalid distance constraint indices %d or %d, max index %d\n", i1, i2, p.pos.size() - 1);
                continue;
            }

            glm::vec2 &p1 = p.pos[i1];
            glm::vec2 &p2 = p.pos[i2];

            float w1 = p.w[i1];
            float w2 = p.w[i2];

            glm::vec2 delta = p1 - p2;
            float len = glm::length(delta);
            if (len < 1e-6f)
                continue;

            float C = len - dc.restDist[i];
            glm::vec2 grad = delta / len;

            float alphaTilde = dc.compliance[i] / (dt * dt);
            float denom = w1 + w2 + alphaTilde;
            if (denom < 1e-6f)
                continue;

            float deltaLambda = (-C - alphaTilde * dc.lambda[i]) / denom;
            dc.lambda[i] += deltaLambda;

            glm::vec2 correction = deltaLambda * grad;
            p1 += w1 * correction;
            p2 -= w2 * correction;
        }
    }
    void apply_distance_constraints_damping(Particles &p, DistanceConstraints &dc, float dt)
    {
        // does not work
        for (size_t i = 0; i < dc.i1.size(); ++i)
        {
            size_t i1 = dc.i1[i];
            size_t i2 = dc.i2[i];
            if (i1 >= p.pos.size() || i2 >= p.pos.size())
                continue;

            glm::vec2 delta = p.pos[i1] - p.pos[i2];
            float len = glm::length(delta);
            if (len < 1e-6f)
                continue;

            glm::vec2 dir = delta / len;
            float rest = dc.restDist[i];
            float C = len - rest;

            float w1 = p.w[i1];
            float w2 = p.w[i2];
            float wSum = w1 + w2;
            if (wSum < 1e-6f)
                continue;

            glm::vec2 relVel = p.vel[i1] - p.vel[i2];
            float relVelAlongConstraint = glm::dot(relVel, dir);

            float damping = (C < 0.0f) ? dc.compressionDamping[i] : dc.extensionDamping[i];

            float impulseMag = -damping * relVelAlongConstraint;
            glm::vec2 impulse = impulseMag * dir;

            p.vel[i1] += w1 / wSum * impulse;
            p.vel[i2] -= w2 / wSum * impulse;
        }
    }

    float compute_polygon_area(const std::vector<glm::vec2> &positions)
    {
        float area = 0.0f;
        const size_t n = positions.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = positions[j];
            const glm::vec2 &p1 = positions[i];
            area += utils::cross_2d(p0, p1);
        }
        return 0.5f * area;
    }
    float compute_polygon_area(const Particles &p, std::vector<size_t> &indices)
    {
        float area = 0.0f;
        const size_t n = indices.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = p.pos[indices[j]];
            const glm::vec2 &p1 = p.pos[indices[i]];
            area += utils::cross_2d(p0, p1);
        }
        return 0.5f * area;
    }

    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance)
    {
        vc.indices.push_back(indices);
        vc.restVolume.push_back(compute_polygon_area(p, indices));
        vc.compliance.push_back(compliance);
        vc.lambda.push_back(0);
    }

    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance, float restPressure)
    {
        vc.indices.push_back(indices);
        vc.restVolume.push_back(compute_polygon_area(p, indices) * restPressure);
        vc.compliance.push_back(compliance);
        vc.lambda.push_back(0);
    }

    void solve_volume_constraints(Particles &p, VolumeConstraints &vc, float dt)
    {
        for (size_t id = 0; id < vc.indices.size(); ++id)
        {
            std::vector<size_t> &indices = vc.indices[id];
            size_t N = indices.size();

            float volume = compute_polygon_area(p, indices);
            float C = volume - vc.restVolume[id];

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

            float alphaTilde = vc.compliance[id] / (dt * dt);
            denom += alphaTilde;

            if (denom < 1e-6f)
                continue;

            float deltaLambda = (-C - alphaTilde * vc.lambda[id]) / denom;
            vc.lambda[id] += deltaLambda;

            for (size_t i = 0; i < N; ++i)
            {
                uint32_t idx = indices[i];
                float wi = p.w[idx];

                p.pos[idx] += wi * deltaLambda * grads[i];
            }
        }
    }

    void reset_constraints_lambdas(std::vector<float> &lambdas)
    {
        std::fill(lambdas.begin(), lambdas.end(), 0.0f);
    }

    std::vector<AABB> generate_particles_aabbs(const Particles &particles, const std::vector<std::vector<size_t>> indices)
    {
        std::vector<AABB> out;

        for (size_t o = 0; o < indices.size(); o++)
        {
            AABB current;
            current.l = std::numeric_limits<float>::max();
            current.r = -std::numeric_limits<float>::max();
            current.b = std::numeric_limits<float>::max();
            current.t = -std::numeric_limits<float>::max();
            for (size_t i = 0; i < indices[o].size(); i++)
            {
                size_t id = indices[o][i];
                glm::vec2 pos = particles.pos[id];
                if (pos.x < current.l)
                    current.l = pos.x;
                if (pos.x > current.r)
                    current.r = pos.x;
                if (pos.y < current.b)
                    current.b = pos.y;
                if (pos.y > current.t)
                    current.t = pos.y;
            }
            out.push_back(current);
        }
        return out;
    }

    void add_polygon_collider(xpbd::ColliderPoints &cc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
    {
        cc.indices.push_back(indices);
        cc.staticFriction.push_back(staticFriction);
        cc.kineticFriction.push_back(kineticFriction);
        cc.compliance.push_back(compliance);
    }
    void add_point_collider(xpbd::ColliderPoints &pc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
    {
        pc.indices.push_back(indices);
        pc.staticFriction.push_back(staticFriction);
        pc.kineticFriction.push_back(kineticFriction);
        pc.compliance.push_back(compliance);
    }

    inline bool aabbs_overlap(const AABB &a, const AABB &b)
    {
        return !(a.r < b.l || a.l > b.r || a.t < b.b || a.b > b.t);
    }
    AABB create_aabb_intersection(const AABB &a, const AABB &b)
    {
        return {
            std::max(a.l, b.l),
            std::min(a.r, b.r),
            std::max(a.b, b.b),
            std::min(a.t, b.t)};
    }
    std::vector<AABBsOverlap> create_aabbs_overlaps(const std::vector<AABB> &aabbs)
    {
        std::vector<AABBsOverlap> out;

        std::vector<std::pair<AABB, size_t>> indexed;
        for (size_t i = 0; i < aabbs.size(); ++i)
            indexed.emplace_back(aabbs[i], i);

        std::sort(indexed.begin(), indexed.end(),
                  [](const auto &a, const auto &b)
                  {
                      return a.first.l < b.first.l;
                  });

        for (size_t i = 0; i < indexed.size(); ++i)
        {
            const auto &[a1, i1] = indexed[i];

            for (size_t j = i + 1; j < indexed.size(); ++j)
            {
                const auto &[a2, i2] = indexed[j];

                if (a2.l > a1.r)
                    break;

                if (aabbs_overlap(a1, a2))
                    out.push_back({i1, i2, create_aabb_intersection(a1, a2)});
            }
        }

        return out;
    }
    std::vector<AABBsOverlap> create_aabbs_overlaps(const std::vector<AABB> &aabbs1, const std::vector<AABB> &aabbs2)
    {
        std::vector<AABBsOverlap> out;

        std::vector<std::pair<AABB, size_t>> indexed1, indexed2;
        for (size_t i = 0; i < aabbs1.size(); ++i)
            indexed1.emplace_back(aabbs1[i], i);
        for (size_t j = 0; j < aabbs2.size(); ++j)
            indexed2.emplace_back(aabbs2[j], j);

        std::sort(indexed1.begin(), indexed1.end(), [](const auto &a, const auto &b)
                  { return a.first.l < b.first.l; });
        std::sort(indexed2.begin(), indexed2.end(), [](const auto &a, const auto &b)
                  { return a.first.l < b.first.l; });

        size_t j_start = 0;

        for (const auto &[a1, i1] : indexed1)
        {
            while (j_start < indexed2.size() && indexed2[j_start].first.r < a1.l)
                ++j_start;

            for (size_t j = j_start; j < indexed2.size(); ++j)
            {
                const auto &[a2, i2] = indexed2[j];
                if (a2.l > a1.r)
                    break;

                if (aabbs_overlap(a1, a2))
                    out.push_back({i1, i2, create_aabb_intersection(a1, a2)});
            }
        }

        return out;
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
    std::vector<PointEdgeCollision> get_point_edge_collisions_of_points_inside_polygon(
        const Particles &p,
        const std::vector<size_t> &pointIndices,
        const std::vector<size_t> &polygonIndices)
    {
        std::vector<PointEdgeCollision> out;
        out.reserve(pointIndices.size());

        std::vector<glm::vec2> polygonPositions;
        polygonPositions.reserve(polygonIndices.size());
        for (size_t pid : polygonIndices)
            polygonPositions.push_back(p.pos[pid]);

        size_t n = polygonIndices.size();

        std::vector<glm::vec2> edges(n);
        std::vector<float> edgesLen2(n);
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            glm::vec2 edge = polygonPositions[i] - polygonPositions[j];
            edges[j] = edge;
            edgesLen2[j] = glm::dot(edge, edge);
        }

        for (size_t pid : pointIndices)
        {
            const glm::vec2 &point = p.pos[pid];

            if (!point_in_polygon(point, polygonPositions))
                continue;

            float minDist = std::numeric_limits<float>::max();
            size_t nearestEdge = 0;

            for (size_t i = 0; i < n; ++i)
            {
                if (edgesLen2[i] < 1e-12f)
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

            out.push_back({pid,
                           polygonIndices[nearestEdge],
                           polygonIndices[(nearestEdge + 1) % n]});
        }

        return out;
    }

    inline bool point_in_aabb(const glm::vec2 &point, const AABB &aabb)
    {
        return point.x >= aabb.l && point.x <= aabb.r &&
               point.y >= aabb.b && point.y <= aabb.t;
    }
    std::vector<PointEdgeCollisionConstraints> get_point_edge_collision_constraints_of_point_to_polygon_colliders(
        const Particles &particles,
        const PointPolygonCollision &collision)
    {
        std::vector<PointEdgeCollisionConstraints> pecc;
        // draw overlap intersection box
        // AABB b = overlap.box;
        // renderer::set_color(sf::Color::White);
        // renderer::draw_axis_aligned_bounding_box(b.l, b.r, b.b, b.t);

        std::vector<size_t> filteredPointIndices;
        filteredPointIndices.reserve(collision.points.size());
        for (auto i : collision.points)
            if (point_in_aabb(particles.pos[i], collision.box))
                filteredPointIndices.push_back(i);

        if (filteredPointIndices.empty())
            return pecc;

        std::vector<PointEdgeCollision> detections = get_point_edge_collisions_of_points_inside_polygon(particles, filteredPointIndices, collision.polygon);

        pecc.reserve(detections.size());
        for (const auto &c : detections)
            pecc.push_back({c.point, c.edge1, c.edge2, collision.staticFriction, collision.kineticFriction, collision.compliance, 0.0f});
        return pecc;
    }

    std::vector<PointEdgeCollisionConstraints> get_point_edge_collision_constraints_of_point_to_polygon_colliders_parallel(
        const Particles &particles,
        const std::vector<PointPolygonCollision> &collisions)
    {
        std::vector<std::vector<PointEdgeCollisionConstraints>> pecc_thread(omp_get_max_threads());
#pragma omp parallel for
        for (size_t c = 0; c < collisions.size(); ++c)
        {
            auto insert = get_point_edge_collision_constraints_of_point_to_polygon_colliders(
                particles, collisions[c]);
            auto &local = pecc_thread[omp_get_thread_num()];
            local.insert(local.end(), insert.begin(), insert.end());
        }

        size_t total_size = 0;
        for (auto &v : pecc_thread)
            total_size += v.size();
        std::vector<PointEdgeCollisionConstraints> pecc;
        pecc.reserve(total_size);

        for (auto &v : pecc_thread)
            pecc.insert(pecc.end(), v.begin(), v.end());
        return pecc;
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
    void World::init()
    {
        particles = Particles();
        distanceConstraints = DistanceConstraints();
        volumeConstraints = VolumeConstraints();
        polygonColliders = ColliderPoints();
        pointColliders = ColliderPoints();

        spawnFromJson("ground", {0, 0});
    }

    void World::spawnFromJson(const std::string &name, const glm::vec2 &position)
    {
        json_body_loader::load(*this, name, position);
    }

    void World::addPolygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
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
            add_particle(particles, dir * radius + pos, mass / segments);

            ids.push_back(start + i);
        }

        for (int i = 0; i < segments; ++i)
            add_distance_constraint_auto_restDist(distanceConstraints, start + i, start + (i + 1) % segments, compliance, particles);

        add_volume_constraint(particles, volumeConstraints, ids, compliance);

        add_polygon_collider(polygonColliders, ids, 0.4f, 0.3f, compliance);
    }

    void World::update(float realDelta)
    {
        if (!paused)
            sec += realDelta;

        float substep_time = deltaTick * timeScale / substeps;
        while (!paused && should_tick(sec, deltaTick))
        {
            for (size_t s = 0; s < substeps; s++)
            {
                iterate(particles, substep_time, gravity);

                reset_constraints_lambdas(distanceConstraints.lambda);
                reset_constraints_lambdas(volumeConstraints.lambda);

                for (size_t i = 0; i < iterations; i++)
                {
                    rmt_BeginCPUSample(solve_distance_constraints, 0);
                    solve_distance_constraints(particles, distanceConstraints, substep_time);
                    rmt_EndCPUSample();
                    rmt_BeginCPUSample(solve_volume_constraints, 0);
                    solve_volume_constraints(particles, volumeConstraints, substep_time);
                    rmt_EndCPUSample();

                    rmt_BeginCPUSample(aabbs_overlaps, 0);
                    // polygon/polygon
                    std::vector<AABB> aabbs_polygons = generate_particles_aabbs(particles, polygonColliders.indices);
                    std::vector<AABBsOverlap> aabbs_polygon_polygon_overlaps = create_aabbs_overlaps(aabbs_polygons);
                    // point/polygon
                    std::vector<AABB> aabbs_points = generate_particles_aabbs(particles, pointColliders.indices);
                    std::vector<AABBsOverlap> aabbs_point_polygon_overlaps = create_aabbs_overlaps(aabbs_points, aabbs_polygons);
                    rmt_EndCPUSample();

                    std::vector<PointPolygonCollision> collisions;
                    collisions.reserve(aabbs_polygon_polygon_overlaps.size() * 2 + aabbs_point_polygon_overlaps.size());

                    for (auto o : aabbs_polygon_polygon_overlaps)
                    {
                        float avgStaticFriction = (polygonColliders.staticFriction[o.i1] + polygonColliders.staticFriction[o.i2]) * 0.5f;
                        float avgKineticFriction = (polygonColliders.kineticFriction[o.i1] + polygonColliders.kineticFriction[o.i2]) * 0.5f;
                        float avgCompliance = (polygonColliders.compliance[o.i1] + polygonColliders.compliance[o.i2]) * 0.5f;
                        collisions.push_back({
                            polygonColliders.indices[o.i1],
                            polygonColliders.indices[o.i2],
                            avgStaticFriction,
                            avgKineticFriction,
                            avgCompliance,
                            o.box,
                        });
                        collisions.push_back({
                            polygonColliders.indices[o.i2],
                            polygonColliders.indices[o.i1],
                            avgStaticFriction,
                            avgKineticFriction,
                            avgCompliance,
                            o.box,
                        });
                    }
                    for (auto o : aabbs_point_polygon_overlaps)
                    {
                        float avgStaticFriction = (pointColliders.staticFriction[o.i1] + polygonColliders.staticFriction[o.i2]) * 0.5f;
                        float avgKineticFriction = (pointColliders.kineticFriction[o.i1] + polygonColliders.kineticFriction[o.i2]) * 0.5f;
                        float avgCompliance = (pointColliders.compliance[o.i1] + polygonColliders.compliance[o.i2]) * 0.5f;
                        collisions.push_back({
                            pointColliders.indices[o.i1],
                            polygonColliders.indices[o.i2],
                            avgStaticFriction,
                            avgKineticFriction,
                            avgCompliance,
                            o.box,
                        });
                    }

                    this->collisions = collisions;

                    rmt_BeginCPUSample(polygon_collision_detection, 0);
                    std::vector<PointEdgeCollisionConstraints> pecc = get_point_edge_collision_constraints_of_point_to_polygon_colliders_parallel(
                        particles, collisions);
                    rmt_EndCPUSample();

                    rmt_BeginCPUSample(solve_point_edge_collision_constraints, 0);
                    solve_point_edge_collision_constraints(particles, pecc, substep_time);
                    rmt_EndCPUSample();
                }

                update_velocities(particles, substep_time);
                // apply_point_edge_collision_constraints_kinetic_friction(particles, pecc, substep_time);
                // apply_distance_constraints_damping(particles, distanceConstraints, substep_time);
            }
        }
    }

}
#include "xpbd.hpp"
#include "utils.hpp"
#include "cstdio"
#include <algorithm>

namespace xpbd
{
    bool should_tick(float &sec, const float dt)
    {
        if (sec < dt)
            return false;
        sec -= dt;
        return true;
    }

    void iterate(Particles &p, float dt, glm::vec2 gravity)
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

    void add_particle(Particles &p, glm::vec2 pos, float mass)
    {
        p.pos.push_back(pos);
        p.prev.push_back(pos);
        p.vel.push_back({0, 0});
        float w = mass == 0 ? 0 : 1 / mass;
        p.w.push_back(w);
    }
    void add_particle(Particles &p, glm::vec2 pos, float mass, glm::vec2 vel)
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
    void solve_distance_constraints_damping(Particles &p, DistanceConstraints &dc, float dt)
    {
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
        size_t N = positions.size();
        for (size_t i = 0; i < N; ++i)
        {
            const glm::vec2 &p0 = positions[i];
            const glm::vec2 &p1 = positions[(i + 1) % N];
            area += utils::cross_2d(p0, p1);
        }
        return 0.5f * area;
    }
    float compute_polygon_area(const Particles &p, std::vector<size_t> &indices)
    {
        float area = 0.0f;
        size_t N = indices.size();
        for (size_t i = 0; i < N; ++i)
        {
            const glm::vec2 &p0 = p.pos[indices[i]];
            const glm::vec2 &p1 = p.pos[indices[(i + 1) % N]];
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
                    current.l = particles.pos[id].x;
                if (pos.x > current.r)
                    current.r = particles.pos[id].x;
                if (pos.y < current.b)
                    current.b = particles.pos[id].y;
                if (pos.y > current.t)
                    current.t = particles.pos[id].y;
            }
            out.push_back(current);
        }
        return out;
    }

    void add_polygon_collider(xpbd::PolygonColliders &cc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
    {
        cc.indices.push_back(indices);
        cc.staticFriction.push_back(staticFriction);
        cc.kineticFriction.push_back(kineticFriction);
        cc.compliance.push_back(compliance);
    }
    void add_point_collider(xpbd::PointColliders &pc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
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
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabbs)
    {
        std::vector<AABBsIntersection> out;

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
                    out.push_back({i1, i2});
            }
        }

        return out;
    }
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabbs1, const std::vector<AABB> &aabbs2)
    {
        std::vector<AABBsIntersection> out;

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
                    out.push_back({i1, i2});
            }
        }

        return out;
    }

    bool point_in_polygon(const glm::vec2 point, const std::vector<glm::vec2> positions)
    {
        int windingNumber = 0;
        size_t n = positions.size();
        for (size_t i = 0; i < n; ++i)
        {
            const glm::vec2 &v1 = positions[i];
            const glm::vec2 &v2 = positions[(i + 1) % n];
            if (v1.y <= point.y)
            {
                if (v2.y > point.y && utils::cross_2d(v2 - v1, point - v1) > 0)
                    ++windingNumber;
            }
            else
            {
                if (v2.y <= point.y && utils::cross_2d(v2 - v1, point - v1) < 0)
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
        std::vector<glm::vec2> polygonPositions;
        for (size_t pid : polygonIndices)
            polygonPositions.push_back(p.pos[pid]);

        for (size_t pid : pointIndices)
        {
            const glm::vec2 &point = p.pos[pid];

            if (!point_in_polygon(point, polygonPositions))
                continue;

            float minDist = std::numeric_limits<float>::max();
            size_t nearestEdge = 0;
            size_t n = polygonIndices.size();

            for (size_t i = 0; i < n; ++i)
            {
                size_t id1 = polygonIndices[i];
                size_t id2 = polygonIndices[(i + 1) % n];
                const glm::vec2 &e1 = p.pos[id1];
                const glm::vec2 &e2 = p.pos[id2];

                glm::vec2 edge = e2 - e1;
                float len = glm::length(edge);
                if (len < 1e-6f)
                    continue;

                glm::vec2 dir = edge / len;
                glm::vec2 toPoint = point - e1;
                float proj = glm::clamp(glm::dot(toPoint, dir), 0.0f, len);
                glm::vec2 closest = e1 + dir * proj;
                float dist = glm::length(point - closest);

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
    void populate_constraints_from_detections(
        const std::vector<PointEdgeCollision> &detections,
        PointEdgeCollisionConstraints &pecc,
        float staticFriction,
        float kineticFriction,
        float compliance)
    {
        for (const auto &c : detections)
        {
            pecc.point.push_back(c.point);
            pecc.edge1.push_back(c.edge1);
            pecc.edge2.push_back(c.edge2);
            pecc.staticFriction.push_back(staticFriction);
            pecc.kineticFriction.push_back(kineticFriction);
            pecc.compliance.push_back(compliance);
            pecc.lambda.push_back(0.0f);
        }
    }
    void add_point_edge_collision_constraints_of_polygon_to_polygon_colliders(
        Particles &p,
        PointEdgeCollisionConstraints &pecc,
        const PolygonColliders &cc,
        const size_t cc_id1,
        const size_t cc_id2)
    {
        if (cc.indices[cc_id1].size() < 2 || cc.indices[cc_id2].size() < 2)
        {
            printf("Error polygon to polygon collision detection. ///if (cc.indices[cc_id1].size() < 2 || cc.indices[cc_id2].size() < 2)");
            return;
        }

        float avgStaticFriction = (cc.staticFriction[cc_id1] + cc.staticFriction[cc_id2]) * 0.5f;
        float avgKineticFriction = (cc.kineticFriction[cc_id1] + cc.kineticFriction[cc_id2]) * 0.5f;
        float avgCompliance = (cc.compliance[cc_id1] + cc.compliance[cc_id2]) * 0.5f;

        std::vector<PointEdgeCollision> detections1 = get_point_edge_collisions_of_points_inside_polygon(p, cc.indices[cc_id1], cc.indices[cc_id2]);
        std::vector<PointEdgeCollision> detections2 = get_point_edge_collisions_of_points_inside_polygon(p, cc.indices[cc_id2], cc.indices[cc_id1]);

        populate_constraints_from_detections(detections1, pecc, avgStaticFriction, avgKineticFriction, avgCompliance);
        populate_constraints_from_detections(detections2, pecc, avgStaticFriction, avgKineticFriction, avgCompliance);
    }
    void add_point_edge_collision_constraints_of_point_to_polygon_colliders(
        Particles &p,
        PointEdgeCollisionConstraints &pecc,
        const PointColliders &pointColliders,
        const PolygonColliders &polygonColliders,
        const size_t pointColliderId,
        const size_t polygonColliderId)
    {
        if (pointColliderId >= pointColliders.indices.size() || polygonColliderId >= polygonColliders.indices.size())
        {
            printf("Error point to polygon collision detection. ///if (pointColliderId >= pointColliders.indices.size() || polygonColliderId >= polygonColliders.indices.size())");
            return;
        }

        const std::vector<size_t> &pointIndices = pointColliders.indices[pointColliderId];
        const std::vector<size_t> &polygonIndices = polygonColliders.indices[polygonColliderId];

        if (pointIndices.empty() || polygonIndices.size() < 3)
        {
            printf("Error point to polygon collision detection. ///if (pointIndices.empty() || polygonIndices.size() < 3)");
            return;
        }

        float avgStaticFriction = (pointColliders.staticFriction[pointColliderId] + polygonColliders.staticFriction[polygonColliderId]) * 0.5f;
        float avgKineticFriction = (pointColliders.kineticFriction[pointColliderId] + polygonColliders.kineticFriction[polygonColliderId]) * 0.5f;
        float avgCompliance = (pointColliders.compliance[pointColliderId] + polygonColliders.compliance[polygonColliderId]) * 0.5f;

        std::vector<PointEdgeCollision> detections =
            get_point_edge_collisions_of_points_inside_polygon(p, pointIndices, polygonIndices);

        populate_constraints_from_detections(detections, pecc, avgStaticFriction, avgKineticFriction, avgCompliance);
    }

    void solve_point_edge_collision_constraints(
        Particles &p,
        PointEdgeCollisionConstraints &pecc,
        float dt)
    {
        for (size_t i = 0; i < pecc.point.size(); ++i)
        {

            size_t i_p = pecc.point[i];
            size_t i_e0 = pecc.edge1[i];
            size_t i_e1 = pecc.edge2[i];

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

            float alpha = pecc.compliance[i] / (dt * dt);
            float deltaLambda = (-C - alpha * pecc.lambda[i]) / (w_sum + alpha);
            pecc.lambda[i] += deltaLambda;

            p_pos += w_p * deltaLambda * grad_p;
            e0_pos += w_e0 * deltaLambda * grad_e0;
            e1_pos += w_e1 * deltaLambda * grad_e1;

            // --- Static friction ---
            glm::vec2 tangent(-n.y, n.x);

            glm::vec2 p_disp = p_pos - p_prev;
            glm::vec2 e0_disp = e0_pos - e0_prev;
            glm::vec2 e1_disp = e1_pos - e1_prev;
            glm::vec2 edge_disp = e0_disp * (1.0f - t) + e1_disp * t;

            glm::vec2 rel_disp = p_disp - edge_disp;
            float tangential_disp = glm::dot(rel_disp, tangent);

            if (fabs(tangential_disp) < pecc.staticFriction[i] * fabs(deltaLambda))
            {
                p_pos -= (w_p / w_sum) * tangential_disp * tangent;
                e0_pos += (w_e0 * (1.0f - t) / w_sum) * tangential_disp * tangent;
                e1_pos += (w_e1 * t / w_sum) * tangential_disp * tangent;
            }
        }
    }
    void apply_point_edge_collision_constraints_kinetic_friction(
        Particles &p,
        const PointEdgeCollisionConstraints &pecc,
        float dt)
    {
        for (size_t i = 0; i < pecc.point.size(); ++i)
        {
            size_t i_p = pecc.point[i];
            size_t i_e0 = pecc.edge1[i];
            size_t i_e1 = pecc.edge2[i];

            glm::vec2 &p_vel = p.vel[i_p];
            glm::vec2 &e0_vel = p.vel[i_e0];
            glm::vec2 &e1_vel = p.vel[i_e1];

            float w_p = p.w[i_p];
            float w_e0 = p.w[i_e0];
            float w_e1 = p.w[i_e1];

            // Положение в текущем кадре
            glm::vec2 p_pos = p.pos[i_p];
            glm::vec2 e0_pos = p.pos[i_e0];
            glm::vec2 e1_pos = p.pos[i_e1];

            glm::vec2 edge = e1_pos - e0_pos;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                continue;

            float t = glm::clamp(glm::dot(p_pos - e0_pos, edge) / edgeLengthSq, 0.0f, 1.0f);
            glm::vec2 closest = e0_pos + t * edge;

            glm::vec2 n = p_pos - closest;
            float dist = glm::length(n);
            if (dist < 1e-6f)
                continue;

            n /= dist;
            glm::vec2 tangent(-n.y, n.x);

            // относительная скорость точки к отрезку
            glm::vec2 edge_vel = e0_vel * (1.0f - t) + e1_vel * t;
            glm::vec2 rel_vel = p_vel - edge_vel;

            float vn = glm::dot(rel_vel, n);
            glm::vec2 vn_vec = vn * n;
            glm::vec2 vt_vec = rel_vel - vn_vec;
            float vt_len = glm::length(vt_vec);

            if (vt_len < 1e-6f)
                continue;

            glm::vec2 vt_dir = vt_vec / vt_len;

            float lambda_n = pecc.lambda[i];
            float mu_k = pecc.kineticFriction[i];

            float impulse_t = -mu_k * fabs(lambda_n / dt);
            glm::vec2 delta_v = impulse_t * vt_dir;

            float w_sum = w_p + w_e0 * (1.0f - t) * (1.0f - t) + w_e1 * t * t;
            if (w_sum < 1e-6f)
                continue;

            // Распределение изменения скорости
            p_vel += (w_p / w_sum) * delta_v;
            e0_vel -= (w_e0 * (1.0f - t) / w_sum) * delta_v;
            e1_vel -= (w_e1 * t / w_sum) * delta_v;
        }
    }

}
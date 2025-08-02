#include "xpbd.hpp"
#include "utils.hpp"
#include "cstdio"

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
    }
    void add_distance_constraint_auto_restDist(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, Particles &p)
    {
        dc.i1.push_back(i1);
        dc.i2.push_back(i2);
        float restDist = glm::distance(p.pos[i1], p.pos[i2]);
        dc.restDist.push_back(restDist);
        dc.compliance.push_back(compliance);
        dc.lambda.push_back(0);
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

            auto &p1 = p.pos[i1];
            auto &p2 = p.pos[i2];

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

            p1 += w1 * deltaLambda * grad;
            p2 -= w2 * deltaLambda * grad;
        }
    }

    float compute_polygon_area(const std::vector<glm::vec2> &p)
    {
        float area = 0.0f;
        size_t N = p.size();
        for (size_t i = 0; i < N; ++i)
        {
            const glm::vec2 &p0 = p[i];
            const glm::vec2 &p1 = p[(i + 1) % N];
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

    std::vector<AABB> generate_colliders_aabbs(const Particles &particles, const std::vector<std::vector<size_t>> indices)
    {
        std::vector<AABB> aabbs;

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
            aabbs.push_back(current);
        }
        return aabbs;
    }

    inline bool aabbs_intersect(const AABB &a, const AABB &b)
    {
        return !(a.r < b.l || a.l > b.r || a.t < b.b || a.b > b.t);
    }

    void add_contour_collider(xpbd::ContourColliders &cc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance)
    {
        cc.indices.push_back(indices);
        cc.staticFriction.push_back(staticFriction);
        cc.kineticFriction.push_back(kineticFriction);
        cc.compliance.push_back(compliance);
    }
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabb)
    {
        std::vector<AABBsIntersection> oc;
        for (size_t i = 0; i < aabb.size(); ++i)
            for (size_t j = i + 1; j < aabb.size(); ++j)
                if (aabbs_intersect(aabb[i], aabb[j]))
                    oc.push_back({i, j});
        return oc;
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

    void add_point_edge_collision_constraints(Particles &p, PointEdgeCollisionConstraints &pecc, ContourColliders &cc, size_t cc_id1, size_t cc_id2)
    {
        if (cc.indices[cc_id1].size() < 2 || cc.indices[cc_id2].size() < 2)
            return;

        auto detect = [&](const size_t &test, const size_t &target)
        {
            std::vector<glm::vec2> targetPositions;
            for (size_t pid : cc.indices[target])
                targetPositions.push_back(p.pos[pid]);

            for (size_t pid : cc.indices[test])
            {
                const glm::vec2 &point = p.pos[pid];

                if (!point_in_polygon(point, targetPositions))
                    continue;

                float minDist = std::numeric_limits<float>::max();
                size_t nearestEdge = 0;
                size_t n = cc.indices[target].size();

                for (size_t i = 0; i < n; ++i)
                {
                    size_t id1 = cc.indices[target][i];
                    size_t id2 = cc.indices[target][(i + 1) % n];
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

                pecc.point.push_back(pid);
                pecc.edge1.push_back(cc.indices[target][nearestEdge]);
                pecc.edge2.push_back(cc.indices[target][(nearestEdge + 1) % cc.indices[target].size()]);
                pecc.staticFriction.push_back((cc.staticFriction[cc_id1] + cc.staticFriction[cc_id2]) * 0.5f);
                pecc.kineticFriction.push_back((cc.kineticFriction[cc_id1] + cc.kineticFriction[cc_id2]) * 0.5f);
                pecc.compliance.push_back((cc.compliance[cc_id1] + cc.compliance[cc_id2]) * 0.5f);
                pecc.lambda.push_back(0);
            }
        };
        detect(cc_id1, cc_id2);
        detect(cc_id2, cc_id1);
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

            // Gradients
            glm::vec2 grad_p = n;
            glm::vec2 grad_e0 = -n * (1.0f - t);
            glm::vec2 grad_e1 = -n * t;

            float w_sum = w_p + w_e0 * (1.0f - t) * (1.0f - t) + w_e1 * t * t;
            if (w_sum < 1e-6f)
                return;

            float alpha = pecc.compliance[i] / (dt * dt);
            float deltaLambda = (-C - alpha * pecc.lambda[i]) / (w_sum + alpha);
            pecc.lambda[i] += deltaLambda;

            // Apply position correction
            p_pos += w_p * deltaLambda * grad_p;
            e0_pos += w_e0 * deltaLambda * grad_e0;
            e1_pos += w_e1 * deltaLambda * grad_e1;

            // --- Friction (static) ---
            glm::vec2 tangent(-n.y, n.x);

            glm::vec2 p_disp = p_pos - p_prev;
            glm::vec2 e0_disp = e0_pos - e0_prev;
            glm::vec2 e1_disp = e1_pos - e1_prev;
            glm::vec2 edge_disp = e0_disp * (1.0f - t) + e1_disp * t;

            glm::vec2 rel_disp = p_disp - edge_disp;
            float tangential_disp = glm::dot(rel_disp, tangent);

            float max_friction = pecc.staticFriction[i] * fabs(deltaLambda);
            float friction_correction = glm::clamp(-tangential_disp, -max_friction, max_friction);

            // Apply frictional correction
            p_pos += (w_p / w_sum) * friction_correction * tangent;
            e0_pos -= (w_e0 * (1.0f - t) / w_sum) * friction_correction * tangent;
            e1_pos -= (w_e1 * t / w_sum) * friction_correction * tangent;
        }
    }

}
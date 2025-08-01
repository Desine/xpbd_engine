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
        for (size_t i = 0; i < p.count; ++i)
        {
            if (p.w[i] == 0.0f)
                continue;

            p.vel[i] += gravity * dt;
            p.prevPos[i] = p.pos[i];
            p.pos[i] += p.vel[i] * dt;
        }
    }

    void update_velocities(Particles &p, float dt)
    {
        for (size_t i = 0; i < p.count; ++i)
            p.vel[i] = (p.pos[i] - p.prevPos[i]) / dt;
    }

    void add_particle(Particles &p, glm::vec2 pos, float mass, size_t object_id)
    {
        p.pos.push_back(pos);
        p.prevPos.push_back(pos);
        p.vel.push_back({0, 0});
        p.w.push_back(1 / mass);
        p.object_id.push_back(object_id);
        p.count++;
    }

    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist)
    {
        dc.i1.push_back(i1);
        dc.i2.push_back(i2);
        dc.restDist.push_back(restDist);
        dc.compliance.push_back(compliance);
        dc.lambda.push_back(0);
        dc.count++;
    }

    float get_distance_between_particles(const Particles &particles, size_t i1, size_t i2)
    {
        return glm::distance(particles.pos[i1], particles.pos[i2]);
    }

    void solve_distance_constraints(Particles &p, DistanceConstraints &dc, float dt)
    {
        for (size_t i = 0; i < dc.count; ++i)
        {
            size_t i1 = dc.i1[i];
            size_t i2 = dc.i2[i];
            if (i1 >= p.count || i2 >= p.count)
            {
                printf("Warning: Invalid distance constraint indices %d of %d\n", i1, i2);
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

    void reset_distance_constraints_lambdas(DistanceConstraints &dc)
    {
        std::fill(dc.lambda.begin(), dc.lambda.end(), 0.0f);
    }

    std::vector<AABB> generate_aabbs(const Particles &particles)
    {
        std::vector<AABB> aabb;

        AABB current;
        current.object_id = -1;
        for (size_t i = 0; i < particles.count; i++)
        {
            if (current.object_id != particles.object_id[i])
            {
                current.object_id = particles.object_id[i];
                current.l = FLT_MAX;
                current.r = -FLT_MAX;
                current.b = FLT_MAX;
                current.t = -FLT_MAX;
                aabb.push_back(current);
            }
            AABB &back = aabb.back();
            glm::vec2 pos = particles.pos[i];
            if (pos.x < back.l)
                back.l = particles.pos[i].x;
            if (pos.x > back.r)
                back.r = particles.pos[i].x;
            if (pos.y < back.b)
                back.b = particles.pos[i].y;
            if (pos.y > back.t)
                back.t = particles.pos[i].y;
        }
        return aabb;
    }

    inline bool aabbs_intersect(const AABB &a, const AABB &b)
    {
        if (&a == &b)
        {
            printf("случилось невозможное\n");
        }
        return false;
        return !(a.r < b.l || a.l > b.r || a.t < b.b || a.b > b.t);
    }
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabb)
    {
        std::vector<AABBsIntersection> oc;
        for (size_t i = 0; i < aabb.size(); ++i)
            for (size_t j = i + 1; j < aabb.size(); ++j)
                if (aabbs_intersect(aabb[i], aabb[j]))
                    oc.push_back({aabb[i].object_id, aabb[j].object_id});
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
    std::vector<PointEdgeCollisionConstraint> generate_contour_contour_collisions(
        Particles &particles,
        ContourCollider &contourA,
        ContourCollider &contourB)
    {
        std::vector<PointEdgeCollisionConstraint> constraints;

        if (contourA.particle_ids.size() < 2 || contourB.particle_ids.size() < 2)
            return constraints;

        auto detect = [&](const ContourCollider &test, const ContourCollider &target)
        {
            std::vector<glm::vec2> targetPositions;
            for (size_t pid : target.particle_ids)
                targetPositions.push_back(particles.pos[pid]);

            for (size_t pid : test.particle_ids)
            {
                const glm::vec2 &p = particles.pos[pid];

                if (!point_in_polygon(p, targetPositions))
                    continue;

                float minDist = std::numeric_limits<float>::max();
                size_t nearestEdge = 0;
                size_t n = target.particle_ids.size();

                for (size_t i = 0; i < n; ++i)
                {
                    size_t id1 = target.particle_ids[i];
                    size_t id2 = target.particle_ids[(i + 1) % n];
                    const glm::vec2 &e1 = particles.pos[id1];
                    const glm::vec2 &e2 = particles.pos[id2];

                    glm::vec2 edge = e2 - e1;
                    float len = glm::length(edge);
                    if (len < 1e-6f)
                        continue;

                    glm::vec2 dir = edge / len;
                    glm::vec2 toPoint = p - e1;
                    float proj = glm::clamp(glm::dot(toPoint, dir), 0.0f, len);
                    glm::vec2 closest = e1 + dir * proj;
                    float dist = glm::length(p - closest);

                    if (dist < minDist)
                    {
                        minDist = dist;
                        nearestEdge = i;
                    }
                }

                PointEdgeCollisionConstraint constraint;
                constraint.point = pid;
                constraint.edge1 = target.particle_ids[nearestEdge];
                constraint.edge2 = target.particle_ids[(nearestEdge + 1) % target.particle_ids.size()];
                constraint.compliance = (contourA.compliance + contourB.compliance) * 0.5f;
                constraint.staticFriction = (contourA.staticFriction + contourB.staticFriction) * 0.5f;
                constraint.kineticFriction = (contourA.kineticFriction + contourB.kineticFriction) * 0.5f;

                constraints.push_back(constraint);
            }
        };
        detect(contourA, contourB);
        detect(contourB, contourA);
        return constraints;
    }
}
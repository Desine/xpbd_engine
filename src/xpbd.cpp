#include "xpbd.hpp"

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
        p.object_ids.push_back(object_id);
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

    float get_distance_between_particles(Particles &particles, size_t i1, size_t i2)
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

    std::vector<AABB> generate_aabbs(Particles &particles)
    {
        std::vector<AABB> aabb;

        AABB current;
        current.object_id = -1;
        for (size_t i = 0; i < particles.count; i++)
        {
            if (current.object_id != particles.object_ids[i])
            {
                current.object_id = particles.object_ids[i];
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
        if (a.object_id == b.object_id) return false;
        return !(a.r < b.l || a.l > b.r || a.t < b.b || a.b > b.t);
    }
    std::vector<AABBsIntersections> find_aabbs_intersections(std::vector<AABB> aabb)
    {
        std::vector<AABBsIntersections> oc;
        for (size_t i = 0; i < aabb.size(); i++)
            for (size_t j = 0; j < aabb.size(); j++)
                if (aabbs_intersect(aabb[i], aabb[j]))
                    oc.push_back({aabb[i].object_id, aabb[j].object_id});
        return oc;
    }
}
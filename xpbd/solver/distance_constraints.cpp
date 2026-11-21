#include "distance_constraints.hpp"
#include "xpbd/core/particles.hpp"
#include "xpbd/core/constraints.hpp"
#include <vector>

namespace xpbd
{

    void solve_distance_constraints(std::vector<Particle> &p, std::vector<DistanceConstraint> &dc, float dt, bool reversePass)
    {
        if (reversePass)
        {
            for (auto it = dc.rbegin(); it != dc.rend(); it++)
                solve_distance_constraint(p, *it, dt);
        }
        else
        {
            for (auto &c : dc)
                solve_distance_constraint(p, c, dt);
        }
    }
    void solve_distance_constraint(std::vector<Particle> &p, DistanceConstraint &dc, float dt)
    {
        size_t i1 = dc.i1;
        size_t i2 = dc.i2;

        glm::vec2 &p1 = p[i1].pos;
        glm::vec2 &p2 = p[i2].pos;

        float w1 = p[i1].w;
        float w2 = p[i2].w;

        glm::vec2 delta = p1 - p2;
        float len = glm::length(delta);
        if (len < 1e-6f)
            return;

        float C = 0.0f;
        if (len < dc.minDist)
            C = len - dc.minDist;
        else if (len > dc.maxDist)
            C = len - dc.maxDist;
        else
        {
            dc.lambda = 0.0f;
            return;
        }

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
    void apply_distance_constraints_damping(Particles &p, std::vector<DistanceConstraint> &dc, float dt)
    {
    }

    void reset_distance_constraints_lambdas(std::vector<DistanceConstraint> &dc)
    {
        for (auto &c : dc)
            c.lambda = 0;
    }
} // namespace xpbd
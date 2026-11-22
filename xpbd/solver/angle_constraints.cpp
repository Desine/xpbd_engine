#include "angle_constraints.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "../utils.hpp"

#include <vector>

namespace xpbd
{
    void solve_angle_constraints(std::vector<Particle> &p, std::vector<AngleConstraint> &ac, float dt)
    {
        for (auto &c : ac)
            solve_angle_constraint(p, c, dt);
    }

    void solve_angle_constraint(std::vector<Particle> &p, AngleConstraint &ac, float dt)
    {
        size_t i1 = ac.i1, i2 = ac.i2;
        size_t i3 = ac.i3, i4 = ac.i4;

        glm::vec2 dir1 = p[i2].pos - p[i1].pos;
        glm::vec2 dir2 = p[i4].pos - p[i3].pos;

        float len1Sq = glm::dot(dir1, dir1);
        float len2Sq = glm::dot(dir2, dir2);
        if (len1Sq < 1e-6f || len2Sq < 1e-6f)
            return;

        constexpr float TWO_PI = M_PI * 2.0f;
        float angle = Utils::relative_angle_rad(dir1, dir2);
        if (angle < 0.0f)
            angle += TWO_PI;

        float minAngle = std::fmod(ac.minAngle, TWO_PI);
        float maxAngle = std::fmod(ac.maxAngle, TWO_PI);

        if (minAngle < 0.0f)
            minAngle += TWO_PI;
        if (maxAngle < 0.0f)
            maxAngle += TWO_PI;

        bool inside;
        if (minAngle <= maxAngle)
            inside = angle >= minAngle && angle <= maxAngle;
        else
            inside = angle >= minAngle || angle <= maxAngle;
        if (inside)
        {
            ac.lambda = 0.0f;
            return;
        }

        auto signedAngleDiff = [&](float a, float b)
        {
            return std::atan2(std::sin(a - b), std::cos(a - b));
        };

        float diffMin = signedAngleDiff(angle, minAngle);
        float diffMax = signedAngleDiff(angle, maxAngle);

        float C;
        if (std::abs(diffMin) < std::abs(diffMax))
            C = diffMin;
        else
            C = diffMax;

        glm::vec2 grad1 = Utils::perpCCW2d(dir1) / len1Sq;
        glm::vec2 grad2 = Utils::perpCCW2d(dir2) / len2Sq;

        float w1 = p[i1].w, w2 = p[i2].w;
        float w3 = p[i3].w, w4 = p[i4].w;

        float alphaTilde = ac.compliance / (dt * dt);
        float denom =
            (w1 + w2) * glm::dot(grad1, grad1) +
            (w3 + w4) * glm::dot(grad2, grad2) +
            alphaTilde;
        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * ac.lambda) / denom;
        ac.lambda += deltaLambda;

        p[i1].pos += w1 * grad1 * deltaLambda;
        p[i2].pos -= w2 * grad1 * deltaLambda;
        p[i3].pos -= w3 * grad2 * deltaLambda;
        p[i4].pos += w4 * grad2 * deltaLambda;
    }
    void reset_angle_constraints_lambdas(std::vector<AngleConstraint> &ac)
    {
        for (auto &c : ac)
            c.lambda = 0;
    }
} // namespace xpbd
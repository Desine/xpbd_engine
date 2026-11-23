#include "angle_constraints.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "xpbd/core/storage/component_storage.hpp"
#include "xpbd/utils.hpp"

#include <cmath>

namespace xpbd
{
void solve_angle_constraints(
    ParticleStorage &p,
    AngleConstraintStorage &ac,
    float dt)
{
    for (size_t i = 0; i < ac.size(); ++i)
        solve_angle_constraint(p, *ac.at(i), dt);
}

void solve_angle_constraint(
    ParticleStorage &p,
    AngleConstraint &ac,
    float dt)
{
    Particle *particle1 = p.get(ac.id1);
    Particle *particle2 = p.get(ac.id2);
    Particle *particle3 = p.get(ac.id3);
    Particle *particle4 = p.get(ac.id4);

    if (!particle1 || !particle2 ||
        !particle3 || !particle4)
        return;

    glm::vec2 dir1 = particle2->pos - particle1->pos;
    glm::vec2 dir2 = particle4->pos - particle3->pos;

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

    auto signedAngleDiff = [](float a, float b)
    {
        return std::atan2(
            std::sin(a - b),
            std::cos(a - b));
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

    float w1 = particle1->w;
    float w2 = particle2->w;
    float w3 = particle3->w;
    float w4 = particle4->w;

    float alphaTilde = ac.compliance / (dt * dt);

    float denom =
        (w1 + w2) * glm::dot(grad1, grad1) +
        (w3 + w4) * glm::dot(grad2, grad2) +
        alphaTilde;

    if (denom < 1e-6f)
        return;

    float deltaLambda =
        (-C - alphaTilde * ac.lambda) / denom;

    ac.lambda += deltaLambda;

    particle1->pos += w1 * grad1 * deltaLambda;
    particle2->pos -= w2 * grad1 * deltaLambda;
    particle3->pos -= w3 * grad2 * deltaLambda;
    particle4->pos += w4 * grad2 * deltaLambda;
}

void reset_angle_constraints_lambdas(
    AngleConstraintStorage &ac)
{
    for (size_t i = 0; i < ac.size(); ++i)
        ac.at(i)->lambda = 0.0f;
}
} // namespace xpbd

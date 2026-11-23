#include "distance_constraints.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "xpbd/core/storage/component_storage.hpp"

namespace xpbd
{
void solve_distance_constraints(
    ParticleStorage &p,
    DistanceConstraintStorage &dc,
    float dt)
{
    for (size_t i = 0; i < dc.size(); ++i)
        solve_distance_constraint(p, *dc.at(i), dt);
}

void solve_distance_constraint(
    ParticleStorage &p,
    DistanceConstraint &dc,
    float dt)
{
    Particle *particle1 = p.get(dc.id1);
    Particle *particle2 = p.get(dc.id2);

    if (!particle1 || !particle2)
        return;

    glm::vec2 &p1 = particle1->pos;
    glm::vec2 &p2 = particle2->pos;

    float w1 = particle1->w;
    float w2 = particle2->w;

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

    float deltaLambda =
        (-C - alphaTilde * dc.lambda) / denom;

    dc.lambda += deltaLambda;

    glm::vec2 correction = deltaLambda * grad;

    p1 += w1 * correction;
    p2 -= w2 * correction;
}

void apply_distance_constraints_damping(
    ParticleStorage &p,
    DistanceConstraintStorage &dc,
    float dt)
{
}

void reset_distance_constraints_lambdas(
    DistanceConstraintStorage &dc)
{
    for (size_t i = 0; i < dc.size(); ++i)
        dc.at(i)->lambda = 0.0f;
}
} // namespace xpbd

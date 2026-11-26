#include "volume_constraints.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "xpbd/core/storage/component_storage.hpp"
#include "../utils.hpp"

#include <vector>

namespace xpbd
{
float polygon_area(
    const ParticleStorage &storage,
    const std::vector<ParticleId> &ids)
{
    if (ids.size() < 3)
        return 0.0f;

    float area = 0.0f;

    for (size_t i = 0; i < ids.size(); ++i)
    {
        const auto *p1 = storage.get(ids[i]);
        const auto *p2 = storage.get(ids[(i + 1) % ids.size()]);

        if (!p1 || !p2)
            return 0.0f;

        area += Utils::cross(p1->pos, p2->pos);
    }

    return area * 0.5f;
}

void solve_volume_constraints(
    ParticleStorage &p,
    VolumeConstraintStorage &vc,
    float dt)
{
    for (size_t i = 0; i < vc.size(); ++i)
        solve_volume_constraint(p, *vc.at(i), dt);
}

void solve_volume_constraint(
    ParticleStorage &p,
    VolumeConstraint &vc,
    float dt)
{
    const std::vector<ParticleId> &indices = vc.indices;
    size_t N = indices.size();

    float volume = polygon_area(p, indices);
    float C = volume - vc.volume * vc.rest_pressure;

    float denom = 0.0f;
    std::vector<glm::vec2> grads(N);

    for (size_t i = 0; i < N; ++i)
    {
        const Particle *particlePrev =
            p.get(indices[(i + N - 1) % N]);

        const Particle *particleNext =
            p.get(indices[(i + 1) % N]);

        const Particle *particle =
            p.get(indices[i]);

        if (!particlePrev || !particleNext || !particle)
            return;

        const glm::vec2 &pi_prev = particlePrev->pos;
        const glm::vec2 &pi_next = particleNext->pos;

        glm::vec2 grad =
            0.5f * Utils::perpCW2d(pi_next - pi_prev);

        grads[i] = grad;

        float wi = particle->w;
        denom += wi * glm::dot(grad, grad);
    }

    float alphaTilde = vc.compliance / (dt * dt);
    denom += alphaTilde;

    if (denom < 1e-6f)
        return;

    float deltaLambda =
        (-C - alphaTilde * vc.lambda) / denom;

    vc.lambda += deltaLambda;

    for (size_t i = 0; i < N; ++i)
    {
        Particle *particle = p.get(indices[i]);

        if (!particle)
            return;

        float wi = particle->w;

        particle->pos +=
            wi * deltaLambda * grads[i];
    }
}

void reset_volume_constraints_lambdas(
    VolumeConstraintStorage &vc)
{
    for (size_t i = 0; i < vc.size(); ++i)
        vc.at(i)->lambda = 0.0f;
}
} // namespace xpbd

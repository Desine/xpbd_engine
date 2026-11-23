#include "particle.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/storage/component_storage.hpp"

namespace xpbd
{
    void iterate(
        ParticleStorage& p,
        float dt,
        const glm::vec2& gravity
    )
    {
        for (size_t i = 0; i < p.size(); ++i)
        {
            Particle* particle = p.at(i);

            if (particle->w == 0.0f)
                continue;

            particle->vel += gravity * dt;
            particle->prev = particle->pos;
            particle->pos += particle->vel * dt;
        }
    }

    void update_velocities(ParticleStorage& p, float dt)
    {
        for (size_t i = 0; i < p.size(); ++i)
        {
            Particle* particle = p.at(i);
            particle->vel = (particle->pos - particle->prev) / dt;
        }
    }
}

#pragma once
#include "glm/glm.hpp"
#include "object_pool.hpp"
#include "handles.hpp"
#include <vector>

namespace xpbd
{
    struct Particle
    {
        glm::vec2 pos;
        glm::vec2 prev;
        glm::vec2 vel;
        float w;
    };
    struct Particles
    {
        std::vector<glm::vec2> pos;
        std::vector<glm::vec2> prev;
        std::vector<glm::vec2> vel;
        std::vector<float> w;
    };

    struct ParticleStorage
    {
        Particles data;
        ObjectPool<ParticleTag> handles;

        ParticleId add(
            const glm::vec2 &pos,
            float mass,
            const glm::vec2 &vel);

        bool remove(ParticleId id);

        std::optional<uint32_t> get_index(ParticleId id) const;

        bool contains(ParticleId id) const;
    };
} // namespace xpbd
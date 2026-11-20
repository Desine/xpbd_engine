#include "particles.hpp"
#include "object_pool.hpp"
#include "handles.hpp"
#include "glm/glm.hpp"
#include <cstdint>
#include <optional>
#include <vector>

namespace xpbd
{
    ParticleId ParticleStorage::add(
        const glm::vec2 &pos,
        float mass,
        const glm::vec2 &vel)
    {
        uint32_t index = static_cast<uint32_t>(particles.pos.size());

        particles.pos.push_back(pos);
        particles.prev.push_back(pos);
        particles.vel.push_back(vel);
        particles.w.push_back(mass > 0.f ? 1.f / mass : 0.f);

        return handles.create(index);
    }

    bool ParticleStorage::remove(ParticleId id)
    {
        auto indexOpt = handles.get_dense_index(id);

        if (!indexOpt)
            return false;

        uint32_t index = *indexOpt;
        uint32_t last = static_cast<uint32_t>(particles.pos.size() - 1);

        if (index != last)
        {
            auto movedId = handles.get_handle(last);

            particles.pos[index] = particles.pos[last];
            particles.prev[index] = particles.prev[last];
            particles.vel[index] = particles.vel[last];
            particles.w[index] = particles.w[last];

            handles.update_dense_index(*movedId, index);
        }

        particles.pos.pop_back();
        particles.prev.pop_back();
        particles.vel.pop_back();
        particles.w.pop_back();

        return handles.remove(id);
    }

    std::optional<uint32_t> ParticleStorage::get_index(ParticleId id) const
    {
        return handles.get_dense_index(id);
    }

    bool ParticleStorage::contains(ParticleId id) const
    {
        return handles.contains(id);
    }
} // namespace xpbd
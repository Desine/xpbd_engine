#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/storage/component_storage.hpp"

namespace xpbd
{
void iterate(ParticleStorage &p, float dt, const glm::vec2 &gravity);
void update_velocities(ParticleStorage &p, float dt);
} // namespace xpbd
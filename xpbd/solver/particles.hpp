#pragma once
#include "xpbd/core/particles.hpp"

namespace xpbd
{
    void iterate(Particles &p, float dt, const glm::vec2 &gravity);
    void update_velocities(Particles &p, float dt);
} // namespace xpbd
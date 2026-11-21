#pragma once
#include "xpbd/core/particles.hpp"

namespace xpbd
{
    void iterate(std::vector<Particle> &p, float dt, const glm::vec2 &gravity);
    void update_velocities(std::vector<Particle> &p, float dt);
} // namespace xpbd
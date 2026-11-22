#pragma once
#include "glm/glm.hpp"

namespace xpbd
{
    struct Particle
    {
        glm::vec2 pos;
        glm::vec2 prev;
        glm::vec2 vel;
        float w;
    };
} // namespace xpbd
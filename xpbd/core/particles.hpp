#pragma once
#include "glm/glm.hpp"
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
} // namespace xpbd
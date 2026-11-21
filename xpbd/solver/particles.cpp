#include "xpbd/core/particles.hpp"
#include <vector>

namespace xpbd
{
    void iterate(std::vector<Particle> &p, float dt, const glm::vec2 &gravity)
    {
        for (size_t i = 0; i < p.size(); ++i)
        {
            if (p[i].w == 0.0f)
                continue;

            p[i].vel += gravity * dt;
            p[i].prev = p[i].pos;
            p[i].pos += p[i].vel * dt;
        }
    }

    void update_velocities(std::vector<Particle> &p, float dt)
    {
        for (size_t i = 0; i < p.size(); ++i)
            p[i].vel = (p[i].pos - p[i].prev) / dt;
    }

} // namespace xpbd
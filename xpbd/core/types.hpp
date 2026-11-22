#pragma once
#include "glm/glm.hpp"

namespace xpbd
{
    struct AABB
    {
        float l, r, b, t;
        size_t layer, mask;

        bool intersects(const AABB &o) const
        {
            return !(r < o.l || l > o.r || t < o.b || b > o.t);
        }

        AABB get_intersection(const AABB &o) const
        {
            return {
                std::max(l, o.l),
                std::min(r, o.r),
                std::max(b, o.b),
                std::min(t, o.t)};
        }

        bool contains_point(const glm::vec2 &point) const
        {
            return point.x >= l && point.x <= r && point.y >= b && point.y <= t;
        }
    };
} // namespace xpbd
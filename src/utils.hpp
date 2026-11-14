#include "glm/glm.hpp"

namespace Utils
{
    inline float cross(const glm::vec2 &a, const glm::vec2 &b)
    {
        return a.x * b.y - a.y * b.x;
    }
    inline glm::vec2 perpCW2d(const glm::vec2 &v)
    {
        return glm::vec2(v.y, -v.x);
    }
    inline glm::vec2 perpCCW2d(const glm::vec2 &v)
    {
        return glm::vec2(-v.y, v.x);
    }
    inline float relative_angle_rad(glm::vec2 v1, glm::vec2 v2)
    {
        return std::atan2(cross(v1, v2), glm::dot(v1, v2));
    }
    inline float oriented_angle_2d(const glm::vec2 &a, const glm::vec2 &b)
    {
        float angle = glm::acos(glm::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0f, 1.0f));
        return (cross(a, b) < 0.0f) ? -angle : angle;
    }

    inline glm::mat2 rotation_matrix_from_angle_2d(float angle)
    {
        float c = cos(angle);
        float s = sin(angle);
        return glm::mat2(c, -s, s, c);
    }
}
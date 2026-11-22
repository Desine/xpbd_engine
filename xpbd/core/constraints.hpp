#pragma once
#include "glm/glm.hpp"
#include <vector>

namespace xpbd
{
    struct DistanceConstraint
    {
        size_t i1;
        size_t i2;
        float minDist;
        float maxDist;
        float compressionDamping;
        float extensionDamping;
        float compliance;
        float lambda = 0.0f;
    };
    struct AngleConstraint
    {
        size_t i1, i2;
        size_t i3, i4;
        float minAngle;
        float maxAngle;
        float compliance;
        float lambda = 0.0f;
    };
    struct VolumeConstraint
    {
        std::vector<size_t> indices;
        float restVolume;
        float compliance;
        float lambda = 0.0f;
    };
} // namespace xpbd
#pragma once
#include "storage/handles.hpp"

#include "glm/glm.hpp"
#include <vector>

namespace xpbd
{
    struct DistanceConstraint
    {
        ParticleId id1;
        ParticleId id2;
        float min_dist;
        float max_dist;
        float compressionDamping;
        float extensionDamping;
        float compliance;
        float lambda = 0.0f;
    };
    struct AngleConstraint
    {
        ParticleId id1, id2;
        ParticleId id3, id4;
        float min_angle;
        float max_angle;
        float compliance;
        float lambda = 0.0f;
    };
    struct VolumeConstraint
    {
        std::vector<ParticleId> indices;
        float volume;
        float rest_pressure = 1.0f;
        float compliance;
        float lambda = 0.0f;
    };
} // namespace xpbd
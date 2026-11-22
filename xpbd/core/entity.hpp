#pragma once
#include "storage/handles.hpp"

#include <string>
#include <vector>
#include <optional>

namespace xpbd
{
    struct Entity
    {
        std::string name;

        EntityId id;
        EntityId parent;

        std::vector<EntityId> children;

        std::vector<ParticleId> particles;

        std::vector<DistanceConstraintId> distanceConstraints;
        std::vector<AngleConstraintId> angleConstraints;
        std::vector<VolumeConstraintId> volumeConstraints;

        std::vector<PointsColliderId> pointsColliders;
        std::vector<PolygonColliderId> polygonColliders;

        std::vector<AngularMotorId> motors;
        std::vector<AngularServoId> servos;
    };
} // namespace xpbd
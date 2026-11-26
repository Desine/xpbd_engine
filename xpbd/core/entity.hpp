#pragma once
#include "storage/handles.hpp"

#include <string>
#include <vector>
#include <optional>

namespace xpbd
{
    struct Entity
    {
        EntityId parent;

        std::vector<EntityId> children;

        std::vector<ParticleId> particles;

        std::vector<DistanceConstraintId> distance_constraints;
        std::vector<AngleConstraintId> angle_constraints;
        std::vector<VolumeConstraintId> volume_constraints;

        std::vector<PointsColliderId> points_colliders;
        std::vector<PolygonColliderId> polygon_colliders;

        std::vector<AngularMotorId> angular_motors;
        std::vector<AngularServoId> angular_servos;
    };
} // namespace xpbd
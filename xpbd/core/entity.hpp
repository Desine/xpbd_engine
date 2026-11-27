#pragma once
#include "storage/handles.hpp"

#include <string>
#include <vector>
#include <optional>

namespace xpbd
{
    struct Entity
    {
        EntityId parent_id;

        std::vector<EntityId> children_ids;

        std::vector<ParticleId> particles_ids;

        std::vector<DistanceConstraintId> distance_constraint_ids;
        std::vector<AngleConstraintId> angle_constraint_ids;
        std::vector<VolumeConstraintId> volume_constraint_ids;

        std::vector<PointsColliderId> points_collider_ids;
        std::vector<PolygonColliderId> polygon_collider_ids;

        std::vector<AngularMotorId> angular_motor_ids;
        std::vector<AngularServoId> angular_servo_ids;
    };
} // namespace xpbd
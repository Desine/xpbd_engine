#pragma once
#include "glm/glm.hpp"
#include <string>
#include <variant>

#include "types.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"

namespace xpbd
{
    struct init_command
    {
    };

    struct set_tickRate_command
    {
        float tickRate;
    };

    struct add_particle_command
    {
        glm::vec2 pos;
        float mass;
        glm::vec2 vel;
    };

    struct add_distance_constraint_command
    {
        DistanceConstraint constraint;
    };
    struct add_angle_constraint_command
    {
        AngleConstraint constraint;
    };
    struct add_volume_constraint_command
    {
        VolumeConstraint constraint;
    };

    struct add_polygon_collider_command
    {
        PointsCollider collider;
    };

    struct add_points_collider_command
    {
        PointsCollider collider;
    };

    struct spawn_from_json_command
    {
        std::string name;
        glm::vec2 position;
    };

    struct spawn_polygon_command
    {
        glm::vec2 pos;
        float radius;
        size_t segments;
        float mass;
        float compliance;
    };
    struct add_angular_motor_command
    {
        AngularMotor motor;
    };
    struct update_angular_motor_command
    {
        size_t index;
        AngularMotor motor;
    };
    struct add_angular_servo_command
    {
        AngularServo servo;
    };
    struct update_angular_servo_command
    {
        size_t index;
        AngularServo servo;
    };

    using Command = std::variant<
        init_command,
        set_tickRate_command,
        add_particle_command,
        add_distance_constraint_command,
        add_volume_constraint_command,
        add_polygon_collider_command,
        add_points_collider_command,
        spawn_from_json_command,
        spawn_polygon_command,
        add_angular_motor_command,
        update_angular_motor_command,
        add_angular_servo_command,
        update_angular_servo_command>;
} // namespace xpbd
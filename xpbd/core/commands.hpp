#pragma once
#include "actuators.hpp"
#include "storage/handles.hpp"

#include "glm/glm.hpp"
#include <string>
#include <variant>
#include <cstddef>

namespace xpbd
{
struct init_command
{
};
struct set_tickRate_command
{
    float tickRate;
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
struct remove_particle_command
{
    ParticleId id;
};
struct update_angular_motor_command
{
    AngularMotorId id;
    float target_rpm;
    float max_torque;
};
struct update_angular_servo_command
{
    AngularMotorId id;
    float target_angle;
    float kp;
    float kd;
    float max_torque;
};
struct remove_entity_command
{
    EntityId id;
};

using Command = std::variant<
    init_command,
    set_tickRate_command,
    spawn_from_json_command,
    spawn_polygon_command,
    remove_particle_command,
    remove_entity_command,
    update_angular_motor_command,
    update_angular_servo_command>;
} // namespace xpbd

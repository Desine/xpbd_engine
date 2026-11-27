#pragma once
#include "storage/component_storage.hpp"
#include "storage/handles.hpp"

#include <cstdint>
#include <string>

namespace xpbd
{
enum class EntityRemoveMode
{
    Cascade,
    Detach,
    ReparentToParent
};
struct Scene
{

    EntityStorage entity_storage;
    ParticleStorage particle_storage;

    DistanceConstraintStorage distance_constraint_storage;
    AngleConstraintStorage angle_constraint_storage;
    VolumeConstraintStorage volume_constraint_storage;

    PointsColliderStorage points_collider_storage;
    PolygonColliderStorage polygon_collider_storage;

    AngularMotorStorage angular_motor_storage;
    AngularServoStorage angular_servo_storage;

    EntityId add_entity();
    bool remove_entity(EntityId entity_id, EntityRemoveMode mode = EntityRemoveMode::Cascade);

    ParticleId add_particle(const glm::vec2 &pos, float mass, const glm::vec2 &vel = {});
    bool remove_particle(ParticleId id);

    DistanceConstraintId add_distance_constraint(ParticleId id1, ParticleId id2, float compliance);
    DistanceConstraintId add_distance_constraint(ParticleId id1, ParticleId id2, float compliance, float min_dist, float max_dist);
    bool remove_distance_constraint(DistanceConstraintId id);
    AngleConstraintId add_angle_constraint(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float compliance);
    AngleConstraintId add_angle_constraint(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float compliance, float min_angle, float max_angle);
    bool remove_angle_constraint(AngleConstraintId id);
    VolumeConstraintId add_volume_constraint(const std::vector<ParticleId> &indices, float compliance, float rest_pressure = 1.0f);
    VolumeConstraintId add_volume_constraint(const std::vector<ParticleId> &indices, float compliance, float rest_pressure, float volume);
    bool remove_volume_constraint(VolumeConstraintId id);

    PointsColliderId add_points_collider(const std::vector<ParticleId> &indices, float static_friction, float kinetic_friction, float compliance, uint32_t layer, uint32_t mask);
    bool remove_points_collider(PointsColliderId id);
    void remove_particle_from_points_collider(ParticleId particleId);
    PolygonColliderId add_polygon_collider(const std::vector<ParticleId> &indices, float static_friction, float kinetic_friction, float compliance, uint32_t layer, uint32_t mask);
    bool remove_polygon_collider(PolygonColliderId id);

    AngularMotorId add_angular_motor(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float target_rpm, float max_torque);
    bool update_angular_motor(AngularMotorId id, float target_rpm, float max_torque);
    bool remove_angular_motor(AngularMotorId id);

    AngularServoId add_angular_servo(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float target_angle, float kp, float kd, float max_torque);
    bool update_angular_servo(AngularServoId id, float target_angle, float kp, float kd, float max_torque);
    bool remove_angular_servo(AngularServoId id);

    void clear();
    EntityId spawn_entity(
        const std::string &file_name,
        const glm::vec2 &position,
        std::string name = {},
        EntityId root_id = {});

    void print_entity_tree();
};
} // namespace xpbd

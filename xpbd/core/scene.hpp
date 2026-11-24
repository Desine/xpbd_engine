#pragma once
#include "storage/component_storage.hpp"
#include "storage/handles.hpp"

#include <string>

namespace xpbd
{
struct Scene
{
    enum class EntityRemoveMode
    {
        Cascade,
        Detach,
        ReparentToParent
    };

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

    ParticleId add_particle(const glm::vec2 &pos, float mass, const glm::vec2 &vel = {});
    bool remove_particle(ParticleId id);

    DistanceConstraintId add_distance_constraint(ParticleId id1, ParticleId id2, float compliance);
    DistanceConstraintId add_distance_constraint(ParticleId id1, ParticleId id2, float compliance, float minDist, float maxDist);
    bool remove_distance_constraint(DistanceConstraintId id);
    AngleConstraintId add_angle_constraint(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float compliance);
    AngleConstraintId add_angle_constraint(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float compliance, float minAngle, float maxAngle);
    bool remove_angle_constraint(AngleConstraintId id);
    VolumeConstraintId add_volume_constraint(const std::vector<ParticleId> &indices, float compliance, float restPressure = 1.0f);
    VolumeConstraintId add_volume_constraint(const std::vector<ParticleId> &indices, float compliance, float restPressure, float volume);
    bool remove_volume_constraint(VolumeConstraintId id);

    PointsColliderId add_points_collider(const std::vector<ParticleId> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);
    bool remove_points_collider(PointsColliderId id);
    void remove_particle_from_points_collider(ParticleId particleId);
    PolygonColliderId add_polygon_collider(const std::vector<ParticleId> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);
    bool remove_polygon_collider(PolygonColliderId id);

    AngularMotorId add_angular_motor(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float targetRPM, float maxTorque);
    bool update_angular_motor(AngularMotorId id, float targetRPM, float maxTorque);
    bool remove_angular_motor(AngularMotorId id);

    AngularServoId add_angular_servo(ParticleId id1, ParticleId id2, ParticleId id3, ParticleId id4, float targetAngle, float kp, float kd, float maxTorque);
    bool update_angular_servo(AngularServoId id, float targetAngle, float kp, float kd, float maxTorque);
    bool remove_angular_servo(AngularServoId id);

    void load_from_json(
        const std::string &name,
        const glm::vec2 &position);
    void clear();
};
} // namespace xpbd

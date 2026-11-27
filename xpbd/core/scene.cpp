#include "scene.hpp"
#include "nlohmann/json.hpp"
#include "particle.hpp"
#include "storage/component_storage.hpp"
#include "../utils.hpp"
#include "../solver/volume_constraints.hpp"
#include "cluster.hpp"

#include "glm/glm.hpp"
#include "storage/handles.hpp"
#include <cstddef>
#include <cstdint>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <print>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace xpbd
{

namespace
{
template <typename Storage, typename Predicate>
void remove_matching(Storage &storage, Predicate predicate)
{
    for (size_t i = storage.size(); i > 0; --i)
    {
        const auto id = storage.handle_at(i - 1);
        const auto *value = storage.at(i - 1);
        if (id && value && predicate(*value))
            storage.remove(*id);
    }
}
template <typename IdRange>
bool contains_particle(const IdRange &ids, ParticleId particleId)
{
    return std::find(ids.begin(), ids.end(), particleId) != ids.end();
}
} // namespace

void Scene::remove_particle_from_points_collider(ParticleId particleId)
{
    for (size_t i = points_collider_storage.size(); i > 0; --i)
    {
        const size_t index = i - 1;

        auto *collider = points_collider_storage.at(index);
        if (!collider)
            continue;

        auto &indices = collider->indices;

        auto it = std::find(indices.begin(), indices.end(), particleId);
        if (it == indices.end())
            continue;

        indices.erase(it);

        if (indices.empty())
        {
            auto id = points_collider_storage.handle_at(index);
            if (id)
                points_collider_storage.remove(*id);
        }
    }
}

ParticleId Scene::add_particle(
    const glm::vec2 &pos,
    float mass,
    const glm::vec2 &vel)
{
    Particle particle{
        pos,
        pos,
        vel,
        mass > 0.0f ? 1.0f / mass : 0.0f};

    return particle_storage.create(std::move(particle));
}
bool Scene::remove_particle(ParticleId id)
{
    if (!particle_storage.contains(id))
        return false;

    remove_matching(distance_constraint_storage, [id](const DistanceConstraint &c)
                    { return c.id1 == id || c.id2 == id; });
    remove_matching(angle_constraint_storage, [id](const AngleConstraint &c)
                    { return c.id1 == id || c.id2 == id || c.id3 == id || c.id4 == id; });
    remove_matching(volume_constraint_storage, [id](const VolumeConstraint &c)
                    { return contains_particle(c.indices, id); });

    remove_particle_from_points_collider(id);
    remove_matching(polygon_collider_storage, [id](const PointsCollider &c)
                    { return contains_particle(c.indices, id); });

    remove_matching(angular_motor_storage, [id](const AngularMotor &m)
                    { return m.id1 == id || m.id2 == id || m.id3 == id || m.id4 == id; });
    remove_matching(angular_servo_storage, [id](const AngularServo &s)
                    { return s.id1 == id || s.id2 == id || s.id3 == id || s.id4 == id; });

    return particle_storage.remove(id);
}

EntityId Scene::add_entity()
{
    return entity_storage.create();
}
bool Scene::remove_entity(EntityId id, EntityRemoveMode mode)
{
    Entity *entity = entity_storage.get(id);
    if (!entity) return false;

    EntityId parent_id = entity->parent_id;
    Entity *parent_entity = entity_storage.get(parent_id);
    std::vector<EntityId> children_ids = entity->children_ids;

    for (auto &p_id : entity->particles_ids)
        remove_particle(p_id);

    for (auto &dc_id : entity->distance_constraint_ids)
        remove_distance_constraint(dc_id);
    for (auto &ac_id : entity->angle_constraint_ids)
        remove_angle_constraint(ac_id);
    for (auto &vc_id : entity->volume_constraint_ids)
        remove_volume_constraint(vc_id);

    for (auto &pc_id : entity->points_collider_ids)
        remove_points_collider(pc_id);
    for (auto &pc_id : entity->polygon_collider_ids)
        remove_polygon_collider(pc_id);

    for (auto &am_id : entity->angular_motor_ids)
        remove_angular_motor(am_id);
    for (auto &as_id : entity->angular_servo_ids)
        remove_angular_servo(as_id);

    if (mode == EntityRemoveMode::Cascade)
    {
        for (auto &child_id : entity->children_ids)
            remove_entity(child_id, EntityRemoveMode::Cascade);
    }
    else
    {
        for (auto &child_id : children_ids)
        {
            auto child_entity = entity_storage.get(child_id);
            if (!child_entity) continue;

            parent_entity = entity_storage.get(parent_id);
            if (mode == EntityRemoveMode::ReparentToParent && parent_entity)
            {
                parent_entity->children_ids.push_back(child_id);

                child_entity->parent_id = parent_id;
            }
            else if (mode == EntityRemoveMode::Detach)
            {
                child_entity->parent_id = {};
            }
        }
    }

    parent_entity = entity_storage.get(parent_id);
    if (parent_entity)
    {
        auto it = std::find(parent_entity->children_ids.begin(), parent_entity->children_ids.end(), id);
        if (it != parent_entity->children_ids.end())
            parent_entity->children_ids.erase(it);
    }

    return entity_storage.remove(id);
}

DistanceConstraintId Scene::add_distance_constraint(
    ParticleId id1,
    ParticleId id2,
    float compliance)
{
    const auto *p1 = particle_storage.get(id1);
    const auto *p2 = particle_storage.get(id2);
    if (!p1 || !p2)
        return {};

    float distance = glm::distance(p1->pos, p2->pos);

    return add_distance_constraint(
        id1, id2, compliance, distance, distance);
}
DistanceConstraintId Scene::add_distance_constraint(
    ParticleId id1,
    ParticleId id2,
    float compliance,
    float min_dist,
    float max_dist)
{
    if (!particle_storage.contains(id1) ||
        !particle_storage.contains(id2))
        return {};

    DistanceConstraint dc{};
    dc.id1 = id1;
    dc.id2 = id2;
    dc.min_dist = min_dist;
    dc.max_dist = max_dist;
    dc.compliance = compliance;
    return distance_constraint_storage.create(std::move(dc));
}
bool Scene::remove_distance_constraint(DistanceConstraintId id)
{
    return distance_constraint_storage.remove(id);
}

AngleConstraintId Scene::add_angle_constraint(
    ParticleId id1,
    ParticleId id2,
    ParticleId id3,
    ParticleId id4,
    float compliance)
{
    const auto *p1 = particle_storage.get(id1);
    const auto *p2 = particle_storage.get(id2);
    const auto *p3 = particle_storage.get(id3);
    const auto *p4 = particle_storage.get(id4);
    if (!p1 || !p2 || !p3 || !p4)
        return {};

    glm::vec2 a = p2->pos - p1->pos;
    glm::vec2 b = p4->pos - p3->pos;

    float angle = Utils::relative_angle_rad(a, b);

    return add_angle_constraint(
        id1, id2, id3, id4,
        compliance, angle, angle);
}
AngleConstraintId Scene::add_angle_constraint(
    ParticleId id1,
    ParticleId id2,
    ParticleId id3,
    ParticleId id4,
    float compliance,
    float min_angle,
    float max_angle)
{
    if (!particle_storage.contains(id1) ||
        !particle_storage.contains(id2) ||
        !particle_storage.contains(id3) ||
        !particle_storage.contains(id4))
        return {};

    AngleConstraint ac{};
    ac.id1 = id1;
    ac.id2 = id2;
    ac.id3 = id3;
    ac.id4 = id4;
    ac.min_angle = min_angle;
    ac.max_angle = max_angle;
    ac.compliance = compliance;
    return angle_constraint_storage.create(std::move(ac));
}
bool Scene::remove_angle_constraint(AngleConstraintId id)
{
    return angle_constraint_storage.remove(id);
}

VolumeConstraintId Scene::add_volume_constraint(
    const std::vector<ParticleId> &indices,
    float compliance,
    float rest_pressure)
{
    if (indices.size() < 3)
        return {};

    for (ParticleId id : indices)
        if (!particle_storage.contains(id))
            return {};

    float volume = polygon_area(particle_storage, indices);

    return add_volume_constraint(indices, compliance, rest_pressure, volume);
}
VolumeConstraintId Scene::add_volume_constraint(
    const std::vector<ParticleId> &indices,
    float compliance,
    float rest_pressure,
    float volume)
{
    if (indices.size() < 3)
        return {};

    for (ParticleId id : indices)
        if (!particle_storage.contains(id))
            return {};

    VolumeConstraint vc{};
    vc.indices = indices;
    vc.volume = volume;
    vc.rest_pressure = rest_pressure;
    vc.compliance = compliance;
    return volume_constraint_storage.create(std::move(vc));
}
bool Scene::remove_volume_constraint(VolumeConstraintId id)
{
    return volume_constraint_storage.remove(id);
}

PointsColliderId Scene::add_points_collider(
    const std::vector<ParticleId> &indices,
    float static_friction,
    float kinetic_friction,
    float compliance,
    uint32_t layer,
    uint32_t mask)
{
    if (indices.empty())
        return {};

    for (ParticleId id : indices)
        if (!particle_storage.contains(id))
            return {};

    PointsCollider collider{};
    collider.indices = indices;
    collider.static_friction = static_friction;
    collider.kinetic_friction = kinetic_friction;
    collider.compliance = compliance;
    collider.layer = layer;
    collider.mask = mask;
    return points_collider_storage.create(std::move(collider));
}
bool Scene::remove_points_collider(PointsColliderId id)
{
    return points_collider_storage.remove(id);
}

PolygonColliderId Scene::add_polygon_collider(
    const std::vector<ParticleId> &indices,
    float static_friction,
    float kinetic_friction,
    float compliance,
    uint32_t layer,
    uint32_t mask)
{
    if (indices.size() < 3)
        return {};

    for (ParticleId id : indices)
    {
        if (!particle_storage.contains(id))
            return {};
    }

    PointsCollider collider{};

    collider.indices = indices;
    collider.static_friction = static_friction;
    collider.kinetic_friction = kinetic_friction;
    collider.compliance = compliance;
    collider.layer = layer;
    collider.mask = mask;

    return polygon_collider_storage.create(std::move(collider));
}
bool Scene::remove_polygon_collider(PolygonColliderId id)
{
    return polygon_collider_storage.remove(id);
}

AngularMotorId Scene::add_angular_motor(
    ParticleId id1,
    ParticleId id2,
    ParticleId id3,
    ParticleId id4,
    float target_rpm,
    float max_torque)
{
    if (!particle_storage.contains(id1) ||
        !particle_storage.contains(id2) ||
        !particle_storage.contains(id3) ||
        !particle_storage.contains(id4))
        return {};

    AngularMotor motor{};
    motor.id1 = id1;
    motor.id2 = id2;
    motor.id3 = id3;
    motor.id4 = id4;
    motor.target_rpm = target_rpm;
    motor.max_torque = max_torque;
    return angular_motor_storage.create(std::move(motor));
}

bool Scene::update_angular_motor(
    AngularMotorId id,
    float target_rpm,
    float max_torque)
{
    auto *motor = angular_motor_storage.get(id);

    if (!motor)
        return false;

    motor->target_rpm = target_rpm;
    motor->max_torque = max_torque;

    return true;
}
bool Scene::remove_angular_motor(AngularMotorId id)
{
    return angular_motor_storage.remove(id);
}

AngularServoId Scene::add_angular_servo(
    ParticleId id1,
    ParticleId id2,
    ParticleId id3,
    ParticleId id4,
    float target_angle,
    float kp,
    float kd,
    float max_torque)
{
    if (!particle_storage.contains(id1) ||
        !particle_storage.contains(id2) ||
        !particle_storage.contains(id3) ||
        !particle_storage.contains(id4))
        return {};

    AngularServo servo{};
    servo.id1 = id1;
    servo.id2 = id2;
    servo.id3 = id3;
    servo.id4 = id4;
    servo.target_angle = target_angle;
    servo.kp = kp;
    servo.kd = kd;
    servo.max_torque = max_torque;
    return angular_servo_storage.create(std::move(servo));
}
bool Scene::update_angular_servo(
    AngularServoId id,
    float target_angle,
    float kp,
    float kd,
    float max_torque)
{
    auto *servo = angular_servo_storage.get(id);

    if (!servo)
        return false;

    servo->target_angle = target_angle;
    servo->kp = kp;
    servo->kd = kd;
    servo->max_torque = max_torque;

    return true;
}
bool Scene::remove_angular_servo(AngularServoId id)
{
    return angular_servo_storage.remove(id);
}

void Scene::clear()
{
    entity_storage.clear();
    particle_storage.clear();
    distance_constraint_storage.clear();
    angle_constraint_storage.clear();
    volume_constraint_storage.clear();
    points_collider_storage.clear();
    polygon_collider_storage.clear();
    angular_motor_storage.clear();
    angular_servo_storage.clear();
}

void Scene::print_branch(EntityId id, size_t depth)
{
    auto *entity = entity_storage.get(id);
    if (!entity)
        return;

    for (size_t i = 0; i < depth; i++)
        std::print(" - ");
    std::print("{}\n", *entity_storage.get_name(id));

    std::print("   ");
    for (size_t i = 0; i < depth; i++)
        std::print(" - ");
    std::print("distance_constraints.size(): {}\n", entity->distance_constraint_ids.size());

    for (EntityId child_id : entity->children_ids)
    {
        print_branch(child_id, depth + 1);
    }
}
void Scene::print_entity_tree()
{
    std::println("\n--- --- Scene::print_entity_tree() --- ---");

    for (size_t i = 0; i < entity_storage.size(); i++)
    {
        Entity *entity = entity_storage.at(i);
        if (entity && !entity->parent_id.valid())
        {
            if (auto id = entity_storage.handle_at(i))
                print_branch(*id);
        }
    }
}

} // namespace xpbd

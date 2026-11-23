#include "scene.hpp"
#include "storage/component_storage.hpp"
#include "../utils.hpp"
#include "../solver/volume_constraints.hpp"
#include "../json_body_loader.hpp"
#include "cluster.hpp"

#include "glm/glm.hpp"
#include <memory>
#include <stdexcept>

namespace xpbd
{
EntityId Scene::add_entity()
{
    EntityId id = entity_storage.create();

    auto *entity = entity_storage.get(id);

    if (!entity)
        return {};

    entity->id = id;

    return id;
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
    float minDist,
    float maxDist)
{
    if (!particle_storage.contains(id1) ||
        !particle_storage.contains(id2))
        return {};

    DistanceConstraint dc{};

    dc.id1 = id1;
    dc.id2 = id2;
    dc.minDist = minDist;
    dc.maxDist = maxDist;
    dc.compliance = compliance;

    return distance_constraint_storage.create(std::move(dc));
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
    float minAngle,
    float maxAngle)
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
    ac.minAngle = minAngle;
    ac.maxAngle = maxAngle;
    ac.compliance = compliance;

    return angle_constraint_storage.create(std::move(ac));
}

VolumeConstraintId Scene::add_volume_constraint(
    const std::vector<ParticleId> &indices,
    float compliance)
{
    if (indices.size() < 3)
        return {};

    for (ParticleId id : indices)
    {
        if (!particle_storage.contains(id))
            return {};
    }

    VolumeConstraint vc{};

    vc.indices = indices;
    vc.volume = polygon_area(particle_storage, indices);
    vc.compliance = compliance;

    return volume_constraint_storage.create(std::move(vc));
}

VolumeConstraintId Scene::add_volume_constraint(
    const std::vector<ParticleId> &indices,
    float compliance,
    float restPressure)
{
    if (indices.size() < 3)
        return {};

    for (ParticleId id : indices)
    {
        if (!particle_storage.contains(id))
            return {};
    }

    VolumeConstraint vc{};

    vc.indices = indices;
    vc.volume = polygon_area(particle_storage, indices);
    vc.restPressure = restPressure;
    vc.compliance = compliance;

    return volume_constraint_storage.create(std::move(vc));
}

VolumeConstraintId Scene::add_volume_constraint(
    const std::vector<ParticleId> &indices,
    float compliance,
    float restPressure,
    float volume)
{
    if (indices.size() < 3)
        return {};

    for (ParticleId id : indices)
    {
        if (!particle_storage.contains(id))
            return {};
    }

    VolumeConstraint vc{};

    vc.indices = indices;
    vc.volume = volume;
    vc.restPressure = restPressure;
    vc.compliance = compliance;

    return volume_constraint_storage.create(std::move(vc));
}

PointsColliderId Scene::add_points_collider(
    const std::vector<ParticleId> &indices,
    float staticFriction,
    float kineticFriction,
    float compliance,
    int8_t layer,
    int8_t mask)
{
    if (indices.empty())
        return {};

    for (ParticleId id : indices)
    {
        if (!particle_storage.contains(id))
            return {};
    }

    PointsCollider collider{};

    collider.indices = indices;
    collider.staticFriction = staticFriction;
    collider.kineticFriction = kineticFriction;
    collider.compliance = compliance;
    collider.layer = layer;
    collider.mask = mask;

    return points_collider_storage.create(std::move(collider));
}

PolygonColliderId Scene::add_polygon_collider(
    const std::vector<ParticleId> &indices,
    float staticFriction,
    float kineticFriction,
    float compliance,
    int8_t layer,
    int8_t mask)
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
    collider.staticFriction = staticFriction;
    collider.kineticFriction = kineticFriction;
    collider.compliance = compliance;
    collider.layer = layer;
    collider.mask = mask;

    return polygon_collider_storage.create(std::move(collider));
}

AngularMotorId Scene::add_angular_motor(
    ParticleId id1,
    ParticleId id2,
    ParticleId id3,
    ParticleId id4,
    float targetRPM,
    float maxTorque)
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
    motor.targetRPM = targetRPM;
    motor.maxTorque = maxTorque;

    return angular_motor_storage.create(std::move(motor));
}

bool Scene::update_angular_motor(
    AngularMotorId id,
    float targetRPM,
    float maxTorque)
{
    auto *motor = angular_motor_storage.get(id);

    if (!motor)
        return false;

    motor->targetRPM = targetRPM;
    motor->maxTorque = maxTorque;

    return true;
}

AngularServoId Scene::add_angular_servo(
    ParticleId id1,
    ParticleId id2,
    ParticleId id3,
    ParticleId id4,
    float targetAngle,
    float kp,
    float kd,
    float maxTorque)
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
    servo.targetAngle = targetAngle;
    servo.kp = kp;
    servo.kd = kd;
    servo.maxTorque = maxTorque;

    return angular_servo_storage.create(std::move(servo));
}

bool Scene::update_angular_servo(
    AngularServoId id,
    float targetAngle,
    float kp,
    float kd,
    float maxTorque)
{
    auto *servo = angular_servo_storage.get(id);

    if (!servo)
        return false;

    servo->targetAngle = targetAngle;
    servo->kp = kp;
    servo->kd = kd;
    servo->maxTorque = maxTorque;

    return true;
}

void Scene::load_from_json(
    const std::string &name,
    const glm::vec2 &position)
{
    std::unique_ptr<cluster::Cluster> cluster(json_body_loader::load(name, position));
    if (!cluster)
        return;

    std::vector<ParticleId> particles;
    particles.reserve(cluster->particles.size());
    for (const auto &particle : cluster->particles)
    {
        Particle value{particle.pos, particle.prev, particle.vel, particle.w};
        particles.push_back(particle_storage.create(std::move(value)));
    }

    auto particle_id = [&particles](size_t index) -> ParticleId
    {
        if (index >= particles.size())
            throw std::out_of_range("Body contains an invalid particle index");
        return particles[index];
    };

    for (const auto &dc : cluster->distanceConstraints)
        add_distance_constraint(particle_id(dc.i1), particle_id(dc.i2),
                                dc.compliance, dc.minDist, dc.maxDist);

    for (const auto &ac : cluster->angleConstraints)
        add_angle_constraint(particle_id(ac.i1), particle_id(ac.i2),
                             particle_id(ac.i3), particle_id(ac.i4),
                             ac.compliance, ac.minAngle, ac.maxAngle);

    for (const auto &vc : cluster->volumeConstraints)
    {
        std::vector<ParticleId> ids;
        ids.reserve(vc.indices.size());
        for (size_t index : vc.indices)
            ids.push_back(particle_id(index));
        add_volume_constraint(ids, vc.compliance, vc.restPressure, vc.volume);
    }

    for (const auto &collider : cluster->polygonColliders)
    {
        std::vector<ParticleId> ids;
        ids.reserve(collider.indices.size());
        for (size_t index : collider.indices)
            ids.push_back(particle_id(index));
        add_polygon_collider(ids, collider.staticFriction, collider.kineticFriction,
                             collider.compliance, collider.layer, collider.mask);
    }
    for (const auto &collider : cluster->pointsColliders)
    {
        std::vector<ParticleId> ids;
        ids.reserve(collider.indices.size());
        for (size_t index : collider.indices)
            ids.push_back(particle_id(index));
        add_points_collider(ids, collider.staticFriction, collider.kineticFriction,
                            collider.compliance, collider.layer, collider.mask);
    }

    for (const auto &motor : cluster->angularMotors)
        add_angular_motor(particle_id(motor.i1), particle_id(motor.i2),
                          particle_id(motor.i3), particle_id(motor.i4),
                          motor.targetRPM, motor.maxTorque);
    for (const auto &servo : cluster->angularServos)
        add_angular_servo(particle_id(servo.i1), particle_id(servo.i2),
                          particle_id(servo.i3), particle_id(servo.i4),
                          servo.targetAngle, servo.kp, servo.kd, servo.maxTorque);
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

} // namespace xpbd

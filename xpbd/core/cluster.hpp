#pragma once
#include "entity.hpp"
#include "particle.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"

#include <cstddef>
#include <vector>
#include <string>

namespace xpbd
{
namespace cluster
{
struct Particle
{
    glm::vec2 pos;
    glm::vec2 prev;
    glm::vec2 vel;
    float w;
};
struct DistanceConstraint
{
    size_t i1;
    size_t i2;
    float minDist;
    float maxDist;
    float compressionDamping;
    float extensionDamping;
    float compliance;
    float lambda = 0.0f;
};
struct AngleConstraint
{
    size_t i1, i2;
    size_t i3, i4;
    float minAngle;
    float maxAngle;
    float compliance;
    float lambda = 0.0f;
};
struct VolumeConstraint
{
    std::vector<size_t> indices;
    float volume;
    float restPressure = 1.0f;
    float compliance;
    float lambda = 0.0f;
};
struct PointsCollider
{
    std::vector<size_t> indices;
    float staticFriction;
    float kineticFriction;
    float compliance;
    int8_t layer;
    int8_t mask;
};
struct PointEdgeCollision
{
    size_t point;
    size_t edge1;
    size_t edge2;
};
struct PolygonCache
{
    std::vector<glm::vec2> positions;
    size_t positions_size;
    std::vector<glm::vec2> edges;
    std::vector<float> edgesLen2;
};
struct PointsPolygonCollision
{
    std::vector<size_t> points;
    std::vector<size_t> polygon;
    float staticFriction;
    float kineticFriction;
    float compliance;
    AABB box;
    PolygonCache cache;
};
struct PointEdgeCollisionConstraints
{
    size_t point;
    size_t edge1;
    size_t edge2;
    float staticFriction;
    float kineticFriction;
    float compliance;
    float lambdaN = 0.0f;
    float lambdaT = 0.0f;
};
struct AngularMotor
{
    size_t i1, i2;
    size_t i3, i4;
    float targetRPM;
    float maxTorque;
};
struct AngularServo
{
    size_t i1, i2;
    size_t i3, i4;
    float targetAngle;
    float kp;
    float kd;
    float maxTorque;
};

struct Entity
{
    std::string name;

    size_t parent = SIZE_MAX;
    std::vector<size_t> children;

    std::vector<size_t> particles;

    std::vector<size_t> distanceConstraints;
    std::vector<size_t> angleConstraints;
    std::vector<size_t> volumeConstraints;

    std::vector<size_t> pointsColliders;
    std::vector<size_t> polygonColliders;

    std::vector<size_t> motors;
    std::vector<size_t> servos;
};

struct Cluster
{
    std::vector<Entity> entities;

    std::vector<std::pair<std::string, size_t>> particleIds;
    std::vector<Particle> particles;

    std::vector<DistanceConstraint> distanceConstraints;
    std::vector<AngleConstraint> angleConstraints;
    std::vector<VolumeConstraint> volumeConstraints;

    std::vector<PointsCollider> polygonColliders;
    std::vector<PointsCollider> pointsColliders;

    std::vector<AngularMotor> angularMotors;
    std::vector<AngularServo> angularServos;
};

void merge_clusters(Cluster &into, Cluster &from);

void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass, glm::vec2 vel);

void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance);
void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float minDist, float maxDist);
void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance);
void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float minAngle, float maxAngle);
void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance);
void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float restPressure);

void add_polygon_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);
void add_points_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);

void add_angular_motor(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float targetRPM, float maxTorque);
void add_angular_servo(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float targetAngle, float kp, float kd, float maxTorque);

void spawn_polygon(Cluster &cluster, glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
} // namespace cluster
} // namespace xpbd
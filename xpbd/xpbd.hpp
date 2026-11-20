#pragma once
#include "glm/glm.hpp"
#include <vector>

#include "core/types.hpp"
#include "core/particles.hpp"
#include "core/constraints.hpp"
#include "core/colliders.hpp"
#include "core/actuators.hpp"
#include "core/cluster.hpp"
#include "core/commands.hpp"
#include "core/world.hpp"

#include "solver/particles.hpp"
#include "solver/distance_constraints.hpp"
#include "solver/angle_constraints.hpp"
#include "solver/volume_constraints.hpp"
#include "solver/actuators.hpp"
#include "solver/collision.hpp"
#include "solver/spatial_hashing.hpp"

namespace xpbd{
    void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass, glm::vec2 vel);

    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance);
    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float minDist, float maxDist);
    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance);
    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float minAngle, float maxAngle);
    void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance);
    void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float restPressure);

    void add_polygon_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);
    void add_points_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);

    void add_angular_motor(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float targetVelocity, float maxTorque);
    void update_angular_motor(Cluster &cluster, size_t index, float targetVelocity, float maxTorque);
    void add_angular_servo(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float targetAngle, float kp, float kd, float maxTorque);
    void update_angular_servo(Cluster &cluster, size_t index, float targetAngle, float kp, float kd, float maxTorque);

    void spawn_from_json(Cluster &cluster, const std::string &name, const glm::vec2 &position);
    void spawn_polygon(Cluster &cluster, glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
}
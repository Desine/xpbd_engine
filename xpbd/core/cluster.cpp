#include "cluster.hpp"
#include "entity.hpp"
#include "particle.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"
#include "../utils.hpp"

#include <vector>
#include <string>

namespace xpbd
{
namespace cluster
{
void merge_clusters(Cluster &into, Cluster &from)
{
    if (&into == &from)
        return;

    const size_t offset = into.particles.size();

    // Particles
    into.particles.insert(into.particles.end(),
                          from.particles.begin(), from.particles.end());

    // Particle IDs
    into.particleIds.reserve(
        into.particleIds.size() + from.particleIds.size());

    for (auto &entry : from.particleIds)
    {
        entry.second += offset;
        into.particleIds.push_back(std::move(entry));
    }

    // Distance constraints
    into.distance_constraints.reserve(into.distance_constraints.size() + from.distance_constraints.size());
    for (auto &dc : from.distance_constraints)
    {
        dc.i1 += offset;
        dc.i2 += offset;
        into.distance_constraints.push_back(dc);
    }

    // Angle constraints
    into.angleConstraints.reserve(into.angleConstraints.size() + from.angleConstraints.size());
    for (auto &ac : from.angleConstraints)
    {
        ac.i1 += offset;
        ac.i2 += offset;
        ac.i3 += offset;
        ac.i4 += offset;
        into.angleConstraints.push_back(ac);
    }

    // Volume constraints
    into.volumeConstraints.reserve(into.volumeConstraints.size() + from.volumeConstraints.size());
    for (auto &vc : from.volumeConstraints)
    {
        for (auto &idx : vc.indices)
            idx += offset;
        into.volumeConstraints.push_back(vc);
    }

    // Polygon colliders
    into.polygonColliders.reserve(into.polygonColliders.size() + from.polygonColliders.size());
    for (auto &pc : from.polygonColliders)
    {
        for (auto &idx : pc.indices)
            idx += offset;
        into.polygonColliders.push_back(pc);
    }

    // Points colliders
    into.pointsColliders.reserve(into.pointsColliders.size() + from.pointsColliders.size());
    for (auto &pc : from.pointsColliders)
    {
        for (auto &idx : pc.indices)
            idx += offset;
        into.pointsColliders.push_back(pc);
    }

    // Angular motors
    into.angularMotors.reserve(into.angularMotors.size() + from.angularMotors.size());
    for (auto &am : from.angularMotors)
    {
        am.i1 += offset;
        am.i2 += offset;
        am.i3 += offset;
        am.i4 += offset;
        into.angularMotors.push_back(am);
    }

    // Angular servos
    into.angularServos.reserve(into.angularServos.size() + from.angularServos.size());
    for (auto &as : from.angularServos)
    {
        as.i1 += offset;
        as.i2 += offset;
        as.i3 += offset;
        as.i4 += offset;
        into.angularServos.push_back(as);
    }

    from = Cluster();
}

void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass, glm::vec2 vel = {0, 0})
{
    cluster.particles.emplace_back(pos, pos, vel, mass == 0 ? 0 : 1 / mass);
}

void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float min_dist, float max_dist)
{
    float compressionDamping = 0.5f;
    float extensionDamping = 0.5f;
    cluster.distance_constraints.emplace_back(i1, i2, min_dist, max_dist, compressionDamping, extensionDamping, compliance);
}
void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance)
{
    float restDist = glm::distance(cluster.particles[i1].pos, cluster.particles[i2].pos);
    add_distance_constraint(cluster, i1, i2, compliance, restDist, restDist);
}

void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float min_angle, float max_angle)
{
    cluster.angleConstraints.emplace_back(i1, i2, i3, i4, min_angle, max_angle, compliance);
}
void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance)
{
    glm::vec2 v1 = cluster.particles[i1].pos;
    glm::vec2 v2 = cluster.particles[i2].pos;
    glm::vec2 v3 = cluster.particles[i3].pos;
    glm::vec2 v4 = cluster.particles[i4].pos;

    float restAngle = Utils::relative_angle_rad(v2 - v1, v4 - v3);
    add_angle_constraint(cluster, i1, i2, i3, i4, compliance, restAngle, restAngle);
}

float polygon_area(const std::vector<Particle> &p, const std::vector<size_t> &indices)
{
    float area = 0.0f;
    size_t n = indices.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++)
    {
        const glm::vec2 &p0 = p[indices[j]].pos;
        const glm::vec2 &p1 = p[indices[i]].pos;
        area += Utils::cross(p0, p1);
    }
    return 0.5f * area;
}
void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float rest_pressure)
{
    float volume = polygon_area(cluster.particles, indices);
    cluster.volumeConstraints.emplace_back(indices, volume, rest_pressure, compliance);
}
void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance)
{
    float rest_pressure = 1.0f;
    add_volume_constraint(cluster, indices, compliance, rest_pressure);
}

void add_polygon_collider(Cluster &cluster, const std::vector<size_t> &indices, float static_friction, float kinetic_friction, float compliance, int8_t layer, int8_t mask)
{
    cluster.polygonColliders.emplace_back(indices, static_friction, kinetic_friction, compliance, layer, mask);
}
void add_points_collider(Cluster &cluster, const std::vector<size_t> &indices, float static_friction, float kinetic_friction, float compliance, int8_t layer, int8_t mask)
{
    cluster.pointsColliders.emplace_back(indices, static_friction, kinetic_friction, compliance, layer, mask);
}
void add_angular_motor(
    Cluster &cluster,
    size_t i1, size_t i2,
    size_t i3, size_t i4,
    float target_rpm,
    float max_torque)
{
    cluster.angularMotors.emplace_back(i1, i2, i3, i4, target_rpm, max_torque);
}
void add_angular_servo(
    Cluster &cluster,
    size_t i1, size_t i2,
    size_t i3, size_t i4,
    float target_angle,
    float kp,
    float kd,
    float max_torque)
{
    cluster.angularServos.emplace_back(i1, i2, i3, i4, target_angle, kp, kd, max_torque);
}

void spawn_polygon(Cluster &cluster, glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
{
    if (segments < 3)
        segments = 3;

    size_t start = cluster.particles.size();
    std::vector<size_t> ids;

    float angleStep = 2.0f * float(M_PI) / segments;
    for (size_t i = 0; i < segments; ++i)
    {
        float angle = i * angleStep;
        glm::vec2 dir = glm::vec2(cosf(angle), sinf(angle));
        add_particle(cluster, dir * radius + pos, mass / segments);

        ids.emplace_back(start + i);
    }

    for (size_t i = 0; i < segments; ++i)
    {
        add_distance_constraint(cluster, start + i, start + (i + 1) % segments, compliance);
        // add_angle_constraint(cluster, start + i, start + (i + 1) % segments, start + (i + 1) % segments, start + (i + 2) % segments, compliance);
    }

    add_volume_constraint(cluster, ids, compliance);

    float static_friction = 0.4f;
    float kinetic_friction = 0.3f;

    int8_t layer = 1;
    int8_t mask = 1;

    add_polygon_collider(cluster, ids, static_friction, kinetic_friction, 0.0f, layer, mask);
}
} // namespace cluster
} // namespace xpbd
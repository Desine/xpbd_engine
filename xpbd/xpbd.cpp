#include "xpbd.hpp"

#include "utils.hpp"

#include "tracy/Tracy.hpp"
#include "json_body_loader.hpp"
namespace xpbd
{
    void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass, glm::vec2 vel = {0, 0})
    {
        cluster.particles.emplace_back(pos, pos, vel, mass == 0 ? 0 : 1 / mass);
    }

    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float minDist, float maxDist)
    {
        float compressionDamping = 0.5f;
        float extensionDamping = 0.5f;
        cluster.distanceConstraints.emplace_back(i1, i2, minDist, maxDist, compressionDamping, extensionDamping, compliance);
    }
    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance)
    {
        float restDist = glm::distance(cluster.particles[i1].pos, cluster.particles[i2].pos);
        add_distance_constraint(cluster, i1, i2, compliance, restDist, restDist);
    }

    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float minAngle, float maxAngle)
    {
        cluster.angleConstraints.emplace_back(i1, i2, i3, i4, minAngle, maxAngle, compliance);
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

    void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float restVolume)
    {
        cluster.volumeConstraints.emplace_back(indices, restVolume, compliance);
    }
    void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance)
    {
        float restVolume = compute_polygon_area(cluster.particles, indices);
        add_volume_constraint(cluster, indices, compliance, restVolume);
    }

    void add_polygon_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask)
    {
        cluster.polygonColliders.emplace_back(indices, staticFriction, kineticFriction, compliance, layer, mask);
    }
    void add_points_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask)
    {
        cluster.pointsColliders.emplace_back(indices, staticFriction, kineticFriction, compliance, layer, mask);
    }
    void add_angular_motor(
        Cluster &cluster,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        float targetVelocity,
        float maxTorque)
    {
        cluster.angularMotors.emplace_back(i1, i2, i3, i4, targetVelocity, maxTorque);
    }
    void update_angular_motor(
        Cluster &cluster,
        size_t index,
        float targetVelocity,
        float maxTorque)
    {
        if (index >= cluster.angularMotors.size())
            return;

        cluster.angularMotors[index].targetRPM = targetVelocity;
        cluster.angularMotors[index].maxTorque = maxTorque;
    }
    void add_angular_servo(
        Cluster &cluster,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        float targetAngle,
        float kp,
        float kd,
        float maxTorque)
    {
        cluster.angularServos.emplace_back(i1, i2, i3, i4, targetAngle, kp, kd, maxTorque);
    }
    void update_angular_servo(
        Cluster &cluster,
        size_t index,
        float targetAngle,
        float kp,
        float kd,
        float maxTorque)
    {
        if (index >= cluster.angularServos.size())
            return;

        cluster.angularServos[index].targetAngle = targetAngle;
        cluster.angularServos[index].kp = kp;
        cluster.angularServos[index].kd = kd;
        cluster.angularServos[index].maxTorque = maxTorque;
    }

    void spawn_from_json(Cluster &cluster, const std::string &name, const glm::vec2 &position)
    {
        json_body_loader::load(cluster, name, position);
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

        float staticFriction = 0.4f;
        float kineticFriction = 0.3f;

        int8_t layer = 1;
        int8_t mask = 1;

        add_polygon_collider(cluster, ids, staticFriction, kineticFriction, 0.0f, layer, mask);
    }
} // namespace xpbd
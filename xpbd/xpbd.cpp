#include "xpbd.hpp"

#include "utils.hpp"

#include "tracy/Tracy.hpp"
#include "json_body_loader.hpp"
namespace xpbd
{
    
    void merge_clusters(Cluster &into, Cluster &from)
    {
        if (&into == &from)
            return;

        const size_t offset = into.particles.pos.size();

        // Particles
        into.particles.pos.insert(into.particles.pos.end(),
                                  from.particles.pos.begin(), from.particles.pos.end());
        into.particles.prev.insert(into.particles.prev.end(),
                                   from.particles.prev.begin(), from.particles.prev.end());
        into.particles.vel.insert(into.particles.vel.end(),
                                  from.particles.vel.begin(), from.particles.vel.end());
        into.particles.w.insert(into.particles.w.end(),
                                from.particles.w.begin(), from.particles.w.end());

        // Distance constraints
        into.distanceConstraints.reserve(into.distanceConstraints.size() + from.distanceConstraints.size());
        for (auto &dc : from.distanceConstraints)
        {
            dc.i1 += offset;
            dc.i2 += offset;
            into.distanceConstraints.push_back(dc);
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
        cluster.particles.pos.emplace_back(pos);
        cluster.particles.prev.emplace_back(pos);
        cluster.particles.vel.emplace_back(vel);
        cluster.particles.w.emplace_back(mass == 0 ? 0 : 1 / mass);
    }

    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float minDist, float maxDist)
    {
        float compressionDamping = 0.5f;
        float extensionDamping = 0.5f;
        cluster.distanceConstraints.emplace_back(i1, i2, minDist, maxDist, compressionDamping, extensionDamping, compliance);
    }
    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance)
    {
        float restDist = glm::distance(cluster.particles.pos[i1], cluster.particles.pos[i2]);
        add_distance_constraint(cluster, i1, i2, compliance, restDist, restDist);
    }

    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float minAngle, float maxAngle)
    {
        cluster.angleConstraints.emplace_back(i1, i2, i3, i4, minAngle, maxAngle, compliance);
    }
    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance)
    {
        glm::vec2 v1 = cluster.particles.pos[i1];
        glm::vec2 v2 = cluster.particles.pos[i2];
        glm::vec2 v3 = cluster.particles.pos[i3];
        glm::vec2 v4 = cluster.particles.pos[i4];

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

        size_t start = cluster.particles.pos.size();
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

#include "json_body_loader.hpp"
#include "utils.hpp"
#include "core/cluster.hpp"

#include "nlohmann/json.hpp"

#include <algorithm>
#include <fstream>
#include <memory>
#include <optional>
#include <print>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

std::string json_body_loader::folder = "";
const std::string json_body_loader::extension = ".body";

inline glm::vec2 get_vec2(
    const nlohmann::json &p,
    const char *key,
    glm::vec2 def = {0, 0})
{
    auto j = p.value(key, nlohmann::json::array({def.x, def.y}));
    return {j[0].get<float>(), j[1].get<float>()};
}

std::optional<size_t> get_index_by_id(
    const xpbd::cluster::Cluster &cluster,
    const std::string &id)
{
    auto result = std::find_if(
        cluster.particleIds.begin(),
        cluster.particleIds.end(),
        [&id](const auto &p)
        {
            return p.first == id;
        });

    if (result != cluster.particleIds.end())
        return result->second;

    return std::nullopt;
}

size_t get_index(
    const xpbd::cluster::Cluster &cluster,
    const nlohmann::json &element,
    size_t base_index)
{
    size_t index;

    if (element.is_string())
    {
        const std::string id = element.get<std::string>();
        auto idFound = get_index_by_id(cluster, id);

        if (!idFound)
            throw std::runtime_error(
                "Unknown particle ID: " + id);

        index = *idFound;
    }
    else
    {
        const size_t localIndex = element.get<size_t>();
        index = localIndex + base_index;
    }

    if (index >= cluster.particles.size())
        throw std::runtime_error(
            "Particle index is out of range: " + std::to_string(index));

    return index;
}

float polygon_area(
    const std::vector<xpbd::cluster::Particle> &p,
    const std::vector<size_t> &indices)
{
    if (indices.size() < 3)
        return 0.0f;

    float area = 0.0f;
    const size_t n = indices.size();

    for (size_t i = 0, j = n - 1; i < n; j = i++)
    {
        const glm::vec2 &p0 = p[indices[j]].pos;
        const glm::vec2 &p1 = p[indices[i]].pos;

        area += Utils::cross(p0, p1);
    }

    return 0.5f * area;
}

xpbd::cluster::Cluster *json_body_loader::load(
    const std::string &filename,
    glm::vec2 offset)
{
    std::ifstream file(folder + filename + extension);

    if (!file.is_open())
    {
        std::println(
            "Failed to open file: {}|{}|{}",
            folder,
            filename,
            extension);

        return nullptr;
    }

    auto *cluster = new xpbd::cluster::Cluster();

    try
    {
        nlohmann::json j;
        file >> j;

        // Load required bodies.
        if (j.contains("Require") && j["Require"].is_array())
        {
            for (const auto &r : j["Require"])
            {
                const std::string requiredName = r["name"].get<std::string>();
                const glm::vec2 requiredOffset = get_vec2(r, "offset") + offset;

                std::unique_ptr<xpbd::cluster::Cluster> required(
                    json_body_loader::load(
                        requiredName,
                        requiredOffset));

                if (!required)
                {
                    throw std::runtime_error(
                        "Failed to load required body: " + requiredName);
                }

                xpbd::cluster::merge_clusters(*cluster, *required);
            }
        }

        // Offset for indices belonging to this body.
        const size_t base_index = cluster->particles.size();

        // Load Particles.
        if (j.contains("Particles") && j["Particles"].is_array())
        {
            const glm::vec2 j_offset = get_vec2(j, "offset");
            const glm::vec2 j_velocity = get_vec2(j, "velocity");

            size_t index = 0;

            for (const auto &p : j["Particles"])
            {
                const glm::vec2 pos = get_vec2(p, "pos") + offset + j_offset;
                const glm::vec2 vel = get_vec2(p, "vel") + j_velocity;
                const float mass = p.value("mass", 0.0f);

                xpbd::cluster::add_particle(
                    *cluster,
                    pos,
                    mass,
                    vel);

                if (p.contains("id"))
                {
                    cluster->particleIds.emplace_back(
                        p["id"].get<std::string>(),
                        index + base_index);
                }

                index++;
            }
        }

        // Load Distance Constraints.
        if (j.contains("DistanceConstraints") &&
            j["DistanceConstraints"].is_array())
        {
            for (const auto &dc : j["DistanceConstraints"])
            {
                const size_t i1 = get_index(*cluster, dc["i1"], base_index);
                const size_t i2 = get_index(*cluster, dc["i2"], base_index);

                const float compliance = dc.value("compliance", 0.0f);

                if (dc.contains("minDist") && dc.contains("maxDist"))
                {
                    const float minDist = dc.value("minDist", 0.0f);
                    const float maxDist = dc.value("maxDist", 0.0f);

                    xpbd::cluster::add_distance_constraint(
                        *cluster,
                        i1,
                        i2,
                        compliance,
                        minDist,
                        maxDist);
                }
                else if (dc.contains("restDist"))
                {
                    const float restDist = dc.value("restDist", 0.0f);

                    xpbd::cluster::add_distance_constraint(
                        *cluster,
                        i1,
                        i2,
                        compliance,
                        restDist,
                        restDist);
                }
                else
                {
                    xpbd::cluster::add_distance_constraint(
                        *cluster,
                        i1,
                        i2,
                        compliance);
                }
            }
        }

        // Load Angle Constraints.
        if (j.contains("AngleConstraints") &&
            j["AngleConstraints"].is_array())
        {
            for (const auto &ac : j["AngleConstraints"])
            {
                const size_t i1 = get_index(*cluster, ac["i1"], base_index);
                const size_t i2 = get_index(*cluster, ac["i2"], base_index);
                const size_t i3 = get_index(*cluster, ac["i3"], base_index);
                const size_t i4 = get_index(*cluster, ac["i4"], base_index);

                const float compliance = ac.value("compliance", 0.0f);

                if (ac.contains("minAngle") && ac.contains("maxAngle"))
                {
                    const float minAngle = glm::radians(ac.value("minAngle", 0.0f));
                    const float maxAngle = glm::radians(ac.value("maxAngle", 0.0f));

                    xpbd::cluster::add_angle_constraint(*cluster, i1, i2, i3, i4, compliance, minAngle, maxAngle);
                }
                else if (ac.contains("restAngle"))
                {
                    const float restAngle = glm::radians(ac.value("restAngle", 0.0f));

                    xpbd::cluster::add_angle_constraint(*cluster, i1, i2, i3, i4, compliance, restAngle, restAngle);
                }
                else

                    xpbd::cluster::add_angle_constraint(*cluster, i1, i2, i3, i4, compliance);
            }
        }

        // Load Volume Constraints.
        if (j.contains("VolumeConstraints") &&
            j["VolumeConstraints"].is_array())
        {
            for (const auto &vc : j["VolumeConstraints"])
            {
                std::vector<size_t> indices;

                for (const auto &idx : vc["indices"])
                    indices.push_back(get_index(*cluster, idx, base_index));

                const float compliance = vc.value("compliance", 0.0f);
                
                if (vc.contains("restPressure"))
                {
                    const float restPressure = vc.value("restPressure", 1.0f);

                    xpbd::cluster::add_volume_constraint(*cluster, indices, compliance, restPressure);
                }
                else
                    xpbd::cluster::add_volume_constraint(*cluster, indices, compliance);
            }
        }

        // Load Polygon Colliders.
        if (j.contains("PolygonColliders") &&
            j["PolygonColliders"].is_array())
        {
            for (const auto &cc : j["PolygonColliders"])
            {
                std::vector<size_t> indices;

                for (const auto &idx : cc["indices"])
                    indices.push_back(get_index(*cluster, idx, base_index));

                const float staticFriction = cc.value("staticFriction", 0.0f);
                const float kineticFriction = cc.value("kineticFriction", 0.0f);
                const float compliance = cc.value("compliance", 0.0f);
                const int8_t layer = cc.value("layer", 1);
                const int8_t mask = cc.value("mask", 1);

                xpbd::cluster::add_polygon_collider(
                    *cluster,
                    indices,
                    staticFriction,
                    kineticFriction,
                    compliance,
                    layer,
                    mask);
            }
        }

        // Load Points Colliders.
        if (j.contains("PointsColliders") &&
            j["PointsColliders"].is_array())
        {
            for (const auto &cc : j["PointsColliders"])
            {
                std::vector<size_t> indices;

                for (const auto &idx : cc["indices"])
                    indices.push_back(get_index(*cluster, idx, base_index));

                const float staticFriction = cc.value("staticFriction", 0.0f);
                const float kineticFriction = cc.value("kineticFriction", 0.0f);
                const float compliance = cc.value("compliance", 0.0f);
                const int8_t layer = cc.value("layer", 1);
                const int8_t mask = cc.value("mask", 1);

                xpbd::cluster::add_points_collider(
                    *cluster,
                    indices,
                    staticFriction,
                    kineticFriction,
                    compliance,
                    layer,
                    mask);
            }
        }

        // Load Angular Motors.
        if (j.contains("AngularMotors") &&
            j["AngularMotors"].is_array())
        {
            for (const auto &m : j["AngularMotors"])
            {
                const size_t i1 = get_index(*cluster, m["i1"], base_index);
                const size_t i2 = get_index(*cluster, m["i2"], base_index);
                const size_t i3 = get_index(*cluster, m["i3"], base_index);
                const size_t i4 = get_index(*cluster, m["i4"], base_index);

                const float targetRPM = m.value("targetRPM", 0.0f);
                const float maxTorque = m.value("maxTorque", 100.0f);

                xpbd::cluster::add_angular_motor(
                    *cluster,
                    i1, i2, i3, i4,
                    targetRPM,
                    maxTorque);
            }
        }

        // Load Angular Servos.
        if (j.contains("AngularServos") &&
            j["AngularServos"].is_array())
        {
            for (const auto &s : j["AngularServos"])
            {
                const size_t i1 = get_index(*cluster, s["i1"], base_index);
                const size_t i2 = get_index(*cluster, s["i2"], base_index);
                const size_t i3 = get_index(*cluster, s["i3"], base_index);
                const size_t i4 = get_index(*cluster, s["i4"], base_index);

                const float targetAngle = glm::radians(s.value("targetAngle", 0.0f));
                const float kp = s.value("kp", 0.0f);
                const float kd = s.value("kd", 0.0f);
                const float maxTorque = s.value("maxTorque", 0.0f);

                xpbd::cluster::add_angular_servo(
                    *cluster,
                    i1, i2, i3, i4,
                    targetAngle,
                    kp,
                    kd,
                    maxTorque);
            }
        }

        return cluster;
    }
    catch (const std::exception &e)
    {
        std::println(
            "Failed to load body '{}': {}",
            filename,
            e.what());

        delete cluster;
        return nullptr;
    }
}

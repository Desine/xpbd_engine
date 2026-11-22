#include "json_body_loader.hpp"
#include "xpbd.hpp"

#include "nlohmann/json.hpp"
#include <fstream>
#include <print>


std::string json_body_loader::folder = "";
const std::string json_body_loader::extension = ".body";

inline glm::vec2 get_vec2(const nlohmann::json &p, const char *key, glm::vec2 def = {0, 0})
{
    auto j = p.value(key, nlohmann::json::array({def.x, def.y}));
    return {j[0].get<float>(), j[1].get<float>()};
}
std::optional<size_t> get_index_by_id(xpbd::Cluster &cluster, std::string id)
{
    auto result = std::find_if(cluster.particleIds.begin(), cluster.particleIds.end(), [&id](const auto &p)
                               { return p.first == id; });
    if (result != cluster.particleIds.end())
        return result->second;
    return std::nullopt;
}
size_t get_index(xpbd::Cluster &cluster, const nlohmann::json &element, size_t base_index)
{
    size_t index = 0;
    if (element.is_string())
    {
        if (auto idFound = get_index_by_id(cluster, element.get<std::string>()))
            index = idFound.value();
    }
    else
        index = element.get<size_t>() + base_index;
    return index;
}

bool json_body_loader::load(xpbd::Cluster &cluster, const std::string &filename, glm::vec2 offset = {0.0f, 0.0f})
{
    std::ifstream file(folder + filename + extension);
    if (!file.is_open())
    {
        std::println("Failed to open file: {}|{}|{}", folder, filename, extension);
        return false;
    }
    else
        std::println("opened file {}", folder + filename + extension);
    nlohmann::json j;
    file >> j;

    // Load Required bodies
    if (j.contains("Require") && j["Require"].is_array())
    {
        for (const auto &r : j["Require"])
        {
            glm::vec2 r_offset = get_vec2(r, "offset") + offset;

            if (r.contains("name"))
                if (load(cluster, r["name"].get<std::string>(), r_offset) == false)
                    return false;
        }
    }

    size_t base_index = cluster.particles.size(); // offset for constraint indices

    // Load Particles
    if (j.contains("Particles") && j["Particles"].is_array())
    {
        glm::vec2 j_offset = get_vec2(j, "offset");
        glm::vec2 j_velocity = get_vec2(j, "velocity");

        for (const auto &&[index, p] : std::views::enumerate(j["Particles"]))
        {
            glm::vec2 pos = get_vec2(p, "pos") + offset + j_offset;
            glm::vec2 vel = get_vec2(p, "vel") + j_velocity;
            float mass = p.value("mass", 0.0f);
            xpbd::add_particle(cluster, pos, mass, vel);

            if (p.contains("id"))
                cluster.particleIds.push_back(std::pair(p["id"].get<std::string>(), index + base_index));
        }
    }

    // Load Distance Constraints
    if (j.contains("DistanceConstraints") && j["DistanceConstraints"].is_array())
    {
        for (const auto &dc : j["DistanceConstraints"])
        {
            size_t i1 = get_index(cluster, dc["i1"], base_index);
            size_t i2 = get_index(cluster, dc["i2"], base_index);
            float compliance = dc.value("compliance", 0.0f);
            if (dc.contains("minDist") && dc.contains("maxDist"))
            {
                float minDist = dc.value("minDist", 0.0f);
                float maxDist = dc.value("maxDist", 0.0f);
                xpbd::add_distance_constraint(cluster, i1, i2, compliance, minDist, maxDist);
            }
            else if (dc.contains("restDist"))
            {
                float restDist = dc.value("restDist", 0.0f);
                xpbd::add_distance_constraint(cluster, i1, i2, compliance, restDist, restDist);
            }
            else
                xpbd::add_distance_constraint(cluster, i1, i2, compliance);
        }
    }

    // Load Angle Constraints
    if (j.contains("AngleConstraints") && j["AngleConstraints"].is_array())
    {
        for (const auto &ac : j["AngleConstraints"])
        {
            size_t i1 = get_index(cluster, ac["i1"], base_index);
            size_t i2 = get_index(cluster, ac["i2"], base_index);
            size_t i3 = get_index(cluster, ac["i3"], base_index);
            size_t i4 = get_index(cluster, ac["i4"], base_index);
            float compliance = ac.value("compliance", 0.0f);
            if (ac.contains("minAngle") && ac.contains("maxAngle"))
            {
                float minAngle = glm::radians(ac.value("minAngle", 0.0f));
                float maxAngle = glm::radians(ac.value("maxAngle", 0.0f));
                xpbd::add_angle_constraint(cluster, i1, i2, i3, i4, compliance, minAngle, maxAngle);
            }
            else if (ac.contains("restAngle"))
            {
                float restAngle = glm::radians(ac.value("restAngle", 0.0f));
                xpbd::add_angle_constraint(cluster, i1, i2, i3, i4, compliance, restAngle, restAngle);
            }
            else
                xpbd::add_angle_constraint(cluster, i1, i2, i3, i4, compliance);
        }
    }

    // Load Volume Constraints
    if (j.contains("VolumeConstraints") && j["VolumeConstraints"].is_array())
    {
        for (const auto &vc : j["VolumeConstraints"])
        {
            std::vector<size_t> indices;
            for (const auto &idx : vc["indices"])
                indices.push_back(get_index(cluster, idx, base_index));

            float compliance = vc.value("compliance", 0.0f);
            if (vc.contains("restPressure"))
            {
                float restVolume = compute_polygon_area(cluster.particles, indices) * vc.value("restPressure", 1.0f);
                xpbd::add_volume_constraint(cluster, indices, compliance, restVolume);
            }
            else
                xpbd::add_volume_constraint(cluster, indices, compliance);
        }
    }

    // Load Polygon Colliders
    if (j.contains("PolygonColliders") && j["PolygonColliders"].is_array())
    {
        for (const auto &cc : j["PolygonColliders"])
        {
            std::vector<size_t> indices;
            for (const auto &idx : cc["indices"])
                indices.push_back(get_index(cluster, idx, base_index));

            float staticFriction = cc.value("staticFriction", 0.0f);
            float kineticFriction = cc.value("kineticFriction", 0.0f);
            float compliance = cc.value("compliance", 0.0f);
            int8_t layer = cc.value("layer", 1);
            int8_t mask = cc.value("mask", 1);

            xpbd::add_polygon_collider(cluster, indices, staticFriction, kineticFriction, compliance, layer, mask);
        }
    }

    // Load Point Colliders
    if (j.contains("PointsColliders") && j["PointsColliders"].is_array())
    {
        for (const auto &cc : j["PointsColliders"])
        {
            std::vector<size_t> indices;
            for (const auto &idx : cc["indices"])
                indices.push_back(get_index(cluster, idx, base_index));

            float staticFriction = cc.value("staticFriction", 0.0f);
            float kineticFriction = cc.value("kineticFriction", 0.0f);
            float compliance = cc.value("compliance", 0.0f);
            int8_t layer = cc.value("layer", 1);
            int8_t mask = cc.value("mask", 1);

            xpbd::add_points_collider(cluster, indices, staticFriction, kineticFriction, compliance, layer, mask);
        }
    }

    // Load Angular Motors
    if (j.contains("AngularMotors") && j["AngularMotors"].is_array())
    {
        for (const auto &m : j["AngularMotors"])
        {
            size_t i1 = get_index(cluster, m["i1"], base_index);
            size_t i2 = get_index(cluster, m["i2"], base_index);
            size_t i3 = get_index(cluster, m["i3"], base_index);
            size_t i4 = get_index(cluster, m["i4"], base_index);

            float targetRPM = m.value("targetRPM", 0.0f);
            float maxTorque = m.value("maxTorque", 100.0f);

            xpbd::add_angular_motor(
                cluster,
                i1, i2, i3, i4,
                targetRPM,
                maxTorque);
        }
    }

    // Load Angular Servos
    if (j.contains("AngularServos") && j["AngularServos"].is_array())
    {
        for (const auto &s : j["AngularServos"])
        {
            size_t i1 = get_index(cluster, s["i1"], base_index);
            size_t i2 = get_index(cluster, s["i2"], base_index);
            size_t i3 = get_index(cluster, s["i3"], base_index);
            size_t i4 = get_index(cluster, s["i4"], base_index);

            float targetAngle = glm::radians(s.value("targetAngle", 0.0f));
            float kp = s.value("kp", 0.0f);
            float kd = s.value("kd", 0.0f);
            float maxTorque = s.value("maxTorque", 0.0f);

            xpbd::add_angular_servo(
                cluster,
                i1, i2, i3, i4,
                targetAngle,
                kp,
                kd,
                maxTorque);
        }
    }

    return true;
}

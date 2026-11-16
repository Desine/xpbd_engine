#include "json_body_loader.hpp"
#include "nlohmann/json.hpp"
#include <fstream>

#include "xpbd.hpp"

std::string json_body_loader::folder = "";
const std::string json_body_loader::extension = ".body";

std::optional<size_t> get_index_by_id(xpbd::Cluster &cluster, std::string id)
{
    auto result = std::find_if(cluster.particleIds.begin(), cluster.particleIds.end(), [&id](const auto &p)
                               { return p.first == id; });
    if (result != cluster.particleIds.end())
        return result->second;
    return NULL;
}

size_t get_index(xpbd::Cluster &cluster, const nlohmann::json_abi_v3_11_3::json &element, size_t base_index)
{
    if (element.is_string())
    {
        if (auto idFound = get_index_by_id(cluster, element.get<std::string>()))
            return idFound.value();
        }
    else
        return element.get<size_t>() + base_index;
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
            glm::vec2 r_offset = offset;
            if (r.contains("offset"))
                r_offset += glm::vec2(r["offset"][0].get<float>(), r["offset"][1].get<float>());

            if (r.contains("name"))
                if (load(cluster, r["name"].get<std::string>(), r_offset) == false)
                    return false;
        }
    }

    size_t base_index = cluster.particles.pos.size(); // offset for constraint indices

    // Load Particles
    if (j.contains("Particles") && j["Particles"].is_array())
    {
        glm::vec2 j_offset = {0, 0};
        if (j.contains("offset") && j["offset"].is_array() && j["offset"].size() == 2)
            j_offset = {j["offset"][0].get<float>(), j["offset"][1].get<float>()};

        for (const auto &&[index, p] : std::views::enumerate(j["Particles"]))
        {
            if (!p.contains("pos") || !p.contains("mass"))
                continue;

            glm::vec2 pos(p["pos"][0].get<float>(), p["pos"][1].get<float>());
            pos += j_offset + offset;

            glm::vec2 vel = {0, 0};
            if (p.contains("vel"))
                vel = {p["vel"][0].get<float>(), p["vel"][1].get<float>()};

            float mass = p["mass"].get<float>();
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
            float compliance = dc["compliance"].get<float>();
            if (dc.contains("restDistance"))
                xpbd::add_distance_constraint(cluster, i1, i2, compliance, dc["restDistance"].get<float>());
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
            float compliance = ac["compliance"].get<float>();
            if (ac.contains("restAngle"))
                xpbd::add_angle_constraint(cluster, i1, i2, i3, i4, compliance, glm::radians(ac["restAngle"].get<float>()));
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

            float compliance = vc["compliance"].get<float>();
            if (vc.contains("restPressure"))
                xpbd::add_volume_constraint(cluster, indices, compliance, vc["restPressure"].get<float>());
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

            float staticFriction = cc["staticFriction"].get<float>();
            float kineticFriction = cc["kineticFriction"].get<float>();
            float compliance = cc["compliance"].get<float>();
            int8_t layer = cc.contains("layer") ? cc["layer"].get<int8_t>() : 1;
            int8_t mask = cc.contains("mask") ? cc["mask"].get<int8_t>() : 1;

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

            float staticFriction = cc["staticFriction"].get<float>();
            float kineticFriction = cc["kineticFriction"].get<float>();
            float compliance = cc["compliance"].get<float>();
            int8_t layer = cc.contains("layer") ? cc["layer"].get<int8_t>() : 1;
            int8_t mask = cc.contains("mask") ? cc["mask"].get<int8_t>() : 1;

            xpbd::add_points_collider(cluster, indices, staticFriction, kineticFriction, compliance, layer, mask);
        }
    }

    return true;
}

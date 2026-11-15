#include "json_body_loader.hpp"
#include "nlohmann/json.hpp"
#include <fstream>

#include "xpbd.hpp"

std::string json_body_loader::folder = "";
const std::string json_body_loader::extension = ".body";

void json_body_loader::load(
    xpbd::Cluster &cluster,
    const std::string &filename,
    glm::vec2 offset = {0.0f, 0.0f})
{
    std::ifstream file(folder + filename + extension);
    if (!file.is_open())
        throw std::runtime_error("Failed to open file: |" + folder + "|" + filename + "|" + extension);
    else
        std::println("opened file {}", folder + filename + extension);

    nlohmann::json j;
    file >> j;

    glm::vec2 read_offset = {0, 0};
    if (j.contains("position") && j["position"].is_array() && j["position"].size() == 2)
    {
        read_offset = {j["position"][0].get<float>(), j["position"][1].get<float>()};
    }

    size_t base_index = cluster.particles.pos.size(); // offset for constraint indices

    // Load Particles
    if (j.contains("Particles") && j["Particles"].is_array())
    {
        for (const auto &p : j["Particles"])
        {
            if (!p.contains("pos") || !p.contains("mass"))
                continue;

            glm::vec2 pos(p["pos"][0].get<float>(), p["pos"][1].get<float>());
            pos += read_offset + offset;

            glm::vec2 vel = {0, 0};
            if (p.contains("vel"))
                vel = {p["vel"][0].get<float>(), p["vel"][1].get<float>()};

            float mass = p["mass"].get<float>();
            xpbd::add_particle(cluster, pos, mass, vel);
        }
    }

    // Load Distance Constraints
    if (j.contains("DistanceConstraints") && j["DistanceConstraints"].is_array())
    {
        for (const auto &dc : j["DistanceConstraints"])
        {
            size_t i1 = dc["i1"].get<size_t>() + base_index;
            size_t i2 = dc["i2"].get<size_t>() + base_index;
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
            size_t i1 = ac["i1"].get<size_t>() + base_index;
            size_t i2 = ac["i2"].get<size_t>() + base_index;
            size_t i3 = ac["i3"].get<size_t>() + base_index;
            size_t i4 = ac["i4"].get<size_t>() + base_index;
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
                indices.push_back(idx.get<size_t>() + base_index);

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
                indices.push_back(idx.get<size_t>() + base_index);

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
                indices.push_back(idx.get<size_t>() + base_index);

            float staticFriction = cc["staticFriction"].get<float>();
            float kineticFriction = cc["kineticFriction"].get<float>();
            float compliance = cc["compliance"].get<float>();
            int8_t layer = cc.contains("layer") ? cc["layer"].get<int8_t>() : 1;
            int8_t mask = cc.contains("mask") ? cc["mask"].get<int8_t>() : 1;

            xpbd::add_points_collider(cluster, indices, staticFriction, kineticFriction, compliance, layer, mask);
        }
    }

    // Load Required bodies
    if (j.contains("Require") && j["Require"].is_array())
    {
        for (const auto &r : j["Require"])
        {
            glm::vec2 r_offset = offset;
            if (r.contains("offset"))
                r_offset += glm::vec2(r["offset"][0].get<float>(), r["offset"][1].get<float>());

            if (r.contains("name"))
                load(cluster, r["name"].get<std::string>(), r_offset);
        }
    }
}

#include "json_body_loader.hpp"
#include <fstream>
#include <nlohmann/json.hpp>

#include "defines.hpp"
#include "xpbd.hpp"

using json = nlohmann::json;

namespace xpbd
{
    class World;
}

namespace json_body_loader
{
    void load(
        xpbd::World &world,
        const std::string &filename,
        glm::vec2 offset = {0.0f, 0.0f})
    {
        std::ifstream file(folder + filename + extension);
        if (!file.is_open())
            throw std::runtime_error("Failed to open file: " + folder + filename + extension);

        json j;
        file >> j;

        glm::vec2 position_offset = {0, 0};
        if (j.contains("position") && j["position"].is_array() && j["position"].size() == 2)
        {
            position_offset = {j["position"][0].get<float>(), j["position"][1].get<float>()};
        }

        size_t base_index = world.particles.pos.size(); // offset for constraint indices

        // Load Particles
        if (j.contains("Particles") && j["Particles"].is_array())
        {
            for (const auto &p : j["Particles"])
            {
                if (!p.contains("pos") || !p.contains("mass"))
                    continue;

                glm::vec2 pos(p["pos"][0].get<float>(), p["pos"][1].get<float>());
                pos += position_offset + offset;

                glm::vec2 vel = {0, 0};
                if (p.contains("vel"))
                    vel = {p["vel"][0].get<float>(), p["vel"][1].get<float>()};

                float mass = p["mass"].get<float>();
                world.add_particle(pos, mass, vel);
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
                float restDist = dc.contains("restDistance") ? dc["restDistance"].get<float>() : FLOAT_DEFAULT;
                world.add_distance_constraint(i1, i2, compliance, restDist);
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
                float restPressure = vc.contains("restPressure") ? vc["restPressure"].get<float>() : FLOAT_DEFAULT;
                world.add_volume_constraint(indices, compliance, restPressure);
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

                world.add_polygon_collider(indices, staticFriction, kineticFriction, compliance);
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

                world.add_points_collider(indices, staticFriction, kineticFriction, compliance);
            }
        }
    }
}

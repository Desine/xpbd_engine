#include <fstream>
#include <string.h>
#include "glm/glm.hpp"
#include <nlohmann/json.hpp>

#include "xpbd.hpp"
using json = nlohmann::json;

namespace json_body_loader
{
    const std::string extension = ".body";

    void load(
        const std::string &filename,
        xpbd::Particles &particles,
        xpbd::DistanceConstraints &distanceConstraints,
        xpbd::VolumeConstraints &volumeConstraints,
        xpbd::ContourColliders &contourColliders)
    {
        std::ifstream file(filename + extension);
        if (!file.is_open())
            throw std::runtime_error("Failed to open file: " + filename + extension);
    
        json j;
        file >> j;
    
        glm::vec2 position_offset = {0, 0};
        if (j.contains("position") && j["position"].is_array() && j["position"].size() == 2)
        {
            position_offset = {j["position"][0].get<float>(), j["position"][1].get<float>()};
        }
    
        size_t base_index = particles.pos.size(); // offset for constraint indices
    
        // Load Particles
        if (j.contains("Particles") && j["Particles"].is_array())
        {
            for (const auto &p : j["Particles"])
            {
                if (!p.contains("pos") || !p.contains("mass"))
                    continue;
    
                glm::vec2 pos(p["pos"][0].get<float>(), p["pos"][1].get<float>());
                pos += position_offset;
    
                glm::vec2 vel = {0, 0};
                if (p.contains("vel"))
                    vel = {p["vel"][0].get<float>(), p["vel"][1].get<float>()};
    
                float mass = p["mass"].get<float>();
                xpbd::add_particle(particles, pos, mass, vel);
            }
        }
    
        // Load Distance Constraints
        if (j.contains("DistanceConstraints") && j["DistanceConstraints"].is_array())
        {
            for (const auto &dc : j["DistanceConstraints"])
            {
                size_t i1 = dc["i1"].get<size_t>() + base_index;
                size_t i2 = dc["i2"].get<size_t>() + base_index;
                float restDist = dc["restDistance"].get<float>();
                float compliance = dc["compliance"].get<float>();
                xpbd::add_distance_constraint(distanceConstraints, i1, i2, compliance, restDist);
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
    
                float restPressure = vc["restPressure"].get<float>();
                float compliance = vc["compliance"].get<float>();
                xpbd::add_volume_constraint(particles, volumeConstraints, indices, compliance, restPressure);
            }
        }
    
        // Load Contour Colliders
        if (j.contains("ContourColliders") && j["ContourColliders"].is_array())
        {
            for (const auto &cc : j["ContourColliders"])
            {
                std::vector<size_t> indices;
                for (const auto &idx : cc["indices"])
                    indices.push_back(idx.get<size_t>() + base_index);
    
                float staticFriction = cc["staticFriction"].get<float>();
                float kineticFriction = cc["kineticFriction"].get<float>();
                float compliance = cc["compliance"].get<float>();
    
                xpbd::add_contour_collider(contourColliders, indices, staticFriction, kineticFriction, compliance);
            }
        }
    }

}
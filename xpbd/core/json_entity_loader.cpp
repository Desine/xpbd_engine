#include "scene.hpp"

#include "nlohmann/json.hpp"
#include <filesystem>
#include <fstream>
#include <print>

namespace xpbd
{

using json = nlohmann::json;
using JsonVariables = std::map<std::string, std::string>;

void read_variables(
    const std::filesystem::path &folder,
    const std::string &file_name,
    const std::string &extension,
    JsonVariables &variables,
    const std::string &prefix = "",
    bool allow_includes = true)
{
    auto path = folder / (file_name + extension);
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + path.string());

    std::string line;

    while (std::getline(file, line))
    {
        if (line.find('{') != std::string::npos)
            break;

        // Удаляем пробелы по краям.
        auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            continue;

        auto last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);

        // Инклюд: $car_body
        if (line[0] == '$')
        {
            if (!allow_includes)
                continue;

            std::string include_name = line.substr(1);

            if (include_name.empty())
                throw std::runtime_error("Empty include name");

            read_variables(
                folder,
                include_name,
                extension,
                variables,
                include_name + ".",
                false);

            continue;
        }

        auto separator = line.find(':');
        if (separator == std::string::npos)
            continue;

        std::string key = line.substr(0, separator);
        std::string value = line.substr(separator + 1);

        if (key.empty())
            throw std::runtime_error("Empty variable name: " + path.string());

        for (char c : key)
        {
            if (std::isspace(static_cast<unsigned char>(c)))
                throw std::runtime_error("Whitespace in variable name: " + key);
        }

        std::string full_key = prefix + key;

        variables[full_key] = value;
    }
}

std::string replace_variables(
    const std::string &source,
    const JsonVariables &variables)
{
    std::string result;

    for (size_t i = 0; i < source.size();)
    {
        if (source[i] != '$')
        {
            result += source[i++];
            continue;
        }

        size_t start = i++;

        while (i < source.size())
        {
            unsigned char c = static_cast<unsigned char>(source[i]);

            if (std::isalnum(c) || c == '_' || c == '.')
                ++i;
            else
                break;
        }

        std::string key = source.substr(start + 1, i - start - 1);

        if (key.empty())
            throw std::runtime_error("Empty variable reference");

        auto it = variables.find(key);
        if (it == variables.end())
            throw std::runtime_error("no such key: $" + key);

        result += it->second;
    }

    return result;
}

json parse_entity_json(
    const std::string &folder,
    const std::string &file_name,
    const std::string &extension)
{
    JsonVariables variables;

    read_variables(
        folder,
        file_name,
        extension,
        variables);

    auto path = std::filesystem::path(folder) /
                (file_name + extension);

    std::ifstream file(path);

    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + path.string());

    std::stringstream buffer;
    std::string line;

    bool json_started = false;

    while (std::getline(file, line))
    {
        if (!json_started)
        {
            if (line.find('{') == std::string::npos)
                continue;

            json_started = true;

            auto pos = line.find('{');
            line = line.substr(pos);
        }

        buffer << line << '\n';
    }

    if (!json_started)
        throw std::runtime_error("JSON object not found: " + path.string());

    std::string json_text = buffer.str();
    json_text = replace_variables(json_text, variables);

    try
    {
        return json::parse(json_text);
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error(
            std::string("JSON parse error: ") + e.what());
    }
}

/*





*/

inline glm::vec2 get_vec2(
    const nlohmann::json &json,
    const std::string &key,
    glm::vec2 def = {0, 0})
{
    if (!json.contains(key))
        return def;

    const auto &j = json[key];

    if (!j.is_array() || j.size() != 2)
        throw std::runtime_error(
            "invalid vec2: " + key);

    return {
        j[0].get<float>(),
        j[1].get<float>()};
}

ParticleId get_particle_id_by_name(const ParticleStorage &storage, const std::vector<ParticleId> &ids, const std::string &name)
{
    for (auto &id : ids)
    {
        const std::string *storage_name = storage.get_name(id);
        if (!storage_name)
        {
            std::println("id: {}, has no storage_name", id.index);
            continue;
        }
        if (*storage_name == name)
            return id;
    }
    return {};
}
ParticleId get_particle_id_by_path(
    const EntityStorage &entity_storage,
    const ParticleStorage &particle_storage,
    const EntityId &root_id,
    const std::string &path)
{
    const Entity *root_entity = entity_storage.get(root_id);
    if (!root_entity) return {};

    const size_t separator = path.find('.');

    if (separator == std::string::npos)
    {
        return get_particle_id_by_name(
            particle_storage,
            root_entity->particles,
            path);
    }

    const std::string entity_name = path.substr(0, separator);
    const std::string particle_name = path.substr(separator + 1);

    for (EntityId child_id : root_entity->children)
    {
        const Entity *child = entity_storage.get(child_id);
        if (!child)
            continue;

        const std::string *child_name = entity_storage.get_name(child_id);

        if (child_name && *child_name == entity_name)
        {
            return get_particle_id_by_name(
                particle_storage,
                child->particles,
                particle_name);
        }
    }

    return {};
}

ParticleId get_particle_id(
    const nlohmann::json &value,
    const EntityId &root_id,
    const EntityStorage &entity_storage,
    const ParticleStorage &particle_storage)
{
    const Entity *root_entity = entity_storage.get(root_id);
    if (!root_entity) return {};

    if (value.is_number_integer())
    {
        const int index = value.get<int>();

        if (index < 0 ||
            static_cast<size_t>(index) >= root_entity->particles.size())
            return {};

        return root_entity->particles[index];
    }

    if (value.is_string())
    {
        return get_particle_id_by_path(
            entity_storage,
            particle_storage,
            root_id,
            value.get<std::string>());
    }

    return {};
}
ParticleId get_particle_id(
    const nlohmann::json &j,
    const char *key,
    const EntityId &root_id,
    const EntityStorage &entity_storage,
    const ParticleStorage &particle_storage)
{
    if (!j.contains(key))
        return {};

    const auto &value = j[key];

    return get_particle_id(value, root_id, entity_storage, particle_storage);
}

EntityId Scene::spawn_entity(
    const std::string &file_name,
    const glm::vec2 &position,
    std::string name,
    EntityId root_id)
{
    const std::string folder = std::filesystem::current_path().generic_string() + "/body/";
    const std::string extension = ".entity";
    if (file_name.empty()) return {};

    EntityId entity_id = add_entity();
    if (!entity_id.valid()) return {};
    if (!root_id.valid()) root_id = entity_id;

    try
    {
        nlohmann::json j = parse_entity_json(folder, file_name, extension);

        if (name.empty())
            name = j.value("name", file_name);
        entity_storage.set_name(entity_id, name);

        if (j.contains("entities") && j["entities"].is_array())
        {
            for (const auto &entity : j["entities"])
            {
                std::string entity_name = entity.value("name", "");
                std::string entity_file_name = entity.value("file_name", "");
                glm::vec2 offset = get_vec2(entity, "offset");

                std::println(
                    "entity_file_name: {}, offset: {}, entity_name: {}",
                    entity_file_name, offset[0], entity_name);

                EntityId childId = spawn_entity(
                    entity_file_name,
                    offset + position,
                    entity_name);

                if (!childId.valid())
                    throw std::runtime_error(
                        "Failed to spawn child entity: " + entity_file_name);

                Entity *parent = entity_storage.get(root_id);
                if (!parent) throw std::runtime_error("!parent");
                parent->children.push_back(childId);
            }
        }

        // Load Particles.
        if (j.contains("Particles") && j["Particles"].is_array())
        {
            const glm::vec2 j_offset = get_vec2(j, "offset");
            const glm::vec2 j_velocity = get_vec2(j, "velocity");

            size_t index = 0;

            for (const auto &p : j["Particles"])
            {
                const glm::vec2 pos = get_vec2(p, "pos") + j_offset + position;
                const glm::vec2 vel = get_vec2(p, "vel") + j_velocity;
                const float mass = p.value("mass", 0.0f);

                ParticleId particle_id = add_particle(pos, mass, vel);
                std::string particle_name = p.value("name", std::to_string(index));
                if (!particle_storage.set_name(particle_id, particle_name))
                    throw std::runtime_error("can't set particle name: " + particle_name + ", index: " + std::to_string(index));

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->particles.push_back(particle_id);

                index++;
            }
        }

        // Load Distance Constraints.
        if (j.contains("DistanceConstraints") && j["DistanceConstraints"].is_array())
        {
            for (const auto &dc : j["DistanceConstraints"])
            {
                const ParticleId p1 = get_particle_id(dc, "p1", root_id, entity_storage, particle_storage);
                const ParticleId p2 = get_particle_id(dc, "p2", root_id, entity_storage, particle_storage);

                const float compliance = dc.value("compliance", 0.0f);

                DistanceConstraintId distance_constraint_id;

                if (dc.contains("min_dist") && dc.contains("max_dist"))
                {
                    const float min_dist = dc.value("min_dist", 0.0f);
                    const float max_dist = dc.value("max_dist", 0.0f);

                    distance_constraint_id = add_distance_constraint(p1, p2, compliance, min_dist, max_dist);
                }
                else if (dc.contains("restDist"))
                {
                    const float restDist = dc.value("restDist", 0.0f);

                    distance_constraint_id = add_distance_constraint(p1, p2, compliance, restDist, restDist);
                }
                else
                {
                    distance_constraint_id = add_distance_constraint(p1, p2, compliance);
                }

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->distance_constraints.push_back(distance_constraint_id);
            }
        }

        // Load Angle Constraints.
        if (j.contains("AngleConstraints") && j["AngleConstraints"].is_array())
        {
            for (const auto &ac : j["AngleConstraints"])
            {
                const ParticleId p1 = get_particle_id(ac, "p1", root_id, entity_storage, particle_storage);
                const ParticleId p2 = get_particle_id(ac, "p2", root_id, entity_storage, particle_storage);
                const ParticleId p3 = get_particle_id(ac, "p3", root_id, entity_storage, particle_storage);
                const ParticleId p4 = get_particle_id(ac, "p4", root_id, entity_storage, particle_storage);

                const float compliance = ac.value("compliance", 0.0f);

                AngleConstraintId angle_constraint_id;

                if (ac.contains("min_angle") && ac.contains("max_angle"))
                {
                    const float min_angle = glm::radians(ac.value("min_angle", 0.0f));
                    const float max_angle = glm::radians(ac.value("max_angle", 0.0f));

                    angle_constraint_id = add_angle_constraint(p1, p2, p3, p4, compliance, min_angle, max_angle);
                }
                else if (ac.contains("restAngle"))
                {
                    const float restAngle = glm::radians(ac.value("restAngle", 0.0f));

                    angle_constraint_id = add_angle_constraint(p1, p2, p3, p4, compliance, restAngle, restAngle);
                }
                else
                {
                    angle_constraint_id = add_angle_constraint(p1, p2, p3, p4, compliance);
                }

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->angle_constraints.push_back(angle_constraint_id);
            }
        }

        // Load Volume Constraints.
        if (j.contains("VolumeConstraints") && j["VolumeConstraints"].is_array())
        {
            for (const auto &vc : j["VolumeConstraints"])
            {
                std::vector<ParticleId> indices;
                for (const auto &idx : vc["indices"])
                    indices.push_back(get_particle_id(idx, root_id, entity_storage, particle_storage));

                const float compliance = vc.value("compliance", 0.0f);

                VolumeConstraintId volume_constraint_id;

                if (vc.contains("rest_pressure"))
                {
                    const float rest_pressure = vc.value("rest_pressure", 1.0f);

                    volume_constraint_id = add_volume_constraint(indices, compliance, rest_pressure);
                }
                else
                {
                    volume_constraint_id = add_volume_constraint(indices, compliance);
                }

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->volume_constraints.push_back(volume_constraint_id);
            }
        }

        // Load Polygon Colliders.
        if (j.contains("PolygonColliders") && j["PolygonColliders"].is_array())
        {
            for (const auto &cc : j["PolygonColliders"])
            {
                std::vector<ParticleId> indices;
                for (const auto &idx : cc["indices"])
                    indices.push_back(get_particle_id(idx, root_id, entity_storage, particle_storage));

                const float static_friction = cc.value("static_friction", 0.0f);
                const float kinetic_friction = cc.value("kinetic_friction", 0.0f);
                const float compliance = cc.value("compliance", 0.0f);
                const uint32_t layer = cc.value("layer", 1);
                const uint32_t mask = cc.value("mask", 1);

                PolygonColliderId polygon_collider_id;

                polygon_collider_id = add_polygon_collider(
                    indices,
                    static_friction,
                    kinetic_friction,
                    compliance,
                    layer,
                    mask);

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->polygon_colliders.push_back(polygon_collider_id);
            }
        }

        // Load Points Colliders.
        if (j.contains("PointsColliders") && j["PointsColliders"].is_array())
        {
            for (const auto &cc : j["PointsColliders"])
            {
                std::vector<ParticleId> indices;
                for (const auto &idx : cc["indices"])
                    indices.push_back(get_particle_id(idx, root_id, entity_storage, particle_storage));

                const float static_friction = cc.value("static_friction", 0.0f);
                const float kinetic_friction = cc.value("kinetic_friction", 0.0f);
                const float compliance = cc.value("compliance", 0.0f);
                const uint32_t layer = cc.value("layer", 1);
                const uint32_t mask = cc.value("mask", 1);

                PointsColliderId points_collider_id;

                points_collider_id = add_points_collider(
                    indices,
                    static_friction,
                    kinetic_friction,
                    compliance,
                    layer,
                    mask);

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->points_colliders.push_back(points_collider_id);
            }
        }

        // Load Angular Motors.
        if (j.contains("AngularMotors") && j["AngularMotors"].is_array())
        {
            for (const auto &m : j["AngularMotors"])
            {
                const ParticleId p1 = get_particle_id(m, "p1", root_id, entity_storage, particle_storage);
                const ParticleId p2 = get_particle_id(m, "p2", root_id, entity_storage, particle_storage);
                const ParticleId p3 = get_particle_id(m, "p3", root_id, entity_storage, particle_storage);
                const ParticleId p4 = get_particle_id(m, "p4", root_id, entity_storage, particle_storage);

                const float target_rpm = m.value("target_rpm", 0.0f);
                const float max_torque = m.value("max_torque", 100.0f);

                AngularMotorId angular_motor_id;

                angular_motor_id = add_angular_motor(
                    p1, p2, p3, p4,
                    target_rpm,
                    max_torque);

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->angular_motors.push_back(angular_motor_id);
            }
        }

        // Load Angular Servos.
        if (j.contains("AngularServos") && j["AngularServos"].is_array())
        {
            for (const auto &s : j["AngularServos"])
            {
                const ParticleId p1 = get_particle_id(s, "p1", root_id, entity_storage, particle_storage);
                const ParticleId p2 = get_particle_id(s, "p2", root_id, entity_storage, particle_storage);
                const ParticleId p3 = get_particle_id(s, "p3", root_id, entity_storage, particle_storage);
                const ParticleId p4 = get_particle_id(s, "p4", root_id, entity_storage, particle_storage);

                const float target_angle = glm::radians(s.value("target_angle", 0.0f));
                const float kp = s.value("kp", 0.0f);
                const float kd = s.value("kd", 0.0f);
                const float max_torque = s.value("max_torque", 0.0f);

                AngularServoId angular_servo_id;

                angular_servo_id = add_angular_servo(
                    p1, p2, p3, p4,
                    target_angle,
                    kp,
                    kd,
                    max_torque);

                Entity *entity = entity_storage.get(entity_id);
                if (!entity) throw std::runtime_error("!entity");
                entity->angular_servos.push_back(angular_servo_id);
            }
        }
    }
    catch (std::runtime_error &e)
    {
        std::println("fail spawn_entity, file_name: {}\neror: {}", file_name, e.what());
        entity_storage.get_name(entity_id);
    }
    return entity_id;
}
} // namespace xpbd
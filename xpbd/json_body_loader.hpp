#pragma once
#include "core/cluster.hpp"

#include "glm/glm.hpp"
#include <string>

namespace json_body_loader
{
    extern std::string folder;
    extern const std::string extension;

    xpbd::cluster::Cluster* load(const std::string &filename, glm::vec2 offset = {0.0f, 0.0f});
}
#pragma once
#include "core/cluster.hpp"

#include "glm/glm.hpp"
#include <string>

namespace json_body_loader
{
    extern std::string folder;
    extern const std::string extension;

    bool load(xpbd::Cluster &cluster, const std::string &filename, glm::vec2 offset);
}
#ifndef JSON_BODY_LOADER
#define JSON_BODY_LOADER
#include <string>
#include "glm/glm.hpp"

namespace xpbd
{
    struct Cluster;
}

namespace json_body_loader
{
    extern std::string folder;
    extern const std::string extension;

    bool load(xpbd::Cluster &cluster, const std::string &filename, glm::vec2 offset);
}
#endif // JSON_BODY_LOADER
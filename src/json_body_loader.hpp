#ifndef JSON_BODY_LOADER
#define JSON_BODY_LOADER
#include <string>
#include "glm/glm.hpp"

namespace xpbd
{
    class World;
}

namespace json_body_loader
{
    extern std::string folder;
    extern const std::string extension;

    void load(xpbd::World &world, const std::string &filename, glm::vec2 offset);
}
#endif // JSON_BODY_LOADER
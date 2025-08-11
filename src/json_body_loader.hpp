#ifndef define_json_body_loader
#define define_json_body_loader
#include "glm/glm.hpp"
#include <string>

namespace xpbd
{
    class World;
}

namespace json_body_loader
{
    const std::string folder = "./body/";
    const std::string extension = ".body";

    void load(xpbd::World &world, const std::string &filename, glm::vec2 offset);
}
#endif // define_json_body_loader
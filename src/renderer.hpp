#include <SFML/Graphics.hpp>
#include "glm/glm.hpp"
#include <vector>

#include "imgui.h"
#include <imgui-SFML.h>

namespace renderer
{
    const int WINDOW_WIDTH = 2000;
    const int WINDOW_HEIGHT = 2000;

    extern sf::RenderWindow window;
    extern sf::View view;

    extern sf::Color color;

    void set_color(sf::Color color);
    void draw_segmented_line(const std::vector<glm::vec2> &points);
    void draw_line(const glm::vec2 &from, const glm::vec2 &to);
    void draw_circle(const glm::vec2 &pos, float radius);
    void draw_axis_aligned_bounding_box(float left, float right, float bottom, float top);

    void setup_window();
    void setup_view();
    bool setup_imgui();
};

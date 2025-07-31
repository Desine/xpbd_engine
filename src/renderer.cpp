#include "renderer.hpp"

namespace renderer
{

    sf::RenderWindow window;
    sf::View view;
    sf::Color color;
    

    void setup_window()
    {
        window.create(sf::VideoMode({WINDOW_WIDTH, WINDOW_HEIGHT}), "Soft Racing");
        window.setFramerateLimit(60);
    }

    void setup_view()
    {
        view = window.getDefaultView();
        view.setSize(WINDOW_WIDTH, -WINDOW_HEIGHT);
        view.zoom(0.5f);
        sf::Vector2f cameraCenter;
        cameraCenter.x = 0;
        cameraCenter.y = 0;
        view.setCenter(cameraCenter);
        window.setView(view);
    }

    bool setup_imgui()
    {
        if (!ImGui::SFML::Init(window))
            return false;

        ImGui::GetStyle().ScaleAllSizes(3.f);
        ImGui::GetIO().FontGlobalScale = 2.f;
        return true;
    }
    void set_color(sf::Color color)
    {
        renderer::color = color;
    }

    void draw_line(const glm::vec2 &from, const glm::vec2 &to)
    {
        sf::Vertex line[] =
            {
                sf::Vertex(sf::Vector2f(from.x, from.y), color),
                sf::Vertex(sf::Vector2f(to.x, to.y), color)};
        window.draw(line, 2, sf::Lines);
    }
    void draw_segmented_line(const std::vector<glm::vec2> &points)
    {
        for (int d = 1; d < points.size(); d++)
        {
            draw_line(points[d], points[d - 1]);
        }
    }
    void draw_circle(const glm::vec2 &pos, float radius)
    {
        sf::CircleShape shape(radius);
        shape.setFillColor(color);
        shape.setOrigin(radius, radius);
        shape.setPosition(pos.x, pos.y);
        window.draw(shape);
    }

    void draw_axis_aligned_bounding_box(float left, float right, float bottom, float top)
    {
        draw_line({left, bottom}, {right, bottom});
        draw_line({right, bottom}, {right, top});
        draw_line({right, top}, {left, top});
        draw_line({left, top}, {left, bottom});
    }
}
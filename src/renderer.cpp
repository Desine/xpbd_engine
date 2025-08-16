#include "renderer.hpp"

namespace renderer
{
    sf::RenderWindow window;
    sf::View view;
    sf::Color color;

    void setup_window()
    {
        window.create(sf::VideoMode({static_cast<unsigned int>(width), static_cast<unsigned int>(height)}), "xpbd_engine");
        window.setFramerateLimit(144);
    }

    void setup_view()
    {
        view = window.getDefaultView();
        view.setSize(width, -height);
        sf::Vector2f cameraCenter(0, 0);
        view.setCenter(cameraCenter);
        window.setView(view);
    }
    
    void zoom(float factor){
        view.zoom(factor);
        window.setView(view);
    }
    void resize(float w, float h)
    {
        width = w;
        height = h;
        view.setSize(width, -height);
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
            draw_line(points[d], points[d - 1]);
    }
    void draw_segmented_loop(const std::vector<glm::vec2> &points)
    {
        std::vector<glm::vec2> loop = points;
        loop.push_back(points[0]);
        draw_segmented_line(loop);
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
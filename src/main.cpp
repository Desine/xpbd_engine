#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>

#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "renderer.hpp"
#include "xpbd.hpp"
#include "json_body_loader.hpp"

#include "Remotery.h"

void draw_world(xpbd::World &world)
{
    renderer::set_color(sf::Color::Red);
    for (auto p : world.particles.pos)
        renderer::draw_circle(p, 3);

    for (size_t i = 0; i < world.distanceConstraints.i1.size(); ++i)
        renderer::draw_line(world.particles.pos[world.distanceConstraints.i1[i]], world.particles.pos[world.distanceConstraints.i2[i]]);

    renderer::set_color(sf::Color::Green);
    for (size_t i = 0; i < world.polygonColliders.indices.size(); ++i)
    {
        std::vector<glm::vec2> points;

        for (auto id : world.polygonColliders.indices[i])
            points.push_back(world.particles.pos[id]);

        renderer::draw_segmented_loop(points);
    }
}

void draw_collisions(std::vector<xpbd::PointPolygonCollision> collisions)
{
    for (auto c : collisions)
        renderer::draw_axis_aligned_bounding_box(c.box.l, c.box.r, c.box.b, c.box.t);
}

int main()
{
    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();

    rmtSettings *settings = rmt_Settings();
    settings->port = 17815; // default 17815
    Remotery *rmt;
    // rmt_CreateGlobalInstance(&rmt);

    xpbd::World world;
    world.init();

    sf::Vector2i mouseDrag;
    sf::Clock clock;
    while (renderer::window.isOpen())
    {
        sf::Time deltaTime = clock.restart();

        sf::Event event;
        while (renderer::window.pollEvent(event))
        {
            ImGui::SFML::ProcessEvent(renderer::window, event);

            if (event.type == sf::Event::Closed)
                renderer::window.close();

            // camera
            sf::Vector2i mouseCurrent = sf::Mouse::getPosition(renderer::window);

            const sf::Vector2f delta =
                renderer::window.mapPixelToCoords(mouseCurrent) -
                renderer::window.mapPixelToCoords(mouseDrag);
            mouseDrag = mouseCurrent;
            if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Right))
            {
                renderer::view.move(-delta);
                renderer::window.setView(renderer::view);
            }
            if (event.type == sf::Event::MouseWheelScrolled && event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel)
            {
                if (event.mouseWheelScroll.delta > 0)
                    renderer::view.zoom(0.95);
                else
                    renderer::view.zoom(1.05);
                renderer::window.setView(renderer::view);
            }

            // spawn
            if (event.type == sf::Event::KeyReleased)
            {
                sf::Vector2f worldPos = renderer::window.mapPixelToCoords(mouseCurrent);
                glm::vec2 position(worldPos.x, worldPos.y);
                if (event.key.code == sf::Keyboard::A)
                    world.spawnFromJson("2boxes", position);
                if (event.key.code == sf::Keyboard::S)
                    world.spawnFromJson("square", position);
                if (event.key.code == sf::Keyboard::D)
                    world.spawnFromJson("balloon", position);
                if (event.key.code == sf::Keyboard::F)
                    world.addPolygon(position, 40, 6, 3, 0.005f);

                if (event.key.code == sf::Keyboard::R)
                    world.init();
            }
        }

        ImGui::NewFrame();
        ImGui::Begin("Main");
        ImGui::Text("FPS: %.1f", 1.0f / deltaTime.asSeconds());
        if (ImGui::Button(world.paused ? "Play" : "Pause"))
            world.paused = !world.paused;
        if (ImGui::Button("stepOnce - todo"))
            world.stepOnce = true;
        ImGui::SliderFloat("timeScale", &world.timeScale, 0.01f, 50, "%.3f");
        size_t min_value = 1;
        size_t max_value = 20;
        ImGui::SliderScalar("substeps", ImGuiDataType_U64, &world.substeps, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderScalar("iterations", ImGuiDataType_U64, &world.iterations, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderFloat("gravity.x", &world.gravity.x, -20, 20);
        ImGui::SliderFloat("gravity.y", &world.gravity.y, -20, 20);
        ImGui::End();
        ImGui::EndFrame();

        renderer::window.clear();
        draw_world(world);
        renderer::set_color(sf::Color::White);
        draw_collisions(world.collisions);
        ImGui::SFML::Render(renderer::window);
        renderer::window.display();

        world.update(deltaTime.asSeconds());
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();
    // rmt_DestroyGlobalInstance(rmt);
}
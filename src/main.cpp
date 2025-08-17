#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>

#include <stdio.h>
#include <string>
#include <omp.h>
#include <thread>

#include "renderer.hpp"
#include "xpbd.hpp"

#include "Remotery.h"

void draw_world(xpbd::World &world)
{
    renderer::set_color(sf::Color::Red);
    for (auto p : world.particles.pos)
        renderer::draw_circle(p, 2);

    for (size_t i = 0; i < world.distanceConstraints.size(); ++i)
        renderer::draw_line(world.particles.pos[world.distanceConstraints[i].i1], world.particles.pos[world.distanceConstraints[i].i2]);

    renderer::set_color(sf::Color::Green);
    for (auto &pc : world.polygonColliders)
    {
        std::vector<glm::vec2> points;
        points.reserve(pc.indices.size());

        for (auto id : pc.indices)
            points.emplace_back(world.particles.pos[id]);

        renderer::draw_segmented_loop(points);
    }

    renderer::set_color(sf::Color::Green);
    for (auto &pc : world.pointsColliders)
        for (auto id : pc.indices)
            renderer::draw_circle(world.particles.pos[id], 2);
}

void draw_collisions(std::vector<xpbd::PointsPolygonCollision> collisions)
{
    for (auto c : collisions)
        renderer::draw_axis_aligned_bounding_box(c.box.l, c.box.r, c.box.b, c.box.t);
}

void ShowFpsGraph(float newFps)
{
    const size_t historyFrames = 120;
    static float fpsHistory[historyFrames] = {0};
    static int offset = 0;
    offset = (offset + 1) % IM_ARRAYSIZE(fpsHistory);

    fpsHistory[offset] = newFps;

    float fpsHistoryToDraw[historyFrames];
    std::copy(std::begin(fpsHistory), std::end(fpsHistory), std::begin(fpsHistoryToDraw));

    float min = 0.0f;
    float max = 144.0f;

    static bool scroll = true;
    float scrollFloat = offset;
    if (!scroll)
    {
        scrollFloat = 0.0f;
        fpsHistoryToDraw[offset] = max;
        fpsHistoryToDraw[(offset + 1) % IM_ARRAYSIZE(fpsHistory)] = min;
    }

    ImVec2 size(0, 150);

    ImGui::PlotLines(
        "FPS",
        fpsHistoryToDraw,
        IM_ARRAYSIZE(fpsHistory),
        scrollFloat,
        nullptr,
        min,
        max,
        size);
    if (ImGui::Button("toggle fps scroll"))
        scroll = !scroll;
}

void spawn_many(
    xpbd::World &world,
    glm::vec2 &pos,
    size_t amountX,
    size_t amountY,
    float spacingX,
    float spacingY,
    std::string name)
{

    for (size_t i = 0; i < amountX; ++i)
    {
        float x = pos.x - amountX * spacingX * 0.5 + spacingX * i;
        for (size_t j = 0; j < amountY; ++j)
        {
            float y = pos.y + spacingY * j;
            world.spawnFromJson(name, {x, y});
        }
    }
}

int main()
{
    rmtSettings *settings = rmt_Settings();
    settings->port = 17815; // default 17815
    Remotery *rmt;
    // rmt_CreateGlobalInstance(&rmt);

    xpbd::World world;
    world.init();
    // world.spawnPolygon({0, 0}, 40, 6, 3, 0.005f);

    // while (true)
    // {
    //     world.update(0.01f);
    // }

    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();

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

            if (event.type == sf::Event::Resized)
                renderer::resize(
                    static_cast<float>(event.size.width),
                    static_cast<float>(event.size.height));

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
                renderer::zoom(event.mouseWheelScroll.delta > 0 ? 0.95 : 1.05);

            // spawn
            if (event.type == sf::Event::KeyReleased)
            {
                sf::Vector2f worldPos = renderer::window.mapPixelToCoords(mouseCurrent);
                glm::vec2 position(worldPos.x, worldPos.y);
                if (event.key.code == sf::Keyboard::Q)
                {
                    size_t amountX = 5;
                    size_t amountY = 5;
                    float spacingX = 300;
                    float spacingY = 200;
                    spawn_many(world, position, amountX, amountY, spacingX, spacingY, "car");
                }
                if (event.key.code == sf::Keyboard::W)
                {
                    size_t amountX = 5;
                    size_t amountY = 5;
                    float spacingX = 300;
                    float spacingY = 200;
                    spawn_many(world, position, amountX, amountY, spacingX, spacingY, "balloon");
                }
                if (event.key.code == sf::Keyboard::E)
                {
                    size_t amountX = 10;
                    size_t amountY = 10;
                    float spacingX = 100;
                    float spacingY = 100;
                    for (size_t i = 0; i < amountX; ++i)
                    {
                        float x = position.x - amountX * spacingX * 0.5 + spacingX * i;
                        for (size_t j = 0; j < amountY; ++j)
                        {
                            float y = position.y + spacingY * j;
                            const float radius = 40;
                            const size_t segments = 6;
                            const float mass = 10;
                            const float compliance = 0.005f;
                            world.spawnPolygon({x, y}, radius, segments, mass, compliance);
                        }
                    }
                }

                if (event.key.code == sf::Keyboard::A)
                    world.spawnFromJson("car", position);
                if (event.key.code == sf::Keyboard::S)
                    world.spawnFromJson("balloon", position);
                if (event.key.code == sf::Keyboard::D)
                {
                    const float radius = 40;
                    const size_t segments = 6;
                    const float mass = 5;
                    const float compliance = 0.005f;
                    world.spawnPolygon(position, radius, segments, mass, compliance);
                }
                if (event.key.code == sf::Keyboard::F)
                    world.spawnFromJson("person", position);
                if (event.key.code == sf::Keyboard::G)
                    world.spawnFromJson("square", position);

                if (event.key.code == sf::Keyboard::R)
                    world.init();
            }
        }

        ImGui::NewFrame();
        ImGui::Begin("Main");
        ShowFpsGraph(1.0f / deltaTime.asSeconds());
        if (ImGui::Button(world.paused ? "Play" : "Pause"))
            world.paused = !world.paused;
        // if (ImGui::Button("stepOnce - todo"))
        //     world.stepOnce = true;
        ImGui::SliderFloat("timeScale", &world.timeScale, 0.01f, 50, "%.3f");
        size_t min_value = 1;
        size_t max_value = 20;
        ImGui::SliderScalar("substeps", ImGuiDataType_U64, &world.substeps, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderScalar("iterations", ImGuiDataType_U64, &world.iterations, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderFloat("gravity.x", &world.gravity.x, -20, 20);
        ImGui::SliderFloat("gravity.y", &world.gravity.y, -20, 20);
        // int cellSize = world.spatialHashAABB.cellSize;
        // ImGui::SliderInt("spatialHashing_cellSize: ", &cellSize, 30, 1000);
        // world.spatialHashAABB.cellSize = cellSize;

        ImGui::End();
        ImGui::EndFrame();

        renderer::window.clear();
        draw_world(world);
        // renderer::set_color(sf::Color::White);
        // draw_collisions(world.collisions);

        world.update(deltaTime.asSeconds());

        ImGui::SFML::Render(renderer::window);
        renderer::window.display();
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();
    // rmt_DestroyGlobalInstance(rmt);
}
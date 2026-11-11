#include <filesystem>
#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>
#include "glm/glm.hpp"

#include <print>

#include "renderer.hpp"
#include "xpbd.hpp"

#include "profiler.hpp"

void draw_worldDraw(xpbd::WorldDraw worldDraw, glm::vec2 offset = {0.0f, 0.0f})
{
    renderer::set_color(sf::Color::Red);
    for (auto p : worldDraw.particles.pos)
        renderer::draw_circle(p, 2);

    for (size_t i = 0; i < worldDraw.distanceConstraints.size(); ++i)
        renderer::draw_line(worldDraw.particles.pos[worldDraw.distanceConstraints[i].i1] + offset, worldDraw.particles.pos[worldDraw.distanceConstraints[i].i2] + offset);

    renderer::set_color(sf::Color::Green);
    for (auto &pc : worldDraw.polygonColliders)
    {
        std::vector<glm::vec2> points;
        points.reserve(pc.indices.size());

        for (auto id : pc.indices)
            points.emplace_back(worldDraw.particles.pos[id] + offset);

        renderer::draw_segmented_loop(points);
    }

    renderer::set_color(sf::Color::Green);
    for (auto &pc : worldDraw.pointsColliders)
        for (auto id : pc.indices)
            renderer::draw_circle(worldDraw.particles.pos[id], 2);
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
            world.enqueue_spawn_from_json(name, {x, y});
        }
    }
}

int main()
{
    //json_body_loader::folder = std::filesystem::current_path().u8string() + "./body/";


    xpbd::World world;
    world.init();
    world.spawn_from_json("ground", {0, -1000});
    world.submit_draw();
    bool paused = world.updateState == xpbd::World::UpdateState::Paused;

    xpbd::WorldDraw worldDraw;
    // world.spawnPolygon({0, 0}, 40, 6, 3, 0.005f);

    // while (true)
    // {
    //     world.step(0.01f);
    // }

    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();

    sf::Vector2i mouseDrag;
    sf::Clock clock;
    while (renderer::window.isOpen())
    {
        FrameMark;
        sf::Time deltaTime = clock.restart();

        while (const auto event = renderer::window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(renderer::window, *event);

            if (event->is<sf::Event::Closed>())
                renderer::window.close();

            if (const auto resized = event->getIf<sf::Event::Resized>())
                renderer::resize(
                    resized->size.x,
                    resized->size.y);

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
            if (const auto scrolled = event->getIf<sf::Event::MouseWheelScrolled>())
                renderer::zoom(scrolled->delta > 0 ? 0.95 : 1.05);

            // keyboard
            if (const auto keyReleased = event->getIf<sf::Event::KeyReleased>())
            {
                // spawn
                sf::Vector2f worldPos = renderer::window.mapPixelToCoords(mouseCurrent);
                glm::vec2 position(worldPos.x, worldPos.y);
                if (keyReleased->code == sf::Keyboard::Key::Q)
                {
                    size_t amountX = 5;
                    size_t amountY = 5;
                    float spacingX = 300;
                    float spacingY = 200;
                    spawn_many(world, position, amountX, amountY, spacingX, spacingY, "car");
                }
                if (keyReleased->code == sf::Keyboard::Key::W)
                {
                    size_t amountX = 5;
                    size_t amountY = 5;
                    float spacingX = 300;
                    float spacingY = 200;
                    spawn_many(world, position, amountX, amountY, spacingX, spacingY, "balloon");
                }
                if (keyReleased->code == sf::Keyboard::Key::E)
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
                            world.enqueue_spawn_polygon({x, y}, radius, segments, mass, compliance);
                        }
                    }
                }
                if (keyReleased->code == sf::Keyboard::Key::C)
                {
                    size_t amountX = 30;
                    size_t amountY = 1;
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
                            world.enqueue_spawn_polygon({x, y}, radius, segments, mass, compliance);
                        }
                    }
                }

                if (keyReleased->code == sf::Keyboard::Key::A)
                    world.enqueue_spawn_from_json("car", position);
                if (keyReleased->code == sf::Keyboard::Key::S)
                    world.enqueue_spawn_from_json("balloon", position);
                if (keyReleased->code == sf::Keyboard::Key::D)
                {
                    float radius = 40;
                    size_t segments = 6;
                    float mass = 5;
                    float compliance = 0.005f;
                    world.enqueue_spawn_polygon(position, radius, segments, mass, compliance);
                }
                if (keyReleased->code == sf::Keyboard::Key::F)
                    world.enqueue_spawn_from_json("person", position);
                if (keyReleased->code == sf::Keyboard::Key::G)
                    world.enqueue_spawn_from_json("square", position);

                // clear world
                if (keyReleased->code == sf::Keyboard::Key::R)
                {
                    world.enqueue_init();
                    world.enqueue_spawn_from_json("ground", {0, -1000});
                }
                // pause / resume
                if (keyReleased->code == sf::Keyboard::Key::Space)
                {
                    if (paused)
                        world.resume_update();
                    else
                        world.pause_update();

                    paused = !paused;
                }
            }
        }

        ImGui::NewFrame();
        ImGui::Begin("Main");
        ShowFpsGraph(1.0f / deltaTime.asSeconds());
        static bool draw = true;
        if (ImGui::Button(draw ? "turn off: draw" : "turn on: draw"))
            draw = !draw;
        static bool draw_thick = true;
        if (ImGui::Button(draw_thick ? "turn off: draw_thick" : "turn on: draw_thick"))
            draw_thick = !draw_thick;

        if (ImGui::Button(paused ? "Resume" : "Pause"))
        {
            if (paused)
                world.resume_update();
            else
                world.pause_update();

            paused = !paused;
        }
        if (ImGui::Button("Step once"))
            world.step_once_update();
        // ImGui::SliderFloat("timeScale", &world.timeScale, 0.01f, 50, "%.3f");
        // size_t min_value = 1;
        // size_t max_value = 20;
        // ImGui::SliderScalar("substeps", ImGuiDataType_U64, &world.substeps, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        // ImGui::SliderScalar("iterations", ImGuiDataType_U64, &world.iterations, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        // ImGui::SliderFloat("gravity.x", &world.gravity.x, -20, 20);
        // ImGui::SliderFloat("gravity.y", &world.gravity.y, -20, 20);
        // int cellSize = world.spatialHashAABB.cellSize;
        // ImGui::SliderInt("spatialHashing_cellSize: ", &cellSize, 30, 1000);
        // world.spatialHashAABB.cellSize = cellSize;

        ImGui::End();
        ImGui::EndFrame();

        renderer::window.clear();
        if (draw)
        {
            ZoneScoped;
            {
                std::lock_guard<std::mutex> lk(world.drawMutex);
                if (world.newFrame)
                {
                    std::swap(worldDraw, world.draw);
                    world.newFrame = false;
                }
            }
            draw_worldDraw(worldDraw);
            if (draw_thick)
            {
                draw_worldDraw(worldDraw, {0.0f, 0.5f});
                draw_worldDraw(worldDraw, {0.0f, 1.0f});
                draw_worldDraw(worldDraw, {0.5f, 0.0f});
                draw_worldDraw(worldDraw, {1.0f, 0.0f});
                draw_worldDraw(worldDraw, {1.0f, 0.5f});
                draw_worldDraw(worldDraw, {0.5f, 1.0f});
                draw_worldDraw(worldDraw, {1.0f, 1.0f});
            }
        }
        // renderer::set_color(sf::Color::White);
        // draw_collisions(world.collisions);

        ImGui::SFML::Render(renderer::window);
        renderer::window.display();
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();
}
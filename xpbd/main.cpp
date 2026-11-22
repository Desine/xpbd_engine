#include "xpbd.hpp"
#include "renderer.hpp"
#include "json_body_loader.hpp"
#include "utils.hpp"

#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>

#include "glm/glm.hpp"
#include <vector>
#include "tracy/Tracy.hpp"


void draw_worldDraw(xpbd::Cluster clusterDraw, glm::vec2 offset = {0.0f, 0.0f})
{
    float circleSize = 0.03f;
    renderer::set_color(sf::Color::Red);
    for (auto p : clusterDraw.particles)
        renderer::draw_circle(p.pos, circleSize);

    for (size_t i = 0; i < clusterDraw.distanceConstraints.size(); ++i)
        renderer::draw_line(clusterDraw.particles[clusterDraw.distanceConstraints[i].i1].pos + offset, clusterDraw.particles[clusterDraw.distanceConstraints[i].i2].pos + offset);

    renderer::set_color(sf::Color::Green);
    for (auto &pc : clusterDraw.polygonColliders)
    {
        std::vector<glm::vec2> points;
        points.reserve(pc.indices.size());

        for (auto id : pc.indices)
            points.emplace_back(clusterDraw.particles[id].pos + offset);

        renderer::draw_segmented_loop(points);
    }

    renderer::set_color(sf::Color::Green);
    for (auto &pc : clusterDraw.pointsColliders)
        for (auto id : pc.indices)
            renderer::draw_circle(clusterDraw.particles[id].pos, circleSize);
}

void draw_collisions(std::vector<xpbd::PointsPolygonCollision> collisions)
{
    for (auto c : collisions)
        renderer::draw_axis_aligned_bounding_box(c.box.l, c.box.r, c.box.b, c.box.t);
}

void ShowFpsGraph(float newFps)
{
    const size_t n = 120;
    static float hist[n] = {0};
    static int writePos = 0;
    static int displayOffset = 0;
    static bool scroll = true;

    writePos = (writePos + 1) % n;
    hist[writePos] = newFps;

    if (scroll)
        displayOffset = (writePos + 1) % n;

    ImGui::PlotLines("FPS", hist, (int)n, displayOffset,
                     nullptr, 0.0f, 60.0f, ImVec2(0, 150));

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
            world.enqueue_command(xpbd::spawn_from_json_command(name, {x, y}));
        }
    }
}

int main()
{
    json_body_loader::folder = std::filesystem::current_path().generic_string() + "/body/";
    // std::println("folder {}", json_body_loader::folder);

    xpbd::World world;
    world.enqueue_command(xpbd::init_command());
    world.enqueue_command(xpbd::spawn_from_json_command("ground", {0, -10}));
    bool paused = world.updateState == xpbd::World::UpdateState::Paused;
    bool reversePass = world.reversePass;
    int drawThickness = 1;
    int targetRPM;


    glm::vec2 spawn_pos(0.0f, -5.0f);
    spawn_many(world, spawn_pos, 5, 2, 3, 2, "car");
    spawn_pos = {0.0f, 0.0f};
    spawn_many(world, spawn_pos, 5, 10, 3, 2, "balloon");



    xpbd::Cluster clusterDraw;

    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();
    renderer::zoom(0.01f);

    sf::Vector2i mouseDrag;
    sf::Clock clock;
    while (renderer::window.isOpen())
    {
        // FrameMark;
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
                    size_t amountY = 1;
                    float spacingX = 3;
                    float spacingY = 2;
                    spawn_many(world, position, amountX, amountY, spacingX, spacingY, "car");
                }
                if (keyReleased->code == sf::Keyboard::Key::W)
                {
                    size_t amountX = 5;
                    size_t amountY = 5;
                    float spacingX = 3;
                    float spacingY = 2;
                    spawn_many(world, position, amountX, amountY, spacingX, spacingY, "balloon");
                }
                if (keyReleased->code == sf::Keyboard::Key::Z)
                    world.enqueue_command(xpbd::spawn_from_json_command("wheel", position));
                if (keyReleased->code == sf::Keyboard::Key::C)
                {
                    size_t amountX = 10;
                    size_t amountY = 1;
                    float spacingX = 0.5f;
                    float spacingY = 1;
                    for (size_t i = 0; i < amountX; ++i)
                    {
                        float x = position.x - amountX * spacingX * 0.5 + spacingX * i;
                        for (size_t j = 0; j < amountY; ++j)
                        {
                            float y = position.y + spacingY * j;
                            const float radius = 0.2f;
                            const size_t segments = 6;
                            const float mass = 10;
                            const float compliance = 0.0f;
                            world.enqueue_command(xpbd::spawn_polygon_command({x, y}, radius, segments, mass, compliance));
                        }
                    }
                }

                if (keyReleased->code == sf::Keyboard::Key::A)
                    world.enqueue_command(xpbd::spawn_from_json_command("car", position));
                if (keyReleased->code == sf::Keyboard::Key::S)
                    world.enqueue_command(xpbd::spawn_from_json_command("balloon", position));
                if (keyReleased->code == sf::Keyboard::Key::D)
                {
                    float radius = 0.2f;
                    size_t segments = 6;
                    float mass = 5;
                    float compliance = 0.0f;
                    world.enqueue_command(xpbd::spawn_polygon_command(position, radius, segments, mass, compliance));
                }
                if (keyReleased->code == sf::Keyboard::Key::F)
                    world.enqueue_command(xpbd::spawn_from_json_command("person", position));
                if (keyReleased->code == sf::Keyboard::Key::G)
                    world.enqueue_command(xpbd::spawn_from_json_command("square", position));

                // clear world
                if (keyReleased->code == sf::Keyboard::Key::R)
                {
                    world.enqueue_command(xpbd::init_command());
                    world.enqueue_command(xpbd::spawn_from_json_command("ground", {0, -10}));
                }
                // pause / resume
                if (keyReleased->code == sf::Keyboard::Key::Space)
                {
                    if (world.updateState == xpbd::World::UpdateState::Paused || world.updateState == xpbd::World::UpdateState::StepOnce)
                        world.resume_update();
                    else
                        world.pause_update();
                    paused = world.updateState == xpbd::World::UpdateState::Paused;
                }
            }
        }

        ImGui::NewFrame();
        ImGui::Begin("Main");

        float dt = world.lastStepTime.count();
        float fps = (dt > 0.0f) ? (1.0f / dt) : 0.0f;
        ShowFpsGraph(fps);

        static bool draw = true;
        if (ImGui::Button(draw ? "turn off: draw" : "turn on: draw"))
            draw = !draw;
        static bool draw_thick = true;
        if (ImGui::Button(draw_thick ? "turn off: draw_thick" : "turn on: draw_thick"))
            draw_thick = !draw_thick;

        if (ImGui::Button(paused ? "Resume" : "Pause"))
        {
            paused ? world.resume_update() : world.pause_update();
            paused = !paused;
        }
        if (ImGui::Button(reversePass ? "Reverse pass is on" : "Reverse pass is off"))
        {
            reversePass ? world.reversePass = false : world.reversePass = true;
            reversePass = !reversePass;
        }
        if (ImGui::Button("Step once"))
            world.step_once_update();
        // ImGui::SliderFloat("gravity.x", &world.gravity.x, -20, 20);
        // ImGui::SliderFloat("gravity.y", &world.gravity.y, -20, 20);
        ImGui::SliderInt("Draw thickness", &drawThickness, 1, 10);
        if (ImGui::SliderInt("targetRPM", &targetRPM, -1000, 1000))
        {
            if (world.cluster.angularMotors.size() > 1)
            {
                xpbd::update_angular_motor_command motor1;
                motor1.index = 0;
                motor1.motor = world.cluster.angularMotors[0];
                motor1.motor.targetRPM = targetRPM;
                world.enqueue_command(motor1);

                xpbd::update_angular_motor_command motor2;
                motor2.index = 1;
                motor2.motor = world.cluster.angularMotors[1];
                motor2.motor.targetRPM = targetRPM;
                world.enqueue_command(motor2);
            }
        }

        ImGui::End();
        ImGui::EndFrame();

        renderer::window.clear();
        if (draw)
        {
            ZoneScoped;
            if (world.newFrame)
            {
                world.newFrame = false;
                std::lock_guard<std::mutex> lk(world.clusterDrawMutex);
                std::swap(clusterDraw, world.clusterDraw);
            }
            float angleStep = 2.0f * float(M_PI) / drawThickness;
            for (size_t i = 0; i < size_t(drawThickness); ++i)
            {
                float angle = i * angleStep;
                glm::vec2 dir = glm::vec2(cosf(angle), sinf(angle));

                draw_worldDraw(clusterDraw, dir * 0.005f);
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
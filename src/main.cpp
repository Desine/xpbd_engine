#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>

#include <stdio.h>
#include <string.h>

#include "renderer.hpp"
#include "xpbd.hpp"
#include "json_body_loader.hpp"

#include "Remotery.h"

void add_poligon(xpbd::Particles &p, xpbd::DistanceConstraints &dc, xpbd::VolumeConstraints &vc, xpbd::ColliderPoints &pc, glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
{
    if (segments < 3)
        segments = 3;

    size_t start = p.pos.size();
    size_t end = start + segments;
    std::vector<size_t> ids;

    float angleStep = 2.0f * float(M_PI) / segments;
    for (int i = 0; i < segments; ++i)
    {
        float angle = i * angleStep;
        glm::vec2 dir = glm::vec2(cosf(angle), sinf(angle));
        xpbd::add_particle(p, dir * radius + pos, mass / segments);

        ids.push_back(start + i);
    }

    for (int i = 0; i < segments; ++i)
        xpbd::add_distance_constraint_auto_restDist(dc, start + i, start + (i + 1) % segments, compliance, p);

    xpbd::add_volume_constraint(p, vc, ids, compliance);

    xpbd::add_polygon_collider(pc, ids, 0.4f, 0.3f, compliance);
}
void draw_aabbs(std::vector<xpbd::AABB> aabb)
{
    for (auto a : aabb)
        renderer::draw_axis_aligned_bounding_box(a.l, a.r, a.b, a.t);
}
void draw_aabbs_intersections(std::vector<xpbd::AABB> aabb, std::vector<xpbd::AABBsOverlap> intersections)
{
    for (auto i : intersections)
    {
        renderer::draw_axis_aligned_bounding_box(aabb[i.i1].l, aabb[i.i1].r, aabb[i.i1].b, aabb[i.i1].t);
        renderer::draw_axis_aligned_bounding_box(aabb[i.i2].l, aabb[i.i2].r, aabb[i.i2].b, aabb[i.i2].t);
    }
}
void draw_point_edge_collision_constraints(xpbd::Particles p, xpbd::PointEdgeCollisionConstraints pecc)
{
    for (size_t i = 0; i < pecc.point.size(); ++i)
    {
        renderer::set_color(sf::Color::Red);
        renderer::draw_circle(p.pos[pecc.point[i]], 5);
        renderer::set_color(sf::Color::Magenta);
        renderer::draw_line(p.pos[pecc.edge1[i]], p.pos[pecc.edge2[i]]);
    }
}

int main()
{
    rmtSettings *settings = rmt_Settings();
    settings->port = 17816;
    Remotery *rmt;
    rmt_CreateGlobalInstance(&rmt);

    size_t substeps = 10;
    size_t iterations = 1;
    float sec = 0.0f;
    float rate = 60.0f;
    float deltaTick = 1.0f / rate;
    float timeScale = 10;
    bool paused = true;
    bool stepOnce = false;
    glm::vec2 gravity = {0, -9.8f};
    xpbd::Particles particles;
    xpbd::DistanceConstraints distanceConstraints;
    xpbd::VolumeConstraints volumeConstraints;
    xpbd::ColliderPoints polygonColliders;
    xpbd::ColliderPoints pointColliders;

    json_body_loader::load("2boxes", particles, distanceConstraints, volumeConstraints, polygonColliders, pointColliders);

    add_poligon(particles, distanceConstraints, volumeConstraints, polygonColliders, {300, 300}, 70, 5, 5, 0.001f);

    json_body_loader::load("ground", particles, distanceConstraints, volumeConstraints, polygonColliders, pointColliders);
    json_body_loader::load("square", particles, distanceConstraints, volumeConstraints, polygonColliders, pointColliders, {0, 400});

    json_body_loader::load("balloon", particles, distanceConstraints, volumeConstraints, polygonColliders, pointColliders, {-150, 700});

    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();

    sf::Vector2i mouseDrag;

    sf::Clock clock;
    while (renderer::window.isOpen())
    {
        sf::Time deltaTime = clock.restart();
        if (!paused)
            sec += deltaTime.asSeconds();

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
            if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left)
            {
                sf::Vector2i pixelPos = sf::Mouse::getPosition(renderer::window);
                sf::Vector2f worldPos = renderer::window.mapPixelToCoords(pixelPos);
                glm::vec2 position(worldPos.x, worldPos.y);
                // json_body_loader::load("square", particles, distanceConstraints, volumeConstraints, polygonColliders, pointColliders, position);
                add_poligon(particles, distanceConstraints, volumeConstraints, polygonColliders, position, 40, 6, 3, 0.005f);
            }
        }

        renderer::window.clear();

        ImGui::NewFrame();
        ImGui::Begin("Main");
        ImGui::Text("FPS: %.1f", 1.0f / deltaTime.asSeconds());
        if (ImGui::Button(paused ? "Play" : "Pause"))
            paused = !paused;
        if (ImGui::Button("stepOnce - todo"))
            stepOnce = true;
        ImGui::SliderFloat("timeScale", &timeScale, 0.01f, 50, "%.3f");
        size_t min_value = 1;
        size_t max_value = 10;
        ImGui::SliderScalar("substeps", ImGuiDataType_U64, &substeps, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderScalar("iterations", ImGuiDataType_U64, &iterations, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderFloat("gravity.x", &gravity.x, -20, 20);
        ImGui::SliderFloat("gravity.y", &gravity.y, -20, 20);

        ImGui::End();
        ImGui::EndFrame();

        float substep_time = deltaTick * timeScale / substeps;
        while (!paused && xpbd::should_tick(sec, deltaTick))
        {
            for (size_t s = 0; s < substeps; s++)
            {
                xpbd::iterate(particles, substep_time, gravity);

                xpbd::reset_constraints_lambdas(distanceConstraints.lambda);
                xpbd::reset_constraints_lambdas(volumeConstraints.lambda);

                xpbd::PointEdgeCollisionConstraints pecc;
                for (size_t i = 0; i < iterations; i++)
                {
                    xpbd::solve_distance_constraints(particles, distanceConstraints, substep_time);
                    xpbd::solve_volume_constraints(particles, volumeConstraints, substep_time);

                    // polygon/polygon collision detection
                    std::vector<xpbd::AABB> aabbs_polygons = xpbd::generate_particles_aabbs(particles, polygonColliders.indices);
                    std::vector<xpbd::AABBsOverlap> aabbs_polygon_polygon_intersections = xpbd::create_aabbs_intersections(aabbs_polygons);
                    for (auto a : aabbs_polygon_polygon_intersections)
                        xpbd::add_point_edge_collision_constraints_of_polygon_to_polygon_colliders(particles, pecc, polygonColliders, a);

                    // point/polygon collisions detection
                    std::vector<xpbd::AABB> aabbs_points = xpbd::generate_particles_aabbs(particles, pointColliders.indices);
                    std::vector<xpbd::AABBsOverlap> aabbs_point_polygon_intersections = xpbd::create_aabbs_intersections(aabbs_points, aabbs_polygons);
                    for (auto a : aabbs_point_polygon_intersections)
                        xpbd::add_point_edge_collision_constraints_of_point_to_polygon_colliders(particles, pecc, pointColliders, polygonColliders, a);

                    xpbd::solve_point_edge_collision_constraints(particles, pecc, substep_time);
                }

                xpbd::update_velocities(particles, substep_time);
                // xpbd::apply_point_edge_collision_constraints_kinetic_friction(particles, pecc, substep_time);
                // xpbd::apply_distance_constraints_damping(particles, distanceConstraints, substep_time);
            }
        }

        renderer::set_color(sf::Color::Green);
        for (auto p : particles.pos)
            renderer::draw_circle(p, 3);

        for (size_t i = 0; i < distanceConstraints.i1.size(); ++i)
            renderer::draw_line(particles.pos[distanceConstraints.i1[i]], particles.pos[distanceConstraints.i2[i]]);

        renderer::set_color(sf::Color::Red);
        for (size_t i = 0; i < polygonColliders.indices.size(); ++i)
        {
            std::vector<glm::vec2> points;

            for (auto id : polygonColliders.indices[i])
                points.push_back(particles.pos[id]);

            renderer::draw_segmented_loop(points);
        }

        ImGui::SFML::Render(renderer::window);
        renderer::window.display();
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();

    rmt_DestroyGlobalInstance(rmt);
}
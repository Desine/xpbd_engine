#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>

#include <stdio.h>
#include <string.h>

#include "renderer.hpp"
#include "xpbd.hpp"
#include "json_body_loader.hpp"

void add_poligon(xpbd::Particles &p, xpbd::DistanceConstraints &dc, xpbd::VolumeConstraints &vc, glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
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
}

int main()
{
    size_t substeps = 10;
    size_t iterations = 1;
    float sec = 0.0f;
    float rate = 60.0f;
    float deltaTick = 1.0f / rate;
    float timeScale = 1;
    bool paused = true;
    bool stepOnce = false;
    glm::vec2 gravity = {0, -9.8f};
    gravity = {0, 0};
    xpbd::Particles particles;
    xpbd::DistanceConstraints distanceConstraints;
    xpbd::VolumeConstraints volumeConstraints;
    xpbd::ContourColliders contourColliders;

    json_body_loader::load("2boxes", particles, distanceConstraints, volumeConstraints, contourColliders);

    add_poligon(particles, distanceConstraints, volumeConstraints, {300, 300}, 70, 5, 5, 0.001f);
    xpbd::add_distance_constraint(distanceConstraints, 7, 8, 0.1f, 0);
    xpbd::add_contour_collider(contourColliders, {8, 9, 10, 11, 12}, 1, 0.5f, 0);

    json_body_loader::load("square", particles, distanceConstraints, volumeConstraints, contourColliders);

    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();

    sf::Clock clock;
    while (renderer::window.isOpen())
    {
        sf::Time deltaTime = clock.restart();
        if (!paused)
            sec += deltaTime.asSeconds();

        sf::Event event;
        while (renderer::window.pollEvent(event))
        {
            ImGui::SFML::ProcessEvent(event);

            if (event.type == sf::Event::Closed)
                renderer::window.close();
        }
        renderer::window.clear();

        ImGui::NewFrame();
        ImGui::Begin("Main");
        if (ImGui::Button(paused ? "Play" : "Pause"))
            paused = !paused;
        if (ImGui::Button("stepOnce - todo"))
            stepOnce = true;
        ImGui::SliderFloat("timeScale", &timeScale, 0.001f, 5, "%.3f");
        size_t min_value = 1;
        size_t max_value = 10;
        ImGui::SliderScalar("substeps", ImGuiDataType_U64, &substeps, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderScalar("iterations", ImGuiDataType_U64, &iterations, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::End();
        ImGui::EndFrame();

        std::vector<xpbd::AABB> aabbs;
        std::vector<xpbd::AABBsIntersection> aabbs_intersections;

        float substep_time = deltaTick * timeScale / substeps;
        while (!paused && xpbd::should_tick(sec, deltaTick))
        {
            for (size_t s = 0; s < substeps; s++)
            {
                xpbd::iterate(particles, substep_time, gravity);
                xpbd::reset_constraints_lambdas(distanceConstraints.lambda);
                xpbd::reset_constraints_lambdas(volumeConstraints.lambda);
                for (size_t i = 0; i < iterations; i++)
                {
                    xpbd::solve_distance_constraints(particles, distanceConstraints, substep_time);
                    xpbd::solve_volume_constraints(particles, volumeConstraints, substep_time);
                    aabbs = xpbd::generate_colliders_aabbs(particles, contourColliders.indices);
                    aabbs_intersections = xpbd::find_aabbs_intersections(aabbs);

                    xpbd::PointEdgeCollisionConstraints pointEdgeCollisionConstraints;
                    renderer::set_color(sf::Color::Blue);
                    for (auto a : aabbs_intersections)
                    {
                        renderer::draw_axis_aligned_bounding_box(aabbs[a.i1].l, aabbs[a.i1].r, aabbs[a.i1].b, aabbs[a.i1].t);
                        renderer::draw_axis_aligned_bounding_box(aabbs[a.i2].l, aabbs[a.i2].r, aabbs[a.i2].b, aabbs[a.i2].t);
                        xpbd::add_point_edge_collision_constraints(particles, pointEdgeCollisionConstraints, contourColliders, a.i1, a.i2);
                    }

                    xpbd::solve_point_edge_collision_constraints(particles, pointEdgeCollisionConstraints, substep_time);

                    for (size_t i = 0; i < pointEdgeCollisionConstraints.point.size(); ++i)
                    {
                        auto &pecc = pointEdgeCollisionConstraints;
                        renderer::set_color(sf::Color::Red);
                        renderer::draw_circle(particles.pos[pecc.point[i]], 5);
                        renderer::set_color(sf::Color::Magenta);
                        renderer::draw_line(particles.pos[pecc.edge1[i]], particles.pos[pecc.edge2[i]]);
                    }
                }
                xpbd::update_velocities(particles, substep_time);
            }
        }

        renderer::set_color(sf::Color::Green);
        for (auto p : particles.pos)
            renderer::draw_circle(p, 3);

        for (size_t i = 0; i < distanceConstraints.i1.size(); ++i)
            renderer::draw_line(particles.pos[distanceConstraints.i1[i]], particles.pos[distanceConstraints.i2[i]]);

        renderer::set_color(sf::Color::Red);
        for (size_t i = 0; i < contourColliders.indices.size(); ++i)
        {
            std::vector<glm::vec2> points;

            for (auto id : contourColliders.indices[i])
                points.push_back(particles.pos[id]);

            renderer::draw_segmented_loop(points);
        }

        ImGui::SFML::Render(renderer::window);
        renderer::window.display();
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();
}
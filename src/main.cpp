#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>

#include <stdio.h>
#include "string.h"

#include "renderer.hpp"
#include "xpbd.hpp"

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
    xpbd::ContourCollider contour1;
    xpbd::ContourCollider contour2;

    add_particle(particles, {0, 0}, 1);
    add_particle(particles, {0, 100}, 1);
    add_particle(particles, {100, 100}, 1);
    add_particle(particles, {100, 0}, 1);

    add_distance_constraint(distanceConstraints, 0, 1, 0.002f, 100);
    add_distance_constraint(distanceConstraints, 1, 2, 0.002f, 100);
    add_distance_constraint(distanceConstraints, 2, 3, 0.002f, 100);
    add_distance_constraint(distanceConstraints, 3, 0, 0.002f, 100);

    add_distance_constraint(distanceConstraints, 0, 2, 0.0f, 140);

    contour1.particle_ids.push_back(0);
    contour1.particle_ids.push_back(1);
    contour1.particle_ids.push_back(2);
    contour1.particle_ids.push_back(3);

    add_particle(particles, {-100, -100}, 1);
    add_particle(particles, {-100, -200}, 1);
    add_particle(particles, {-200, -200}, 1);
    add_particle(particles, {-200, -100}, 1);

    add_distance_constraint(distanceConstraints, 4, 5, 0.002f, 100);
    add_distance_constraint(distanceConstraints, 5, 6, 0.002f, 100);
    add_distance_constraint(distanceConstraints, 6, 7, 0.002f, 100);
    add_distance_constraint(distanceConstraints, 7, 4, 0.002f, 100);

    add_distance_constraint(distanceConstraints, 4, 6, 0.005f, 140);
    add_distance_constraint(distanceConstraints, 5, 7, 0.005f, 140);

    add_distance_constraint(distanceConstraints, 0, 5, 0.1f, 0);

    contour2.particle_ids.push_back(4);
    contour2.particle_ids.push_back(5);
    contour2.particle_ids.push_back(6);
    contour2.particle_ids.push_back(7);

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
        ImGui::SliderFloat("timeScale", &timeScale, 0, 3, "%.3f");
        size_t min_value = 1;
        size_t max_value = 10;
        ImGui::SliderScalar("substeps", ImGuiDataType_U64, &substeps, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::SliderScalar("iterations", ImGuiDataType_U64, &iterations, &min_value, &max_value, "%zu", ImGuiSliderFlags_None);
        ImGui::End();
        ImGui::EndFrame();

        std::vector<xpbd::AABB> aabbs;
        std::vector<xpbd::AABBsIntersection> aabbs_intersections;
        std::vector<xpbd::PointEdgeCollisionConstraint> pointEdgeCollisionConstraints;

        float substep_time = deltaTick * timeScale / substeps;
        while (!paused && xpbd::should_tick(sec, deltaTick))
        {
            for (size_t s = 0; s < substeps; s++)
            {
                xpbd::iterate(particles, substep_time, gravity);
                xpbd::reset_constraints_lambdas(distanceConstraints.lambda);
                for (size_t i = 0; i < iterations; i++)
                {
                    xpbd::solve_distance_constraints(particles, distanceConstraints, substep_time);
                    std::vector<std::vector<size_t>> collider_particles_ids;
                    collider_particles_ids.push_back(contour1.particle_ids);
                    collider_particles_ids.push_back(contour2.particle_ids);
                    aabbs = xpbd::generate_colliders_aabbs(particles, collider_particles_ids);
                    aabbs_intersections = xpbd::find_aabbs_intersections(aabbs);

                    pointEdgeCollisionConstraints = xpbd::generate_contour_contour_collisions(particles, contour1, contour2);
                    for (auto pecc : pointEdgeCollisionConstraints)
                    {
                        xpbd::solve_point_edge_collision_constraint(particles, pecc, substep_time);
                        renderer::set_color(sf::Color::Red);
                        renderer::draw_circle(particles.pos[pecc.point], 5);
                        renderer::set_color(sf::Color::Magenta);
                        renderer::draw_line(particles.pos[pecc.edge1], particles.pos[pecc.edge2]);
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

        renderer::set_color(sf::Color::Blue);
        for (auto a : aabbs)
            renderer::draw_axis_aligned_bounding_box(a.l, a.r, a.b, a.t);

        // for some reason this printf turns off collisions
        // for (auto a : aabbs_intersections)
        //     printf("aabbs intersect index %d with %d\n", a.i1, a.i2);

        printf("\n");

        ImGui::SFML::Render(renderer::window);
        renderer::window.display();
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();
}
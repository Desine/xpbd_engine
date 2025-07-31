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
    glm::vec2 gravity = {0, -9.8f};
    xpbd::Particles particles;
    xpbd::DistanceConstraints distanceConstraints;

    add_particle(particles, {0, 100}, 1, 0);
    add_particle(particles, {100, 300}, 2, 0);
    add_particle(particles, {210, 500}, 3, 0);
    add_particle(particles, {330, 270}, 4, 0);

    add_distance_constraint(distanceConstraints, 0, 1, 0.005f, 200);
    add_distance_constraint(distanceConstraints, 1, 2, 0.002f, 150);
    add_distance_constraint(distanceConstraints, 2, 3, 0.003f, 200);
    add_distance_constraint(distanceConstraints, 2, 0, 0.001f, 200);


    add_particle(particles, {0, 100}, 1, 1);
    add_particle(particles, {100, 100}, 2, 1);
    add_particle(particles, {210, 100}, 3, 1);
    add_particle(particles, {230, 170}, 3, 1);

    add_distance_constraint(distanceConstraints, 4, 5, 0.005f, 200);
    add_distance_constraint(distanceConstraints, 5, 6, 0.002f, 150);
    add_distance_constraint(distanceConstraints, 6, 7, 0.003f, 200);
    add_distance_constraint(distanceConstraints, 6, 4, 0.001f, 200);

    renderer::setup_window();
    renderer::setup_view();
    renderer::setup_imgui();

    sf::Clock clock;
    while (renderer::window.isOpen())
    {
        sf::Time deltaTime = clock.restart();
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
        ImGui::End();
        ImGui::EndFrame();

        std::vector<xpbd::AABB> aabb;
        std::vector<xpbd::AABBsIntersections> aabbs_intersections;
        float substep_time = deltaTick / substeps;
        float iteration_time = substep_time / iterations;
        while (xpbd::should_tick(sec, deltaTick))
        {
            for (size_t s = 0; s < substeps; s++)
            {
                xpbd::iterate(particles, substep_time, gravity);
                xpbd::reset_distance_constraints_lambdas(distanceConstraints);
                for (size_t i = 0; i < iterations; i++)
                {
                    xpbd::solve_distance_constraints(particles, distanceConstraints, iteration_time);
                    aabb = xpbd::generate_aabbs(particles);
                    aabbs_intersections = xpbd::find_aabbs_intersections(aabb);
                }
                xpbd::update_velocities(particles, substep_time);
            }
        }

        sf::Color color = sf::Color::Green;

        for (size_t i = 0; i < particles.count; ++i)
        {
            renderer::draw_circle(particles.pos[i], 3);
        }
        for (size_t i = 0; i < distanceConstraints.count; ++i)
        {
            size_t p1 = distanceConstraints.i1[i];
            size_t p2 = distanceConstraints.i2[i];
            renderer::draw_line(particles.pos[p1], particles.pos[p2]);
        }
        for (auto a : aabb)
        {
            renderer::set_color(sf::Color::Blue);
            renderer::draw_axis_aligned_bounding_box(a.l, a.r, a.b, a.t);
        }
        for (auto a : aabbs_intersections)
        {
            printf("aabbs intersect index %d with %d\n", a.i1, a.i2);
        }
        printf("\n");
        

        ImGui::SFML::Render(renderer::window);
        renderer::window.display();
    }

    ImGui::SFML::Shutdown();
    renderer::window.close();
}
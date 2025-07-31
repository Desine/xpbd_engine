#include "glm/glm.hpp"
#include <vector>
#include <stdio.h>

namespace xpbd
{
#ifndef xpbd_define
#define xpbd_define
    struct Particles
    {
        std::vector<glm::vec2> pos;
        std::vector<glm::vec2> prevPos;
        std::vector<glm::vec2> vel;
        std::vector<float> w;
        std::vector<size_t> object_ids;
        size_t count = 0;
    };

    struct DistanceConstraints
    {
        std::vector<size_t> i1;
        std::vector<size_t> i2;
        std::vector<float> restDist;
        std::vector<float> compliance;
        std::vector<float> lambda;
        size_t count = 0;
    };

    struct AABB
    {
        float l, r, b, t;
        size_t object_id;
    };
    struct AABBsIntersections
    {
        size_t i1, i2;
    };
#endif // xpbd_define

    bool should_tick(float &sec, const float dt);
    void iterate(Particles &p, float dt, glm::vec2 gravity);
    void update_velocities(Particles &p, float dt);
    void add_particle(Particles &p, glm::vec2 pos, float mass, size_t object_id);
    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist);
    float get_distance_between_particles(Particles &particles, size_t i1, size_t i2);
    void solve_distance_constraints(Particles &p, DistanceConstraints &dc, float dt);
    void reset_distance_constraints_lambdas(DistanceConstraints &dc);
    std::vector<AABB> generate_aabbs(Particles &particles);
    std::vector<AABBsIntersections> find_aabbs_intersections(std::vector<AABB> aabb);
}
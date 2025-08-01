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
        std::vector<size_t> object_id;
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
    struct AABBsIntersection
    {
        size_t i1, i2;
    };
    struct ContourCollider
    {
        std::vector<size_t> particle_ids;
        float staticFriction;
        float kineticFriction;
        float compliance;
    };
    struct PointEdgeCollisionConstraint
    {
        size_t point;
        size_t edge1;
        size_t edge2;
        float staticFriction = 1.0f;
        float kineticFriction = 0.5f;
        float compliance = 0.0f;
    };
#endif // xpbd_define

    bool should_tick(float &sec, const float dt);
    void iterate(Particles &p, float dt, glm::vec2 gravity);
    void update_velocities(Particles &p, float dt);
    void add_particle(Particles &p, glm::vec2 pos, float mass, size_t object_id);
    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist);
    float get_distance_between_particles(const Particles &particles, size_t i1, size_t i2);
    void solve_distance_constraints( Particles &p,  DistanceConstraints &dc, float dt);
    void reset_distance_constraints_lambdas(DistanceConstraints &dc);
    std::vector<AABB> generate_aabbs(const Particles &particles);
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabb);
    std::vector<PointEdgeCollisionConstraint> generate_contour_contour_collisions(Particles &particles, ContourCollider &contourA, ContourCollider &contourB);
}
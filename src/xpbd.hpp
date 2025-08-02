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
        std::vector<glm::vec2> prev;
        std::vector<glm::vec2> vel;
        std::vector<float> w;
    };
    struct DistanceConstraints
    {
        std::vector<size_t> i1;
        std::vector<size_t> i2;
        std::vector<float> restDist;
        std::vector<float> compliance;
        std::vector<float> lambda;
    };
    struct VolumeConstraints
    {
        std::vector<std::vector<size_t>> indices;
        std::vector<float> restVolume;
        std::vector<float> compliance;
        std::vector<float> lambda;
    };
    struct AABB
    {
        float l, r, b, t;
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
    // struct ContourColliders
    // {
    //     std::vector<std::vector<size_t>> particle_ids;
    //     std::vector<float> staticFriction;
    //     std::vector<float> kineticFriction;
    //     std::vector<float> compliance;
    //     std::vector<size_t> collider_id;
    // };
    struct PointEdgeCollisionConstraint
    {
        size_t point;
        size_t edge1;
        size_t edge2;
        float staticFriction = 1.0f;
        float kineticFriction = 0.5f;
        float compliance = 0.0f;
        float lambda = 0.0f;
    };
#endif // xpbd_define

    bool should_tick(float &sec, const float dt);
    void iterate(Particles &p, float dt, glm::vec2 gravity);
    void update_velocities(Particles &p, float dt);
    void add_particle(Particles &p, glm::vec2 pos, float mass, glm::vec2 vel = {0, 0});
    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist);
    void add_distance_constraint_auto_restDist(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, Particles &p);
    void solve_distance_constraints(Particles &p, DistanceConstraints &dc, float dt);
    float compute_polygon_area(const std::vector<glm::vec2> &p);
    float compute_polygon_area(const Particles &p, std::vector<size_t> &indices);
    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance, float restPressure = 1);
    void solve_volume_constraints(Particles &p, VolumeConstraints &vc, float dt);
    void reset_constraints_lambdas(std::vector<float> &lambdas);
    std::vector<AABB> generate_colliders_aabbs(const Particles &p, const std::vector<std::vector<size_t>> particles_ids);
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabb);
    std::vector<PointEdgeCollisionConstraint> generate_contour_contour_collisions(Particles &p, ContourCollider &cA, ContourCollider &cB);
    void solve_point_edge_collision_constraint(Particles &p, PointEdgeCollisionConstraint &pecc, float dt);
}
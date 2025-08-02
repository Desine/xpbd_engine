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
    struct PolygonColliders
    {
        std::vector<std::vector<size_t>> indices;
        std::vector<float> staticFriction;
        std::vector<float> kineticFriction;
        std::vector<float> compliance;
    };
    struct PointColliders
    {
        std::vector<std::vector<size_t>> indices;
        std::vector<float> staticFriction;
        std::vector<float> kineticFriction;
        std::vector<float> compliance;
    };
    struct PointEdgeCollision
    {
        size_t point;
        size_t edge1;
        size_t edge2;
    };
    struct PointEdgeCollisionConstraints
    {
        std::vector<size_t> point;
        std::vector<size_t> edge1;
        std::vector<size_t> edge2;
        std::vector<float> staticFriction;
        std::vector<float> kineticFriction;
        std::vector<float> compliance;
        std::vector<float> lambda;
    };
#endif // xpbd_define

    bool should_tick(float &sec, const float dt);
    void iterate(Particles &p, float dt, glm::vec2 gravity);
    void update_velocities(Particles &p, float dt);
    void add_particle(Particles &p, glm::vec2 pos, float mass);
    void add_particle(Particles &p, glm::vec2 pos, float mass, glm::vec2 vel);
    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist);
    void add_distance_constraint_auto_restDist(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, Particles &p);
    void solve_distance_constraints(Particles &p, DistanceConstraints &dc, float dt);
    float compute_polygon_area(const std::vector<glm::vec2> &p);
    float compute_polygon_area(const Particles &p, std::vector<size_t> &indices);
    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance);
    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance, float restPressure);
    void solve_volume_constraints(Particles &p, VolumeConstraints &vc, float dt);
    void reset_constraints_lambdas(std::vector<float> &lambdas);
    std::vector<AABB> generate_particles_aabbs(const Particles &p, const std::vector<std::vector<size_t>> particles_ids);
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabbs);
    std::vector<AABBsIntersection> find_aabbs_intersections(const std::vector<AABB> &aabbs1, const std::vector<AABB> &aabbs2);
    void add_polygon_collider(xpbd::PolygonColliders &cc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance);
    void add_point_collider(xpbd::PointColliders &pc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance);
    void add_point_edge_collision_constraints_of_polygon_to_polygon_colliders(Particles &p, PointEdgeCollisionConstraints &pecc, const PolygonColliders &pc, const size_t cc_id1, const size_t cc_id2);
    void add_point_edge_collision_constraints_of_point_to_polygon_colliders(Particles &p, PointEdgeCollisionConstraints &pecc, const PointColliders &pointColliders, const PolygonColliders &polygonColliders, const size_t pointColliderId, const size_t polygonColliderId);
    void solve_point_edge_collision_constraints(Particles &p, PointEdgeCollisionConstraints &pecc, float dt);
}
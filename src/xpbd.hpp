#ifndef define_xpbd
#define define_xpbd
#include "glm/glm.hpp"
#include <vector>
#include <stdio.h>
#include <string>

#include "json_body_loader.hpp"

#include "Remotery.h"

namespace xpbd
{
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
        std::vector<float> compressionDamping;
        std::vector<float> extensionDamping;
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
    struct AABBsOverlap
    {
        size_t i1, i2;
        AABB box;
    };
    struct ColliderPoints
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
    struct PointPolygonCollision
    {
        std::vector<size_t> points;
        std::vector<size_t> polygon;
        float staticFriction;
        float kineticFriction;
        float compliance;
        AABB box;
    };
    struct PointEdgeCollisionConstraints
    {
        size_t point;
        size_t edge1;
        size_t edge2;
        float staticFriction;
        float kineticFriction;
        float compliance;
        float lambda;
    };

    class World
    {
    public:
        size_t substeps = 10;
        size_t iterations = 1;
        glm::vec2 gravity = {0, -9.8f};
        float timeScale = 10.f;
        float deltaTick = 1.f / 30.f;
        float sec = 0.0f;
        bool paused = false;
        bool stepOnce = false;

        Particles particles;
        DistanceConstraints distanceConstraints;
        VolumeConstraints volumeConstraints;
        ColliderPoints polygonColliders;
        ColliderPoints pointColliders;
        std::vector<PointPolygonCollision> collisions;

        void init();
        void spawnFromJson(const std::string &name, const glm::vec2 &position);
        void addPolygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
        void update(float realDelta);
    };

    bool should_tick(float &sec, const float &dt);
    void iterate(Particles &p, const float &dt, const glm::vec2 &gravity);
    void update_velocities(Particles &p, const float &dt);

    void add_particle(Particles &p, const glm::vec2 &pos, const float &mass);
    void add_particle(Particles &p, const glm::vec2 &pos, const float &mass, const glm::vec2 &vel);

    void reset_constraints_lambdas(std::vector<float> &lambdas);

    void add_distance_constraint(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, float restDist);
    void add_distance_constraint_auto_restDist(DistanceConstraints &dc, size_t i1, size_t i2, float compliance, Particles &p);
    void solve_distance_constraints(Particles &p, DistanceConstraints &dc, float dt);
    void apply_distance_constraints_damping(Particles &p, DistanceConstraints &dc, float dt);

    float compute_polygon_area(const std::vector<glm::vec2> &positions);
    float compute_polygon_area(const Particles &p, std::vector<size_t> &indices);
    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance);
    void add_volume_constraint(Particles &p, VolumeConstraints &vc, std::vector<size_t> indices, float compliance, float restPressure);
    void solve_volume_constraints(Particles &p, VolumeConstraints &vc, float dt);

    void add_polygon_collider(ColliderPoints &cc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance);
    void add_point_collider(ColliderPoints &pc, std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance);

    std::vector<AABB> generate_particles_aabbs(const Particles &p, const std::vector<std::vector<size_t>> particles_ids);
    std::vector<AABBsOverlap> create_aabbs_overlaps(const std::vector<AABB> &aabbs);
    std::vector<AABBsOverlap> create_aabbs_overlaps(const std::vector<AABB> &aabbs1, const std::vector<AABB> &aabbs2);

    std::vector<PointEdgeCollisionConstraints> get_point_edge_collision_constraints_of_point_to_polygon_colliders(const Particles &p, const PointPolygonCollision &collision);
    std::vector<PointEdgeCollisionConstraints> get_point_edge_collision_constraints_of_point_to_polygon_colliders_parallel(const Particles &particles, const std::vector<PointPolygonCollision> &collisions);
    void solve_point_edge_collision_constraints(Particles &p, std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
    void apply_point_edge_collision_constraints_kinetic_friction(Particles &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
}
#endif // define_xpbd
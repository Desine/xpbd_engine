#ifndef define_xpbd
#define define_xpbd
#include "glm/glm.hpp"
#include <vector>
#include <thread>
#include <condition_variable>
#include <chrono>
#include <queue>
#include <functional>

#include "json_body_loader.hpp"
#include "spatial_hashing.hpp"

namespace xpbd
{
    struct Particles
    {
        std::vector<glm::vec2> pos;
        std::vector<glm::vec2> prev;
        std::vector<glm::vec2> vel;
        std::vector<float> w;
    };
    struct DistanceConstraint
    {
        size_t i1;
        size_t i2;
        float restDist;
        float compressionDamping;
        float extensionDamping;
        float compliance;
        float lambda;
    };
    struct VolumeConstraint
    {
        std::vector<size_t> indices;
        float restVolume;
        float compliance;
        float lambda;
    };
    struct ColliderPoints
    {
        std::vector<size_t> indices;
        float staticFriction;
        float kineticFriction;
        float compliance;
    };
    struct PointEdgeCollision
    {
        size_t point;
        size_t edge1;
        size_t edge2;
    };
    struct PolygonCache
    {
        std::vector<glm::vec2> positions;
        size_t positions_size;
        std::vector<glm::vec2> edges;
        std::vector<float> edgesLen2;
    };
    struct PointsPolygonCollision
    {
        std::vector<size_t> points;
        std::vector<size_t> polygon;
        float staticFriction;
        float kineticFriction;
        float compliance;
        AABB box;
        PolygonCache cache;
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

    struct WorldDraw
    {
        Particles particles;
        std::vector<DistanceConstraint> distanceConstraints;
        std::vector<VolumeConstraint> volumeConstraints;
        std::vector<ColliderPoints> polygonColliders;
        std::vector<ColliderPoints> pointsColliders;
    };

    class World
    {
    public:
        size_t substeps = 10;
        size_t iterations = 1;
        glm::vec2 gravity = {0, -9.8f};
        float timeScale = 10.f;
        float tickRate = 30.f;
        std::chrono::duration<double> deltaTick;

        Particles particles;
        std::vector<DistanceConstraint> distanceConstraints;
        std::vector<VolumeConstraint> volumeConstraints;
        std::vector<ColliderPoints> polygonColliders;
        std::vector<ColliderPoints> pointsColliders;

        // std::vector<AABB> polygonColliders_aabbs;
        // std::vector<AABB> pointsColliders_aabbs;
        SpatialHashAABB spatialHashAABB;
        // std::vector<PolygonCache> polygonsHash;
        // std::vector<ColliderPoints> combinedColliderPoints;
        // std::vector<AABB> combinedAABB;
        // std::vector<PointEdgeCollisionConstraints> pecc;

        std::mutex drawMutex;
        WorldDraw draw;
        bool newFrame = false;

        enum class UpdateState
        {
            Running,
            Paused,
            StepOnce,
            Stopped
        };
        std::thread updateThread;
        std::mutex updateMutex;
        std::condition_variable updateCv;
        UpdateState updateState = UpdateState::Paused;
        std::chrono::time_point<std::chrono::_V2::steady_clock, std::chrono::duration<double, std::chrono::_V2::system_clock::period>>
            nextUpdateTimePoint;

        std::queue<std::function<void()>> queue;
        std::mutex functionQueueMutex;

        World();
        ~World();

        void submit_draw();

        void step();
        void update_loop();
        void pause_update();
        void resume_update();
        void step_once_update();

        void enqueue_function(std::function<void()> f);
        void execute_enqueued_functions();

        void init();
        void enqueue_init();

        void set_tickRate(size_t tickRate);
        void enqueue_set_tickRate(size_t tickRate);

        void add_particle(const glm::vec2 &pos, const float mass, const glm::vec2 &vel);
        void enqueue_add_particle(const glm::vec2 &pos, const float mass, const glm::vec2 &vel);

        void add_distance_constraint(size_t i1, size_t i2, float compliance, float restDist);
        void enqueue_add_distance_constraint(size_t i1, size_t i2, float compliance, float restDist);
        void add_volume_constraint(const std::vector<size_t> &indices, float compliance, float restPressure);
        void enqueue_add_volume_constraint(std::vector<size_t> indices, float compliance, float restPressure);

        void add_polygon_collider(const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance);
        void enqueue_add_polygon_collider(std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance);
        void add_points_collider(const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance);
        void enqueue_add_points_collider( std::vector<size_t> indices, float staticFriction, float kineticFriction, float compliance);

        void spawn_from_json(const std::string &name, const glm::vec2 &position);
        void enqueue_spawn_from_json( std::string name,  glm::vec2 position);
        void spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
        void enqueue_spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
    };

    void iterate(Particles &p, float dt, const glm::vec2 &gravity);
    void update_velocities(Particles &p, float dt);

    void reset_distance_constraints_lambdas(std::vector<DistanceConstraint> &dc);
    void reset_volume_constraints_lambdas(std::vector<VolumeConstraint> &vc);

    void solve_distance_constraint(Particles &p, DistanceConstraint &dc, float dt);
    void solve_distance_constraints(Particles &p, std::vector<DistanceConstraint> &dc, float dt);
    void apply_distance_constraints_damping(Particles &p, std::vector<DistanceConstraint> &dc, float dt);

    float compute_polygon_area(const std::vector<glm::vec2> &positions);
    float compute_polygon_area(const Particles &p, const std::vector<size_t> &indices);
    void solve_volume_constraint(Particles &p, VolumeConstraint &vc, float dt);
    void solve_volume_constraints(Particles &p, std::vector<VolumeConstraint> &vc, float dt);

    void generate_collider_points_aabbs(std::vector<AABB> &out, const Particles &p, const std::vector<ColliderPoints> &colliderPoints);

    void add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(std::vector<PointEdgeCollisionConstraints> &pecc, const Particles &p, const PointsPolygonCollision &collision);
    void solve_point_edge_collision_constraints(Particles &p, std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
    void apply_point_edge_collision_constraints_kinetic_friction(Particles &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
}
#endif // define_xpbd
#ifndef XPBD
#define XPBD
#include "glm/glm.hpp"
#include <vector>
#include <condition_variable>
#include <chrono>
#include <queue>
#include <functional>

#include "spatial_hashing.hpp"

namespace Commands
{
    struct init
    {
    };

    struct set_tickRate
    {
        float tickRate;
    };

    struct add_particle
    {
        glm::vec2 pos;
        float mass;
        glm::vec2 vel;
    };

    struct add_distance_constraint
    {
        size_t i1;
        size_t i2;
        float compliance;
        float restDist;
    };

    struct add_volume_constraint
    {
        std::vector<size_t> indices;
        float compliance;
        float restPressure;
    };

    struct add_polygon_collider
    {
        std::vector<size_t> indices;
        float staticFriction;
        float kineticFriction;
        float compliance;
        int8_t layer;
        int8_t mask;
    };

    struct add_points_collider
    {
        std::vector<size_t> indices;
        float staticFriction;
        float kineticFriction;
        float compliance;
        int8_t layer;
        int8_t mask;
    };

    struct spawn_from_json
    {
        std::string name;
        glm::vec2 position;
    };

    struct spawn_polygon
    {
        glm::vec2 pos;
        float radius;
        size_t segments;
        float mass;
        float compliance;
    };

    using Command = std::variant<
        init,
        set_tickRate,
        add_particle,
        add_distance_constraint,
        add_volume_constraint,
        add_polygon_collider,
        add_points_collider,
        spawn_from_json,
        spawn_polygon>;
}

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
    struct AngleConstraint{
        size_t i1;
        size_t i2;
        size_t i3;
        size_t i4;
        float compliance;
        float restAngle;
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
        int8_t layer;
        int8_t mask;
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

    struct Cluster
    {
        Particles particles;
        std::vector<DistanceConstraint> distanceConstraints;
        std::vector<AngleConstraint> angleConstraints;
        std::vector<VolumeConstraint> volumeConstraints;
        std::vector<ColliderPoints> polygonColliders;
        std::vector<ColliderPoints> pointsColliders;
    };

    class World
    {
    public:
        size_t substeps = 20;
        size_t iterations = 1;
        glm::vec2 gravity = {0, -9.8f};
        float timeScale = 1.f;
        float tickRate = 30.f;
        std::chrono::duration<double> deltaTick;
        bool reversePass = false;

        Cluster cluster;
        std::mutex clusterMutex;
        Cluster clusterAdd;
        std::mutex clusterAddMutex;

        // std::vector<AABB> polygonColliders_aabbs;
        // std::vector<AABB> pointsColliders_aabbs;
        SpatialHashAABB spatialHashAABB;
        // std::vector<PolygonCache> polygonsHash;
        // std::vector<ColliderPoints> combinedColliderPoints;
        // std::vector<AABB> combinedAABB;
        // std::vector<PointEdgeCollisionConstraints> pecc;

        Cluster clusterDraw;
        std::mutex clusterDrawMutex;
        bool newFrame = false;
        std::chrono::duration<float> lastStepTime;

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
        UpdateState updateState = UpdateState::Running;
        std::chrono::time_point<std::chrono::_V2::steady_clock, std::chrono::duration<double, std::chrono::_V2::system_clock::period>>
            nextUpdateTimePoint;

        std::queue<Commands::Command> commandQueue;
        std::mutex commandQueueMutex;
        void enqueue_command(Commands::Command variant);
        void execute_enqueued_commands();

        World();
        ~World();

        void step();
        void update_loop();
        void pause_update();
        void resume_update();
        void step_once_update();

        void init();

    private:
        void set_tickRate(size_t tickRate);

        void spawn_from_json(const std::string &name, const glm::vec2 &position);
        void spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
    };

    void merge_clusters(Cluster &into, Cluster &from);

    void add_particle(Cluster &cluster, const glm::vec2 &pos, const float mass, const glm::vec2 &vel);

    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance);
    void add_distance_constraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float restDist);
    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2,size_t i3, size_t i4, float compliance);
    void add_angle_constraint(Cluster &cluster, size_t i1, size_t i2,size_t i3, size_t i4, float compliance, float restAngle);
    void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance);
    void add_volume_constraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float restPressure);

    void add_polygon_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);
    void add_points_collider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);

    void spawn_from_json(Cluster &cluster, const std::string &name, const glm::vec2 &position);
    void spawn_polygon(Cluster &cluster, glm::vec2 pos, float radius, size_t segments, float mass, float compliance);

    void iterate(Particles &p, float dt, const glm::vec2 &gravity);
    void update_velocities(Particles &p, float dt);

    void reset_distance_constraints_lambdas(std::vector<DistanceConstraint> &dc);
    void reset_volume_constraints_lambdas(std::vector<VolumeConstraint> &vc);

    void solve_distance_constraint(Particles &p, DistanceConstraint &dc, float dt);
    void solve_distance_constraints(Particles &p, std::vector<DistanceConstraint> &dc, float dt, bool reversePass);
    void apply_distance_constraints_damping(Particles &p, std::vector<DistanceConstraint> &dc, float dt);

    void solve_angle_constraint(Particles &p, AngleConstraint &ac, float dt);
    void solve_angle_constraints(Particles &p, std::vector<AngleConstraint> &ac, float dt);

    float compute_polygon_area(const std::vector<glm::vec2> &positions);
    float compute_polygon_area(const Particles &p, const std::vector<size_t> &indices);
    void solve_volume_constraint(Particles &p, VolumeConstraint &vc, float dt);
    void solve_volume_constraints(Particles &p, std::vector<VolumeConstraint> &vc, float dt);

    void generate_collider_points_aabbs(std::vector<AABB> &out, const Particles &p, const std::vector<ColliderPoints> &colliderPoints);

    void add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(std::vector<PointEdgeCollisionConstraints> &pecc, const Particles &p, const PointsPolygonCollision &collision);
    void solve_point_edge_collision_constraints(Particles &p, std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
    void apply_point_edge_collision_constraints_kinetic_friction(Particles &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
}
#endif // XPBD
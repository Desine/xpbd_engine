#ifndef XPBD
#define XPBD
#include "glm/glm.hpp"
#include <vector>
#include <condition_variable>
#include <chrono>
#include <queue>
#include <functional>
#include <variant>
#include <thread>
#include <queue>

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
        float minDist;
        float maxDist;
        float compressionDamping;
        float extensionDamping;
        float compliance;
        float lambda = 0.0f;
    };
    struct AngleConstraint
    {
        size_t i1, i2;
        size_t i3, i4;
        float minAngle;
        float maxAngle;
        float compliance;
        float lambda = 0.0f;
    };
    struct VolumeConstraint
    {
        std::vector<size_t> indices;
        float restVolume;
        float compliance;
        float lambda = 0.0f;
    };

    struct PointsCollider
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
        float lambdaN = 0.0f;
        float lambdaT = 0.0f;
    };

    struct AngularData
    {
        glm::vec2 grad1;
        glm::vec2 grad2;
        float invI1;
        float invI2;
        float invI;
        float angle;
        float omega;
    };
    struct AngularMotor
    {
        size_t i1, i2;
        size_t i3, i4;
        float targetRPM;
        float maxTorque;
    };
    struct AngularServo
    {
        size_t i1, i2;
        size_t i3, i4;
        float targetAngle;
        float kp;
        float kd;
        float maxTorque;
    };

    struct Cluster
    {
        std::vector<std::pair<std::string, size_t>> particleIds;
        Particles particles;
        std::vector<DistanceConstraint> distanceConstraints;
        std::vector<AngleConstraint> angleConstraints;
        std::vector<VolumeConstraint> volumeConstraints;
        std::vector<PointsCollider> polygonColliders;
        std::vector<PointsCollider> pointsColliders;
        std::vector<AngularMotor> angularMotors;
        std::vector<AngularServo> angularServos;
    };

    struct init_command
    {
    };

    struct set_tickRate_command
    {
        float tickRate;
    };

    struct add_particle_command
    {
        glm::vec2 pos;
        float mass;
        glm::vec2 vel;
    };

    struct add_distanceConstraint_command
    {
        DistanceConstraint constraint;
    };
    struct add_angleConstraint_command
    {
        AngleConstraint constraint;
    };
    struct add_volumeConstraint_command
    {
        VolumeConstraint constraint;
    };

    struct add_polygonCollider_command
    {
        PointsCollider collider;
    };

    struct add_pointsCollider_command
    {
        PointsCollider collider;
    };

    struct spawn_from_json_command
    {
        std::string name;
        glm::vec2 position;
    };

    struct spawn_polygon_command
    {
        glm::vec2 pos;
        float radius;
        size_t segments;
        float mass;
        float compliance;
    };
    struct add_angularMotor_command
    {
        AngularMotor motor;
    };
    struct update_angularMotor_command
    {
        size_t index;
        AngularMotor motor;
    };
    struct add_angularServo_command
    {
        AngularServo servo;
    };
    struct update_angularServo_command
    {
        size_t index;
        AngularServo servo;
    };

    using Command = std::variant<
        init_command,
        set_tickRate_command,
        add_particle_command,
        add_distanceConstraint_command,
        add_volumeConstraint_command,
        add_polygonCollider_command,
        add_pointsCollider_command,
        spawn_from_json_command,
        spawn_polygon_command,
        add_angularMotor_command,
        update_angularMotor_command,
        add_angularServo_command,
        update_angularServo_command>;

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

        std::queue<Command> commandQueue;
        std::mutex commandQueueMutex;
        void enqueue_command(Command variant);
        void execute_enqueued_commands();

        World();
        ~World();

        void step();
        void update_loop();
        void pause_update();
        void resume_update();
        void step_once_update();
        void stop_update();

        void init();

    private:
        void set_tickRate(size_t tickRate);

        void spawn_from_json(const std::string &name, const glm::vec2 &position);
        void spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
    };

    void merge_clusters(Cluster &into, Cluster &from);

    void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass, glm::vec2 vel);

    void add_distanceConstraint(Cluster &cluster, size_t i1, size_t i2, float compliance);
    void add_distanceConstraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float minDist, float maxDist);
    void add_angleConstraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance);
    void add_angleConstraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float minAngle, float maxAngle);
    void add_volumeConstraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance);
    void add_volumeConstraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float restPressure);

    void add_polygonCollider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);
    void add_pointsCollider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask);

    void add_angularMotor(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float targetVelocity, float maxTorque);
    void add_angularServo(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float targetAngle, float kp, float kd, float maxTorque);

    void spawn_from_json(Cluster &cluster, const std::string &name, const glm::vec2 &position);
    void spawn_polygon(Cluster &cluster, glm::vec2 pos, float radius, size_t segments, float mass, float compliance);

    void iterate(Particles &p, float dt, const glm::vec2 &gravity);
    void update_velocities(Particles &p, float dt);

    void reset_distanceConstraints_lambdas(std::vector<DistanceConstraint> &dc);
    void reset_volumeConstraints_lambdas(std::vector<VolumeConstraint> &vc);

    void solve_distanceConstraint(Particles &p, DistanceConstraint &dc, float dt);
    void solve_distanceConstraints(Particles &p, std::vector<DistanceConstraint> &dc, float dt, bool reversePass);
    void apply_distanceConstraints_damping(Particles &p, std::vector<DistanceConstraint> &dc, float dt);

    void solve_angleConstraint(Particles &p, AngleConstraint &ac, float dt);
    void solve_angleConstraints(Particles &p, std::vector<AngleConstraint> &ac, float dt);

    float compute_polygon_area(const std::vector<glm::vec2> &positions);
    float compute_polygon_area(const Particles &p, const std::vector<size_t> &indices);
    void solve_volumeConstraint(Particles &p, VolumeConstraint &vc, float dt);
    void solve_volumeConstraints(Particles &p, std::vector<VolumeConstraint> &vc, float dt);

    bool get_angularData(const Particles &p, size_t i1, size_t i2, size_t i3, size_t i4, AngularData &data);
    void apply_angular_velocity_correction(Particles &p, size_t i1, size_t i2, size_t i3, size_t i4, float deltaL);
    void apply_torque(Particles &p, size_t i1, size_t i2, size_t i3, size_t i4, float torque, float dt);
    void solve_angularMotor(Particles &p, const AngularMotor &motor, float dt);
    void solve_angularServo(Particles &p, const AngularServo &servo, float dt);

    void generate_colliderPoints_aabbs(std::vector<AABB> &out, const Particles &p, const std::vector<PointsCollider> &colliderPoints);

    void add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(std::vector<PointEdgeCollisionConstraints> &pecc, const Particles &p, const PointsPolygonCollision &collision);
    void solve_pointEdgeCollisionConstraints(Particles &p, std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
    void apply_pointEdgeCollisionConstraints_kinetic_friction(Particles &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
}
#endif // XPBD
#pragma once
#include "types.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"
#include "cluster.hpp"
#include "commands.hpp"
#include "../solver/spatial_hashing.hpp"

#include "glm/glm.hpp"
#include <vector>
#include <condition_variable>
#include <chrono>
#include <queue>
#include <thread>

namespace xpbd
{
    class World
    {
    public:
        size_t substeps = 10;
        size_t iterations = 1;
        glm::vec2 gravity = {0, -9.8f};
        float timeScale = 1.f;
        float tickRate = 60.f;
        std::chrono::duration<double> deltaTick;
        bool reversePass = false;

        Cluster cluster;
        std::mutex clusterMutex;
        Cluster clusterAdd;
        std::mutex clusterAddMutex;

        std::vector<AABB> polygonColliders_aabbs;
        std::vector<AABB> pointsColliders_aabbs;
        SpatialHashAABB spatialHashAABB;
        std::vector<PointsCollider> combinedPointsColliders;
        std::vector<AABB> combinedAABB;
        std::vector<PointEdgeCollisionConstraints> pecc;
        std::vector<size_t> overlapping_aabb_ids;

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
} // namespace xpbd
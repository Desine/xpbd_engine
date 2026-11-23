#pragma once
#include "scene.hpp"
#include "commands.hpp"
#include "cluster.hpp"
#include "../solver/spatial_hashing.hpp"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace xpbd
{
class World
{
public:
    size_t substeps = 10;
    size_t iterations = 1;
    glm::vec2 gravity = {0.0f, -9.8f};
    float timeScale = 1.0f;
    float tickRate = 60.0f;
    std::chrono::duration<double> deltaTick;

    Scene scene;
    std::mutex sceneMutex;

    std::vector<AABB> polygonColliders_aabbs;
    std::vector<AABB> pointsColliders_aabbs;
    SpatialHashAABB spatialHashAABB;
    std::vector<PointEdgeCollisionConstraints> pecc;
    std::vector<size_t> overlapping_aabb_ids;

    cluster::Cluster clusterDraw;
    std::mutex clusterDrawMutex;
    std::atomic_bool newFrame{false};
    std::chrono::duration<double> lastStepTime{};

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
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> nextUpdateTimePoint{};

    std::queue<Command> commandQueue;
    std::mutex commandQueueMutex;

    World();
    ~World();
    World(const World &) = delete;
    World &operator=(const World &) = delete;

    void enqueue_command(Command command);
    void execute_enqueued_commands();
    void step();
    void update_loop();
    void pause_update();
    void resume_update();
    void step_once_update();
    void stop_update();
    void init();

private:
    void set_tickRate(float rate);
    void spawn_from_json(const std::string &name, const glm::vec2 &position);
    void spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance);
    cluster::Cluster make_draw_snapshot() const;
};
} // namespace xpbd

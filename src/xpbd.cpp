#include "xpbd.hpp"
#include "cstdio"
#include <algorithm>
#include <omp.h>
#include <thread>
#include <thread>
#include <print>

#include "utils.hpp"
#include "renderer.hpp"

#include "profiler.hpp"
#include "json_body_loader.hpp"

namespace xpbd
{
    void iterate(Particles &p, float dt, const glm::vec2 &gravity)
    {
        for (size_t i = 0; i < p.pos.size(); ++i)
        {
            if (p.w[i] == 0.0f)
                continue;

            p.vel[i] += gravity * dt;
            p.prev[i] = p.pos[i];
            p.pos[i] += p.vel[i] * dt;
        }
    }

    void update_velocities(Particles &p, float dt)
    {
        for (size_t i = 0; i < p.pos.size(); ++i)
            p.vel[i] = (p.pos[i] - p.prev[i]) / dt;
    }

    void reset_distanceConstraints_lambdas(std::vector<DistanceConstraint> &dc)
    {
        for (auto &c : dc)
            c.lambda = 0;
    }
    void reset_angleConstraints_lambdas(std::vector<AngleConstraint> &ac)
    {
        for (auto &c : ac)
            c.lambda = 0;
    }
    void reset_volumeConstraints_lambdas(std::vector<VolumeConstraint> &vc)
    {
        for (auto &c : vc)
            c.lambda = 0;
    }

    void solve_distanceConstraints(Particles &p, std::vector<DistanceConstraint> &dc, float dt, bool reversePass)
    {
        if (reversePass)
        {
            for (auto it = dc.rbegin(); it != dc.rend(); it++)
                solve_distanceConstraint(p, *it, dt);
        }
        else
        {
            for (auto &c : dc)
                solve_distanceConstraint(p, c, dt);
        }
    }
    void solve_distanceConstraint(Particles &p, DistanceConstraint &dc, float dt)
    {
        size_t i1 = dc.i1;
        size_t i2 = dc.i2;

        glm::vec2 &p1 = p.pos[i1];
        glm::vec2 &p2 = p.pos[i2];

        float w1 = p.w[i1];
        float w2 = p.w[i2];

        glm::vec2 delta = p1 - p2;
        float len = glm::length(delta);
        if (len < 1e-6f)
            return;

        float C = len - dc.restDist;
        glm::vec2 grad = delta / len;

        float alphaTilde = dc.compliance / (dt * dt);
        float denom = w1 + w2 + alphaTilde;
        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * dc.lambda) / denom;
        dc.lambda += deltaLambda;

        glm::vec2 correction = deltaLambda * grad;
        p1 += w1 * correction;
        p2 -= w2 * correction;
    }
    void apply_distanceConstraints_damping(Particles &p, std::vector<DistanceConstraint> &dc, float dt)
    {
        // does not work
        for (size_t i = 0; i < dc.size(); ++i)
        {
            size_t i1 = dc[i].i1;
            size_t i2 = dc[i].i2;
            if (i1 >= p.pos.size() || i2 >= p.pos.size())
                continue;

            glm::vec2 delta = p.pos[i1] - p.pos[i2];
            float len = glm::length(delta);
            if (len < 1e-6f)
                continue;

            glm::vec2 dir = delta / len;
            float rest = dc[i].restDist;
            float C = len - rest;

            float w1 = p.w[i];
            float w2 = p.w[i];
            float wSum = w1 + w2;
            if (wSum < 1e-6f)
                continue;

            glm::vec2 relVel = p.vel[i1] - p.vel[i2];
            float relVelAlongConstraint = glm::dot(relVel, dir);

            float damping = (C < 0.0f) ? dc[i].compressionDamping : dc[i].extensionDamping;

            float impulseMag = -damping * relVelAlongConstraint;
            glm::vec2 impulse = impulseMag * dir;

            p.vel[i1] += w1 / wSum * impulse;
            p.vel[i2] -= w2 / wSum * impulse;
        }
    }

    void solve_angleConstraints(Particles &p, std::vector<AngleConstraint> &ac, float dt)
    {
        for (auto &c : ac)
            solve_angleConstraint(p, c, dt);
    }

    void solve_angleConstraint(Particles &p, AngleConstraint &ac, float dt)
    {
        size_t i1 = ac.i1, i2 = ac.i2;
        size_t i3 = ac.i3, i4 = ac.i4;

        glm::vec2 dir1 = p.pos[i2] - p.pos[i1];
        glm::vec2 dir2 = p.pos[i4] - p.pos[i3];

        float len1Sq = glm::dot(dir1, dir1);
        float len2Sq = glm::dot(dir2, dir2);
        if (len1Sq < 1e-6f || len2Sq < 1e-6f)
            return;

        constexpr float TWO_PI = M_PI * 2.0f;
        float angle = Utils::relative_angle_rad(dir1, dir2);
        if (angle < 0.0f)
            angle += TWO_PI;

        float minAngle = std::fmod(ac.minAngle, TWO_PI);
        float maxAngle = std::fmod(ac.maxAngle, TWO_PI);

        if (minAngle < 0.0f)
            minAngle += TWO_PI;
        if (maxAngle < 0.0f)
            maxAngle += TWO_PI;

        bool inside;
        if (minAngle <= maxAngle)
            inside = angle >= minAngle && angle <= maxAngle;
        else
            inside = angle >= minAngle || angle <= maxAngle;
        if (inside)
        {
            ac.lambda = 0.0f;
            return;
        }

        auto signedAngleDiff = [&](float a, float b)
        {
            return std::atan2(std::sin(a - b), std::cos(a - b));
        };

        float diffMin = signedAngleDiff(angle, minAngle);
        float diffMax = signedAngleDiff(angle, maxAngle);

        float C;
        if (std::abs(diffMin) < std::abs(diffMax))
            C = diffMin;
        else
            C = diffMax;

        glm::vec2 grad1 = Utils::perpCCW2d(dir1) / len1Sq;
        glm::vec2 grad2 = Utils::perpCCW2d(dir2) / len2Sq;

        float w1 = p.w[i1], w2 = p.w[i2];
        float w3 = p.w[i3], w4 = p.w[i4];

        float alphaTilde = ac.compliance / (dt * dt);
        float denom =
            (w1 + w2) * glm::dot(grad1, grad1) +
            (w3 + w4) * glm::dot(grad2, grad2) +
            alphaTilde;
        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * ac.lambda) / denom;
        ac.lambda += deltaLambda;

        p.pos[i1] += w1 * grad1 * deltaLambda;
        p.pos[i2] -= w2 * grad1 * deltaLambda;
        p.pos[i3] -= w3 * grad2 * deltaLambda;
        p.pos[i4] += w4 * grad2 * deltaLambda;
    }

    float compute_polygon_area(const std::vector<glm::vec2> &positions)
    {
        float area = 0.0f;
        size_t n = positions.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = positions[j];
            const glm::vec2 &p1 = positions[i];
            area += Utils::cross(p0, p1);
        }
        return 0.5f * area;
    }
    float compute_polygon_area(const Particles &p, const std::vector<size_t> &indices)
    {
        float area = 0.0f;
        size_t n = indices.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = p.pos[indices[j]];
            const glm::vec2 &p1 = p.pos[indices[i]];
            area += Utils::cross(p0, p1);
        }
        return 0.5f * area;
    }

    void solve_volumeConstraints(Particles &p, std::vector<VolumeConstraint> &vc, float dt)
    {
        for (auto &c : vc)
            solve_volumeConstraint(p, c, dt);
    }
    void solve_volumeConstraint(Particles &p, VolumeConstraint &vc, float dt)
    {
        std::vector<size_t> &indices = vc.indices;
        size_t N = indices.size();

        float volume = compute_polygon_area(p, indices);
        float C = volume - vc.restVolume;

        float denom = 0.0f;
        std::vector<glm::vec2> grads(N);

        for (size_t i = 0; i < N; ++i)
        {
            const glm::vec2 &pi_prev = p.pos[indices[(i + N - 1) % N]];
            const glm::vec2 &pi_next = p.pos[indices[(i + 1) % N]];

            glm::vec2 grad = 0.5f * Utils::perpCW2d(pi_next - pi_prev);
            grads[i] = grad;

            float wi = p.w[indices[i]];
            denom += wi * glm::dot(grad, grad);
        }

        float alphaTilde = vc.compliance / (dt * dt);
        denom += alphaTilde;
        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * vc.lambda) / denom;
        vc.lambda += deltaLambda;

        for (size_t i = 0; i < N; ++i)
        {
            uint32_t idx = indices[i];
            float wi = p.w[idx];

            p.pos[idx] += wi * deltaLambda * grads[i];
        }
    }

    void generate_colliderPoints_aabbs(
        std::vector<AABB> &out,
        const Particles &p,
        const std::vector<PointsCollider> &colliderPoints)
    {
        out.clear();
        out.reserve(colliderPoints.size());

        for (auto &c : colliderPoints)
        {
            AABB current;
            current.l = std::numeric_limits<float>::max();
            current.r = -std::numeric_limits<float>::max();
            current.b = std::numeric_limits<float>::max();
            current.t = -std::numeric_limits<float>::max();
            current.layer = c.layer;
            current.mask = c.mask;

            for (auto i : c.indices)
            {
                glm::vec2 pos = p.pos[i];
                if (pos.x < current.l)
                    current.l = pos.x;
                if (pos.x > current.r)
                    current.r = pos.x;
                if (pos.y < current.b)
                    current.b = pos.y;
                if (pos.y > current.t)
                    current.t = pos.y;
            }

            out.emplace_back(current);
        }
    }

    bool point_in_polygon(const glm::vec2 &point, const std::vector<glm::vec2> &positions)
    {
        int windingNumber = 0;
        const size_t n = positions.size();

        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &v1 = positions[i];
            const glm::vec2 &v2 = positions[j];
            float cross = Utils::cross(v2 - v1, point - v1);

            if (v1.y <= point.y)
            {
                if (v2.y > point.y && cross > 0)
                    ++windingNumber;
            }
            else
            {
                if (v2.y <= point.y && cross < 0)
                    --windingNumber;
            }
        }
        return windingNumber != 0;
    }
    void get_PointEdgeCollisions_of_points_inside_polygon(
        std::vector<PointEdgeCollision> &out,
        const Particles &p,
        const std::vector<size_t> &pointIndices,
        const std::vector<size_t> &polygonIndices,
        const PolygonCache &cache)
    {
        out.reserve(pointIndices.size());

        const auto &polygonPositions = cache.positions;
        const auto &edges = cache.edges;
        const auto &edgesLen2 = cache.edgesLen2;
        size_t n = cache.positions_size;

        for (size_t pid : pointIndices)
        {
            const glm::vec2 &point = p.pos[pid];

            if (!point_in_polygon(point, polygonPositions))
                continue;

            float minDist = std::numeric_limits<float>::max();
            size_t nearestEdge = 0;

            for (size_t i = 0; i < n; ++i)
            {
                if (cache.edgesLen2[i] < 1e-12f)
                    continue;

                glm::vec2 toPoint = point - polygonPositions[i];
                float t = glm::dot(toPoint, edges[i]) / edgesLen2[i];
                t = glm::clamp(t, 0.0f, 1.0f);
                glm::vec2 closest = polygonPositions[i] + edges[i] * t;

                float dist = glm::dot(point - closest, point - closest);
                if (dist < minDist)
                {
                    minDist = dist;
                    nearestEdge = i;
                }
            }

            out.emplace_back(pid,
                             polygonIndices[nearestEdge],
                             polygonIndices[(nearestEdge + 1) % n]);
        }
    }
    void add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(
        std::vector<PointEdgeCollisionConstraints> &pecc,
        const Particles &particles,
        const PointsPolygonCollision &collision)
    {
        std::vector<size_t> filteredPointIndices;
        filteredPointIndices.reserve(collision.points.size());
        for (auto i : collision.points)
            if (collision.box.contains_point(particles.pos[i]))
                filteredPointIndices.emplace_back(i);

        if (filteredPointIndices.empty())
            return;

        std::vector<PointEdgeCollision> detections;
        detections.reserve(100);
        get_PointEdgeCollisions_of_points_inside_polygon(detections, particles, filteredPointIndices, collision.polygon, collision.cache);

        pecc.reserve(pecc.size() + detections.size());
        for (const auto &c : detections)
            pecc.emplace_back(c.point, c.edge1, c.edge2, collision.staticFriction, collision.kineticFriction, collision.compliance, 0.0f, 0.0f);
    }

    void solve_pointEdgeCollisionConstraints(
        Particles &p,
        std::vector<PointEdgeCollisionConstraints> &pecc,
        float dt)
    {
        for (auto &c : pecc)
        {
            size_t ip = c.point;
            size_t ie0 = c.edge1;
            size_t ie1 = c.edge2;

            glm::vec2 &P = p.pos[ip];
            glm::vec2 &Pprev = p.prev[ip];
            glm::vec2 &E0 = p.pos[ie0];
            glm::vec2 &E1 = p.pos[ie1];
            glm::vec2 &E0prev = p.prev[ie0];
            glm::vec2 &E1prev = p.prev[ie1];

            float wp = p.w[ip];
            float we0 = p.w[ie0];
            float we1 = p.w[ie1];

            glm::vec2 edge = E1 - E0;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                continue;

            float t = glm::clamp(
                glm::dot(P - E0, edge) / edgeLengthSq,
                0.0f,
                1.0f);
            glm::vec2 closest = E0 + t * edge;
            glm::vec2 n = P - closest;
            float C = glm::length(n);
            if (C < 1e-6f)
                continue;
            n /= C;

            glm::vec2 gradP = n;
            glm::vec2 gradE0 = -n * (1.0f - t);
            glm::vec2 gradE1 = -n * t;

            float w =
                wp +
                we0 * (1.0f - t) * (1.0f - t) +
                we1 * t * t;
            if (w < 1e-6f)
                continue;

            float alpha = c.compliance / (dt * dt);
            float deltaLambdaN =
                (-C - alpha * c.lambdaN) / (w + alpha);
            c.lambdaN += deltaLambdaN;

            P += wp * deltaLambdaN * gradP;
            E0 += we0 * deltaLambdaN * gradE0;
            E1 += we1 * deltaLambdaN * gradE1;

            // Static friction
            glm::vec2 tangent = Utils::perpCCW2d(n);
            glm::vec2 dp = P - Pprev;
            glm::vec2 de0 = E0 - E0prev;
            glm::vec2 de1 = E1 - E1prev;

            glm::vec2 edgeDisp =
                de0 * (1.0f - t) +
                de1 * t;

            glm::vec2 relDisp = dp - edgeDisp;
            float Ct = glm::dot(relDisp, tangent);
            float deltaLambdaT = -Ct / w;
            float maxLambdaT = c.staticFriction * fabs(c.lambdaN);
            float oldLambdaT = c.lambdaT;

            c.lambdaT = glm::clamp(
                c.lambdaT + deltaLambdaT,
                -maxLambdaT,
                maxLambdaT);
            deltaLambdaT = c.lambdaT - oldLambdaT;

            P += wp * deltaLambdaT * tangent;
            E0 -= we0 * (1.0f - t) * deltaLambdaT * tangent;
            E1 -= we1 * t * deltaLambdaT * tangent;
        }
    }

    void apply_pointEdgeCollisionConstraints_kinetic_friction(
        Particles &p,
        std::vector<PointEdgeCollisionConstraints> &pecc,
        float dt)
    {
        for (auto &c : pecc)
        {
            size_t ip = c.point;
            size_t ie0 = c.edge1;
            size_t ie1 = c.edge2;

            glm::vec2 &P = p.pos[ip];
            glm::vec2 &E0 = p.pos[ie0];
            glm::vec2 &E1 = p.pos[ie1];

            glm::vec2 &Vp = p.vel[ip];
            glm::vec2 &Ve0 = p.vel[ie0];
            glm::vec2 &Ve1 = p.vel[ie1];

            float wp = p.w[ip];
            float we0 = p.w[ie0];
            float we1 = p.w[ie1];

            // Closest point on edge
            glm::vec2 edge = E1 - E0;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                continue;

            float t = glm::clamp(
                glm::dot(P - E0, edge) / edgeLengthSq,
                0.0f,
                1.0f);
            glm::vec2 closest = E0 + t * edge;
            glm::vec2 n = P - closest;
            float nLen = glm::length(n);
            if (nLen < 1e-6f)
                continue;
            n /= nLen;

            float w =
                wp +
                we0 * (1.0f - t) * (1.0f - t) +
                we1 * t * t;
            if (w < 1e-6f)
                continue;

            glm::vec2 edgeVelocity =
                Ve0 * (1.0f - t) +
                Ve1 * t;
            glm::vec2 relVelocity = Vp - edgeVelocity;
            glm::vec2 vt =
                relVelocity - glm::dot(relVelocity, n) * n;
            float vtLength = glm::length(vt);
            if (vtLength < 1e-6f)
                continue;
            glm::vec2 tangentVelocity = vt / vtLength;

            float maxImpulse =
                c.kineticFriction * fabs(c.lambdaN) / dt;

            float impulse =
                glm::min(vtLength / w, maxImpulse);
            glm::vec2 deltaV = -impulse * tangentVelocity;

            Vp += wp * deltaV;
            Ve0 -= we0 * (1.0f - t) * deltaV;
            Ve1 -= we1 * t * deltaV;
        }
    }

    bool get_angularData(
        const Particles &p,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        AngularData &data)
    {
        glm::vec2 dir1 = p.pos[i2] - p.pos[i1];
        glm::vec2 dir2 = p.pos[i4] - p.pos[i3];

        float len1Sq = glm::dot(dir1, dir1);
        float len2Sq = glm::dot(dir2, dir2);

        if (len1Sq < 1e-6f || len2Sq < 1e-6f)
            return false;

        data.grad1 = Utils::perpCCW2d(dir1) / len1Sq;
        data.grad2 = Utils::perpCCW2d(dir2) / len2Sq;

        data.invI1 =
            (p.w[i1] + p.w[i2]) *
            glm::dot(data.grad1, data.grad1);

        data.invI2 =
            (p.w[i3] + p.w[i4]) *
            glm::dot(data.grad2, data.grad2);

        data.invI = data.invI1 + data.invI2;

        if (data.invI < 1e-6f)
            return false;

        data.angle = Utils::relative_angle_rad(dir1, dir2);

        data.omega =
            glm::dot(p.vel[i4] - p.vel[i3], data.grad2) -
            glm::dot(p.vel[i2] - p.vel[i1], data.grad1);

        return true;
    }
    void apply_angular_velocity_correction(
        Particles &p,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        float deltaL)
    {
        AngularData data;

        if (!get_angularData(
                p, i1, i2, i3, i4, data))
            return;

        float lambda = deltaL / data.invI;

        glm::vec2 correction1 = data.grad1 * lambda;
        glm::vec2 correction2 = data.grad2 * lambda;

        p.vel[i1] += p.w[i1] * correction1;
        p.vel[i2] -= p.w[i2] * correction1;

        p.vel[i3] -= p.w[i3] * correction2;
        p.vel[i4] += p.w[i4] * correction2;
    }
    void apply_torque(
        Particles &p,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        float torque,
        float dt)
    {
        apply_angular_velocity_correction(
            p, i1, i2, i3, i4,
            torque * dt);
    }

    void solve_angularMotor(
        Particles &p,
        const AngularMotor &motor,
        float dt)
    {
        AngularData data;

        if (!get_angularData(
                p,
                motor.i1, motor.i2,
                motor.i3, motor.i4,
                data))
            return;

        // Желаемое изменение относительной угловой скорости
        float deltaOmega = motor.targetVelocity - data.omega;

        // Необходимый угловой импульс
        float deltaL = deltaOmega / data.invI;

        // Ограничение максимального момента
        float maxImpulse = motor.maxTorque * dt;

        deltaL = std::clamp(
            deltaL,
            -maxImpulse,
            maxImpulse);

        apply_angular_velocity_correction(
            p,
            motor.i1, motor.i2,
            motor.i3, motor.i4,
            deltaL);
    }
    void solve_angularServo(
        Particles &p,
        const AngularServo &servo,
        float dt)
    {
        AngularData data;

        if (!get_angularData(
                p,
                servo.i1, servo.i2,
                servo.i3, servo.i4,
                data))
            return;

        // Ошибка угла с учетом перехода через +/- PI
        float angleError = servo.targetAngle - data.angle;

        angleError = std::atan2(
            std::sin(angleError),
            std::cos(angleError));

        // PD-регулятор
        float torque =
            servo.kp * angleError -
            servo.kd * data.omega;

        torque = std::clamp(
            torque,
            -servo.maxTorque,
            servo.maxTorque);

        apply_torque(
            p,
            servo.i1, servo.i2,
            servo.i3, servo.i4,
            torque,
            dt);
    }

    // Class World
    World::World()
        : deltaTick(std::chrono::duration<double>(1.0 / tickRate))
    {
        updateThread = std::thread(&World::update_loop, this);
    }
    World::~World()
    {
        {
            std::lock_guard<std::mutex> lk(updateMutex);
            updateState = UpdateState::Stopped;
        }
        updateCv.notify_all();
        if (updateThread.joinable())
            updateThread.join();
    }

    void World::pause_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::Paused;
        updateCv.notify_all();
    }

    void World::resume_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::Running;
        nextUpdateTimePoint = std::chrono::steady_clock::now();
        updateCv.notify_all();
    }

    void World::step_once_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::StepOnce;
        nextUpdateTimePoint = std::chrono::steady_clock::now();
        updateCv.notify_all();
    }

    void World::stop_update()
    {
        std::lock_guard<std::mutex> lk(updateMutex);
        updateState = UpdateState::Stopped;
        updateCv.notify_all();
    }

    void World::update_loop()
    {
        std::unique_lock<std::mutex> update_lk(updateMutex);
        nextUpdateTimePoint = std::chrono::steady_clock::now();

        while (updateState != UpdateState::Stopped)
        {
            updateCv.wait(update_lk, [this]
                          { return updateState != UpdateState::Paused &&
                                   updateState != UpdateState::Stopped; });
            if (updateState == UpdateState::Stopped)
                break;

            update_lk.unlock();
            auto start = std::chrono::steady_clock::now();
            step();
            lastStepTime = std::chrono::steady_clock::now() - start;
            {
                // submit draw
                std::lock_guard<std::mutex> draw_lk(clusterDrawMutex);
                clusterDraw = cluster;
                newFrame = true;
            }
            execute_enqueued_commands();
            update_lk.lock();
            if (updateState == UpdateState::StepOnce)
            {
                updateState = UpdateState::Paused;
                continue;
            }

            nextUpdateTimePoint += deltaTick;
            updateCv.wait_until(update_lk, nextUpdateTimePoint, [this]
                                { return updateState != UpdateState::Running &&
                                         updateState != UpdateState::Stopped; });
        }
    }

    void World::step()
    {
        ZoneScoped;

        size_t polygons_size = cluster.polygonColliders.size();
        size_t points_size = cluster.pointsColliders.size();

        float substep_time = 1.0f / tickRate * timeScale / substeps;

        bool reversePass_switch = 0;

        for (size_t sub = 0; sub < substeps; ++sub)
        {
            iterate(cluster.particles, substep_time, gravity);

            reset_distanceConstraints_lambdas(cluster.distanceConstraints);
            reset_angleConstraints_lambdas(cluster.angleConstraints);
            reset_volumeConstraints_lambdas(cluster.volumeConstraints);

            // for (size_t iter = 0; iter < iterations; ++iter) ITERATION
            //{                                                 ITERATION
            if (reversePass)
                reversePass_switch = !reversePass_switch;

            solve_distanceConstraints(cluster.particles, cluster.distanceConstraints, substep_time, reversePass_switch);
            solve_angleConstraints(cluster.particles, cluster.angleConstraints, substep_time);
            solve_volumeConstraints(cluster.particles, cluster.volumeConstraints, substep_time);

            std::vector<AABB> polygonColliders_aabbs(polygons_size);
            generate_colliderPoints_aabbs(polygonColliders_aabbs, cluster.particles, cluster.polygonColliders);
            std::vector<AABB> pointsColliders_aabbs(points_size);
            generate_colliderPoints_aabbs(pointsColliders_aabbs, cluster.particles, cluster.pointsColliders);

            spatialHashAABB.clear();
            for (size_t i = 0; i < polygons_size; ++i)
                spatialHashAABB.add_aabb(polygonColliders_aabbs[i], i);

            std::vector<PolygonCache> polygonsCache(polygons_size);
            for (size_t polygonId = 0; polygonId < polygons_size; ++polygonId)
            {
                PolygonCache &hash = polygonsCache[polygonId];
                std::vector<size_t> &polygonIndices = cluster.polygonColliders[polygonId].indices;
                size_t n = polygonIndices.size();
                hash.positions_size = n;

                hash.positions.resize(polygonIndices.size());
                for (size_t i = 0; i < n; ++i)
                    hash.positions[i] = cluster.particles.pos[polygonIndices[i]];

                hash.edges.resize(n);
                hash.edgesLen2.resize(n);
                for (size_t i = 0, j = n - 1; i < n; j = i++)
                {
                    glm::vec2 edge = hash.positions[i] - hash.positions[j];
                    hash.edges[j] = edge;
                    hash.edgesLen2[j] = glm::dot(edge, edge);
                }
            }

            // until polygons_size: polygon/polygon // after: point/polygon
            std::vector<PointsCollider> combinedColliderPoints;
            combinedColliderPoints.reserve(polygons_size + points_size);
            combinedColliderPoints.insert(combinedColliderPoints.end(), cluster.polygonColliders.begin(), cluster.polygonColliders.end());
            combinedColliderPoints.insert(combinedColliderPoints.end(), cluster.pointsColliders.begin(), cluster.pointsColliders.end());

            std::vector<AABB> combinedAABB;
            combinedAABB.reserve(polygons_size + points_size);
            combinedAABB.insert(combinedAABB.end(), polygonColliders_aabbs.begin(), polygonColliders_aabbs.end());
            combinedAABB.insert(combinedAABB.end(), pointsColliders_aabbs.begin(), pointsColliders_aabbs.end());

            std::vector<std::vector<PointEdgeCollisionConstraints>> pecc_thread(omp_get_max_threads());

            // #pragma omp parallel for num_threads(6)
            for (size_t pid = 0; pid < combinedColliderPoints.size(); ++pid)
            {
                std::vector<size_t> overlapping_aabb_ids;
                spatialHashAABB.get_overlapping_aabb_ids_excludeId(overlapping_aabb_ids, combinedAABB[pid], pid);

                bool intersection_brute_force = false;
                if (intersection_brute_force)
                {
                    overlapping_aabb_ids.clear();
                    for (size_t i = 0; i < polygons_size; i++)
                        if (i != pid && combinedAABB[pid].intersects(polygonColliders_aabbs[i]))
                            overlapping_aabb_ids.emplace_back(i);
                }

                auto &local_pecc = pecc_thread[omp_get_thread_num()];
                local_pecc.reserve(overlapping_aabb_ids.size() * 10);

                for (size_t i : overlapping_aabb_ids)
                {
                    float avgStaticFriction = (combinedColliderPoints[pid].staticFriction + cluster.polygonColliders[i].staticFriction) * 0.5f;
                    float avgKineticFriction = (combinedColliderPoints[pid].kineticFriction + cluster.polygonColliders[i].kineticFriction) * 0.5f;
                    float avgCompliance = (combinedColliderPoints[pid].compliance + cluster.polygonColliders[i].compliance) * 0.5f;

                    PointsPolygonCollision pointsPolygonCollision{
                        combinedColliderPoints[pid].indices,
                        cluster.polygonColliders[i].indices,
                        avgStaticFriction,
                        avgKineticFriction,
                        avgCompliance,
                        combinedAABB[pid].get_intersection(polygonColliders_aabbs[i]),
                        polygonsCache[i],
                    };
                    add_PointEdgeCollisionConstraints_from_PointsPolygonCollision(local_pecc, cluster.particles, pointsPolygonCollision);
                }
            }

            // combine pecc from pecc_thread
            size_t total_size = 0;
            for (auto &v : pecc_thread)
                total_size += v.size();

            std::vector<PointEdgeCollisionConstraints> pecc;
            pecc.reserve(total_size);
            for (auto &v : pecc_thread)
                pecc.insert(pecc.end(), v.begin(), v.end());
            solve_pointEdgeCollisionConstraints(cluster.particles, pecc, substep_time);
            //} ITERATION

            update_velocities(cluster.particles, substep_time);
            apply_pointEdgeCollisionConstraints_kinetic_friction(cluster.particles, pecc, substep_time);
            // apply_distance_constraints_damping(particles, distanceConstraints, substep_time);

            for (auto &motor : cluster.angularMotors)
            {
                solve_angularMotor(
                    cluster.particles,
                    motor,
                    substep_time);
            }
            for (auto &servo : cluster.angularServos)
            {
                solve_angularServo(
                    cluster.particles,
                    servo,
                    substep_time);
            }
        }
    }

    void World::init()
    {
        cluster.particles.pos.clear();
        cluster.particles.pos.reserve(5000);
        cluster.particles.prev.clear();
        cluster.particles.prev.reserve(5000);
        cluster.particles.vel.clear();
        cluster.particles.vel.reserve(5000);
        cluster.particles.w.clear();
        cluster.particles.w.reserve(5000);

        cluster.distanceConstraints.clear();
        cluster.distanceConstraints.reserve(5000);
        cluster.angleConstraints.clear();
        cluster.angleConstraints.reserve(5000);
        cluster.volumeConstraints.clear();
        cluster.volumeConstraints.reserve(1000);

        cluster.polygonColliders.clear();
        cluster.polygonColliders.reserve(1000);
        cluster.pointsColliders.clear();
        cluster.pointsColliders.reserve(1000);

        cluster.angularMotors.clear();
        cluster.angularMotors.reserve(100);
        cluster.angularServos.clear();
        cluster.angularServos.reserve(100);

        // polygonColliders_aabbs.reserve(1000);
        // pointsColliders_aabbs.reserve(1000);

        // polygonsHash.reserve(1000);

        // combinedColliderPoints.reserve(2000);
        // combinedAABB.reserve(2000);

        // pecc.reserve(5000);
    }

    void World::set_tickRate(size_t tickRate)
    {
        deltaTick = std::chrono::duration<double>(1.0 / tickRate);
    }
    void World::spawn_from_json(const std::string &name, const glm::vec2 &position)
    {
        std::lock_guard<std::mutex> clusterAdd_lk(clusterAddMutex);
        json_body_loader::load(clusterAdd, name, position);

        std::lock_guard<std::mutex> cluster_lk(clusterMutex);
        merge_clusters(cluster, clusterAdd);
    }
    void World::spawn_polygon(glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
    {
        std::lock_guard<std::mutex> clusterAdd_lk(clusterAddMutex);
        xpbd::spawn_polygon(clusterAdd, pos, radius, segments, mass, compliance);

        std::lock_guard<std::mutex> cluster_lk(clusterMutex);
        merge_clusters(cluster, clusterAdd);
    }

    void World::enqueue_command(Command variant)
    {
        std::lock_guard<std::mutex> lk(commandQueueMutex);
        commandQueue.push(variant);
    }
    void World::execute_enqueued_commands()
    {
        std::queue<Command> local;
        {
            std::lock_guard<std::mutex> lk(commandQueueMutex);
            std::swap(local, commandQueue);
        }

        while (!local.empty())
        {
            Command cmd = std::move(local.front());
            local.pop();

            std::visit([this](auto &&arg)
                       {
                    using T = std::decay_t<decltype(arg)>;

                    if constexpr (std::is_same_v<T, init_command>)
                        init();
                    else if constexpr (std::is_same_v<T, set_tickRate_command>)
                        set_tickRate(arg.tickRate);
                    else if constexpr (std::is_same_v<T, add_particle_command>)
                        add_particle(cluster, arg.pos, arg.mass, arg.vel);
                    else if constexpr (std::is_same_v<T, add_distanceConstraint_command>)
                        add_distanceConstraint(cluster, arg.constraint.i1, arg.constraint.i2, arg.constraint.compliance, arg.constraint.restDist);
                    else if constexpr (std::is_same_v<T, add_angleConstraint_command>)
                        add_angleConstraint(cluster, arg.constraint.i1, arg.constraint.i2, arg.constraint.i3, arg.constraint.i4, arg.constraint.compliance, arg.constraint.minAngle, arg.constraint.maxAngle);
                    else if constexpr (std::is_same_v<T, add_volumeConstraint_command>)
                        add_volumeConstraint(cluster, arg.constraint.indices, arg.constraint.compliance, arg.constraint.restVolume);
                    
                    else if constexpr (std::is_same_v<T, add_polygonCollider_command>)
                        add_polygonCollider(cluster, arg.collider.indices, arg.collider.staticFriction, arg.collider.kineticFriction,
                                            arg.collider.compliance, arg.collider.layer, arg.collider.mask);
                    else if constexpr (std::is_same_v<T, add_pointsCollider_command>)
                        add_pointsCollider(cluster, arg.collider.indices, arg.collider.staticFriction, arg.collider.kineticFriction,
                                            arg.collider.compliance, arg.collider.layer, arg.collider.mask);

                    else if constexpr (std::is_same_v<T, add_angularMotor_command>)
                        add_angularMotor(cluster, arg.motor.i1, arg.motor.i2, arg.motor.i3, arg.motor.i4, arg.motor.targetVelocity, arg.motor.maxTorque);
                    else if constexpr (std::is_same_v<T, update_angularMotor_command>)
                        update_angularMotor(cluster, arg.index, arg.motor.targetVelocity, arg.motor.maxTorque);
                    
                    else if constexpr (std::is_same_v<T, add_angularServo_command>)
                        add_angularServo(cluster, arg.servo.i1, arg.servo.i2, arg.servo.i3, arg.servo.i4, arg.servo.targetAngle, arg.servo.kp, arg.servo.kd, arg.servo.maxTorque);            
                    else if constexpr (std::is_same_v<T, update_angularServo_command>)
                        update_angularServo(cluster, arg.index, arg.servo.targetAngle, arg.servo.kp, arg.servo.kd, arg.servo.maxTorque);
                        
                    else if constexpr (std::is_same_v<T, spawn_from_json_command>)
                        spawn_from_json(arg.name, arg.position);
                    else if constexpr (std::is_same_v<T, spawn_polygon_command>)
                        spawn_polygon(arg.pos, arg.radius, arg.segments, arg.mass, arg.compliance); }, cmd);
        }
    }
    void merge_clusters(Cluster &into, Cluster &from)
    {
        if (&into == &from)
            return;

        const size_t offset = into.particles.pos.size();

        // Particles
        into.particles.pos.insert(into.particles.pos.end(),
                                  from.particles.pos.begin(), from.particles.pos.end());
        into.particles.prev.insert(into.particles.prev.end(),
                                   from.particles.prev.begin(), from.particles.prev.end());
        into.particles.vel.insert(into.particles.vel.end(),
                                  from.particles.vel.begin(), from.particles.vel.end());
        into.particles.w.insert(into.particles.w.end(),
                                from.particles.w.begin(), from.particles.w.end());

        // Distance constraints
        into.distanceConstraints.reserve(into.distanceConstraints.size() + from.distanceConstraints.size());
        for (auto &dc : from.distanceConstraints)
        {
            dc.i1 += offset;
            dc.i2 += offset;
            into.distanceConstraints.push_back(dc);
        }

        // Angle constraints
        into.angleConstraints.reserve(into.angleConstraints.size() + from.angleConstraints.size());
        for (auto &ac : from.angleConstraints)
        {
            ac.i1 += offset;
            ac.i2 += offset;
            ac.i3 += offset;
            ac.i4 += offset;
            into.angleConstraints.push_back(ac);
        }

        // Volume constraints
        into.volumeConstraints.reserve(into.volumeConstraints.size() + from.volumeConstraints.size());
        for (auto &vc : from.volumeConstraints)
        {
            for (auto &idx : vc.indices)
                idx += offset;
            into.volumeConstraints.push_back(vc);
        }

        // Polygon colliders
        into.polygonColliders.reserve(into.polygonColliders.size() + from.polygonColliders.size());
        for (auto &pc : from.polygonColliders)
        {
            for (auto &idx : pc.indices)
                idx += offset;
            into.polygonColliders.push_back(pc);
        }

        // Points colliders
        into.pointsColliders.reserve(into.pointsColliders.size() + from.pointsColliders.size());
        for (auto &pc : from.pointsColliders)
        {
            for (auto &idx : pc.indices)
                idx += offset;
            into.pointsColliders.push_back(pc);
        }

        // Angular motors
        into.angularMotors.reserve(into.angularMotors.size() + from.angularMotors.size());
        for (auto &am : from.angularMotors)
        {
            am.i1 += offset;
            am.i2 += offset;
            am.i3 += offset;
            am.i4 += offset;
            into.angularMotors.push_back(am);
        }

        // Angular servos
        into.angularServos.reserve(into.angularServos.size() + from.angularServos.size());
        for (auto &as : from.angularServos)
        {
            as.i1 += offset;
            as.i2 += offset;
            as.i3 += offset;
            as.i4 += offset;
            into.angularServos.push_back(as);
        }

        from = Cluster();
    }

    void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass)
    {
        add_particle(cluster, pos, mass, {0, 0});
    }
    void add_particle(Cluster &cluster, const glm::vec2 &pos, float mass, const glm::vec2 &vel)
    {
        cluster.particles.pos.emplace_back(pos);
        cluster.particles.prev.emplace_back(pos);
        cluster.particles.vel.emplace_back(vel);
        float w = mass == 0 ? 0 : 1 / mass;
        cluster.particles.w.emplace_back(w);
    }

    void add_distanceConstraint(Cluster &cluster, size_t i1, size_t i2, float compliance)
    {
        add_distanceConstraint(cluster, i1, i2, compliance, glm::distance(cluster.particles.pos[i1], cluster.particles.pos[i2]));
    }
    void add_distanceConstraint(Cluster &cluster, size_t i1, size_t i2, float compliance, float restDist)
    {
        float compressionDamping = 0.5f;
        float extensionDamping = 0.5f;
        float lambda = 0;
        cluster.distanceConstraints.emplace_back(i1, i2, restDist, compressionDamping, extensionDamping, compliance, lambda);
    }

    void add_angleConstraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance)
    {
        glm::vec2 v1 = cluster.particles.pos[i1];
        glm::vec2 v2 = cluster.particles.pos[i2];
        glm::vec2 v3 = cluster.particles.pos[i3];
        glm::vec2 v4 = cluster.particles.pos[i4];

        float restAngle = Utils::relative_angle_rad(v2 - v1, v4 - v3);
        float lambda = 0;
        cluster.angleConstraints.emplace_back(i1, i2, i3, i4, restAngle, restAngle, compliance, lambda);
    }
    void add_angleConstraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float restAngle)
    {
        float lambda = 0;
        cluster.angleConstraints.emplace_back(i1, i2, i3, i4, restAngle, restAngle, compliance, lambda);
    }
    void add_angleConstraint(Cluster &cluster, size_t i1, size_t i2, size_t i3, size_t i4, float compliance, float minAngle, float maxAngle)
    {
        float lambda = 0;
        cluster.angleConstraints.emplace_back(i1, i2, i3, i4, minAngle, maxAngle, compliance, lambda);
    }

    void add_volumeConstraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance)
    {
        float restVolume = compute_polygon_area(cluster.particles, indices);
        float lambda = 0;
        cluster.volumeConstraints.emplace_back(indices, restVolume, compliance, lambda);
    }
    void add_volumeConstraint(Cluster &cluster, const std::vector<size_t> &indices, float compliance, float restVolume)
    {
        float lambda = 0;
        cluster.volumeConstraints.emplace_back(indices, restVolume, compliance, lambda);
    }
    void add_polygonCollider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask)
    {
        cluster.polygonColliders.emplace_back(indices, staticFriction, kineticFriction, compliance, layer, mask);
    }
    void add_pointsCollider(Cluster &cluster, const std::vector<size_t> &indices, float staticFriction, float kineticFriction, float compliance, int8_t layer, int8_t mask)
    {
        cluster.pointsColliders.emplace_back(indices, staticFriction, kineticFriction, compliance, layer, mask);
    }
    void add_angularMotor(
        Cluster &cluster,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        float targetVelocity,
        float maxTorque)
    {
        cluster.angularMotors.emplace_back(i1, i2, i3, i4, targetVelocity, maxTorque);
    }
    void update_angularMotor(
        Cluster &cluster,
        size_t index,
        float targetVelocity,
        float maxTorque)
    {
        if (index >= cluster.angularMotors.size())
            return;

        cluster.angularMotors[index].targetVelocity = targetVelocity;
        cluster.angularMotors[index].maxTorque = maxTorque;
    }
    void add_angularServo(
        Cluster &cluster,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        float targetAngle,
        float kp,
        float kd,
        float maxTorque)
    {
        cluster.angularServos.emplace_back(i1, i2, i3, i4, targetAngle, kp, kd, maxTorque);
    }
    void update_angularServo(
        Cluster &cluster,
        size_t index,
        float targetAngle,
        float kp,
        float kd,
        float maxTorque)
    {
        if (index >= cluster.angularServos.size())
            return;

        cluster.angularServos[index].targetAngle = targetAngle;
        cluster.angularServos[index].kp = kp;
        cluster.angularServos[index].kd = kd;
        cluster.angularServos[index].maxTorque = maxTorque;
    }

    void spawn_from_json(Cluster &cluster, const std::string &name, const glm::vec2 &position)
    {
        json_body_loader::load(cluster, name, position);
    }
    void spawn_polygon(Cluster &cluster, glm::vec2 pos, float radius, size_t segments, float mass, float compliance)
    {
        if (segments < 3)
            segments = 3;

        size_t start = cluster.particles.pos.size();
        size_t end = start + segments;
        std::vector<size_t> ids;

        float angleStep = 2.0f * float(M_PI) / segments;
        for (int i = 0; i < segments; ++i)
        {
            float angle = i * angleStep;
            glm::vec2 dir = glm::vec2(cosf(angle), sinf(angle));
            add_particle(cluster, dir * radius + pos, mass / segments);

            ids.emplace_back(start + i);
        }

        for (int i = 0; i < segments; ++i)
        {
            add_distanceConstraint(cluster, start + i, start + (i + 1) % segments, compliance);
            // add_angleConstraint(cluster, start + i, start + (i + 1) % segments, start + (i + 1) % segments, start + (i + 2) % segments, compliance);
        }

        add_volumeConstraint(cluster, ids, compliance);

        float staticFriction = 0.4f;
        float kineticFriction = 0.3f;

        int8_t layer = 1;
        int8_t mask = 1;

        add_polygonCollider(cluster, ids, staticFriction, kineticFriction, 0.0f, layer, mask);
    }

}
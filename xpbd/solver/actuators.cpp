#include "actuators.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/actuators.hpp"
#include "../utils.hpp"

#include "glm/glm.hpp"
#include <vector>

namespace xpbd
{
    bool get_angularData(
        const std::vector<Particle> &p,
        size_t i1, size_t i2,
        size_t i3, size_t i4,
        AngularData &data)
    {
        glm::vec2 dir1 = p[i2].pos - p[i1].pos;
        glm::vec2 dir2 = p[i4].pos - p[i3].pos;

        float len1Sq = glm::dot(dir1, dir1);
        float len2Sq = glm::dot(dir2, dir2);

        if (len1Sq < 1e-6f || len2Sq < 1e-6f)
            return false;

        data.grad1 = Utils::perpCCW2d(dir1) / len1Sq;
        data.grad2 = Utils::perpCCW2d(dir2) / len2Sq;

        data.invI1 =
            (p[i1].w + p[i2].w) *
            glm::dot(data.grad1, data.grad1);

        data.invI2 =
            (p[i3].w + p[i4].w) *
            glm::dot(data.grad2, data.grad2);

        data.invI = data.invI1 + data.invI2;

        if (data.invI < 1e-6f)
            return false;

        data.angle = Utils::relative_angle_rad(dir1, dir2);

        data.omega =
            glm::dot(p[i4].vel - p[i3].vel, data.grad2) -
            glm::dot(p[i2].vel - p[i1].vel, data.grad1);

        return true;
    }
    void apply_angular_velocity_correction(
        std::vector<Particle> &p,
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

        p[i1].vel += p[i1].w * correction1;
        p[i2].vel -= p[i2].w * correction1;

        p[i3].vel -= p[i3].w * correction2;
        p[i4].vel += p[i4].w * correction2;
    }
    void apply_torque(
        std::vector<Particle> &p,
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
        std::vector<Particle> &p,
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

        constexpr float RPM_TO_RAD = 2.0f * 3.14159265359f / 60.0f;

        float targetOmega = motor.targetRPM * RPM_TO_RAD;

        float deltaOmega = targetOmega - data.omega;
        float deltaL = deltaOmega / data.invI;
        float maxImpulse = motor.maxTorque * dt;

        deltaL = glm::clamp(
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
        std::vector<Particle> &p,
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

        float angleError = servo.targetAngle - data.angle;

        angleError = std::atan2(
            std::sin(angleError),
            std::cos(angleError));

        float torque =
            servo.kp * angleError -
            servo.kd * data.omega;

        torque = glm::clamp(
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
} // namespace xpbd
#include "actuators.hpp"
#include "xpbd/core/particles.hpp"
#include "xpbd/core/actuators.hpp"
#include "../utils.hpp"
#include "glm/glm.hpp"

namespace xpbd
{
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
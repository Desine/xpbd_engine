#include "actuators.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/actuators.hpp"
#include "xpbd/core/storage/component_storage.hpp"
#include "../utils.hpp"

#include "glm/glm.hpp"
#include <cmath>

namespace xpbd
{
bool get_angularData(
    const ParticleStorage &p,
    ParticleId i1, ParticleId i2,
    ParticleId i3, ParticleId i4,
    AngularData &data)
{
    const Particle *p1 = p.get(i1);
    const Particle *p2 = p.get(i2);
    const Particle *p3 = p.get(i3);
    const Particle *p4 = p.get(i4);

    if (!p1 || !p2 || !p3 || !p4)
        return false;

    glm::vec2 dir1 = p2->pos - p1->pos;
    glm::vec2 dir2 = p4->pos - p3->pos;

    float len1Sq = glm::dot(dir1, dir1);
    float len2Sq = glm::dot(dir2, dir2);

    if (len1Sq < 1e-6f || len2Sq < 1e-6f)
        return false;

    data.grad1 = Utils::perpCCW2d(dir1) / len1Sq;
    data.grad2 = Utils::perpCCW2d(dir2) / len2Sq;

    data.invI1 =
        (p1->w + p2->w) *
        glm::dot(data.grad1, data.grad1);

    data.invI2 =
        (p3->w + p4->w) *
        glm::dot(data.grad2, data.grad2);

    data.invI = data.invI1 + data.invI2;

    if (data.invI < 1e-6f)
        return false;

    data.angle = Utils::relative_angle_rad(dir1, dir2);

    data.omega =
        glm::dot(p4->vel - p3->vel, data.grad2) -
        glm::dot(p2->vel - p1->vel, data.grad1);

    return true;
}

void apply_angular_velocity_correction(
    ParticleStorage &p,
    ParticleId i1, ParticleId i2,
    ParticleId i3, ParticleId i4,
    float deltaL)
{
    AngularData data;

    if (!get_angularData(p, i1, i2, i3, i4, data))
        return;

    Particle *p1 = p.get(i1);
    Particle *p2 = p.get(i2);
    Particle *p3 = p.get(i3);
    Particle *p4 = p.get(i4);

    if (!p1 || !p2 || !p3 || !p4)
        return;

    float lambda = deltaL / data.invI;

    glm::vec2 correction1 = data.grad1 * lambda;
    glm::vec2 correction2 = data.grad2 * lambda;

    p1->vel += p1->w * correction1;
    p2->vel -= p2->w * correction1;

    p3->vel -= p3->w * correction2;
    p4->vel += p4->w * correction2;
}

void apply_torque(
    ParticleStorage &p,
    ParticleId i1, ParticleId i2,
    ParticleId i3, ParticleId i4,
    float torque,
    float dt)
{
    apply_angular_velocity_correction(
        p, i1, i2, i3, i4,
        torque * dt);
}

void solve_angularMotor(
    ParticleStorage &p,
    const AngularMotor &motor,
    float dt)
{
    AngularData data;

    if (!get_angularData(
            p,
            motor.id1, motor.id2,
            motor.id3, motor.id4,
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
        motor.id1, motor.id2,
        motor.id3, motor.id4,
        deltaL);
}

void solve_angularServo(
    ParticleStorage &p,
    const AngularServo &servo,
    float dt)
{
    AngularData data;

    if (!get_angularData(
            p,
            servo.id1, servo.id2,
            servo.id3, servo.id4,
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
        servo.id1, servo.id2,
        servo.id3, servo.id4,
        torque,
        dt);
}
} // namespace xpbd
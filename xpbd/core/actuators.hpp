#pragma once
#include "storage/handles.hpp"

#include "glm/glm.hpp"

namespace xpbd
{
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
        ParticleId id1, id2;
        ParticleId id3, id4;
        float targetRPM;
        float maxTorque;
    };
    struct AngularServo
    {
        ParticleId id1, id2;
        ParticleId id3, id4;
        float targetAngle;
        float kp;
        float kd;
        float maxTorque;
    };
} // namespace xpbd
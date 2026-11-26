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
        float target_rpm;
        float max_torque;
    };
    struct AngularServo
    {
        ParticleId id1, id2;
        ParticleId id3, id4;
        float target_angle;
        float kp;
        float kd;
        float max_torque;
    };
} // namespace xpbd
#pragma once
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
} // namespace xpbd
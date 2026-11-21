#pragma once
#include "xpbd/core/particles.hpp"
#include "xpbd/core/actuators.hpp"
#include <vector>

namespace xpbd
{
    bool get_angularData(const std::vector<Particle> &p, size_t i1, size_t i2, size_t i3, size_t i4, AngularData &data);
    void apply_angular_velocity_correction(std::vector<Particle> &p, size_t i1, size_t i2, size_t i3, size_t i4, float deltaL);
    void apply_torque(std::vector<Particle> &p, size_t i1, size_t i2, size_t i3, size_t i4, float torque, float dt);
    void solve_angularMotor(std::vector<Particle> &p, const AngularMotor &motor, float dt);
    void solve_angularServo(std::vector<Particle> &p, const AngularServo &servo, float dt);
} // namespace xpbd
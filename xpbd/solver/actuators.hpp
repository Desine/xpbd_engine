#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/actuators.hpp"
#include "xpbd/core/storage/component_storage.hpp"

namespace xpbd
{
bool get_angularData(const ParticleStorage &p, ParticleId i1, ParticleId i2, ParticleId i3, ParticleId i4, AngularData &data);
void apply_angular_velocity_correction(ParticleStorage &p, ParticleId i1, ParticleId i2, ParticleId i3, ParticleId i4, float deltaL);
void apply_torque(ParticleStorage &p, ParticleId i1, ParticleId i2, ParticleId i3, ParticleId i4, float torque, float dt);
void solve_angularMotor(ParticleStorage &p, const AngularMotor &motor, float dt);
void solve_angularServo(ParticleStorage &p, const AngularServo &servo, float dt);
} // namespace xpbd
#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "xpbd/core/storage/component_storage.hpp"

#include <vector>

namespace xpbd
{
void solve_angle_constraints(ParticleStorage &p, AngleConstraintStorage &ac, float dt);
void solve_angle_constraint(ParticleStorage &p, AngleConstraint &ac, float dt);
void reset_angle_constraints_lambdas(AngleConstraintStorage &ac);
} // namespace xpbd
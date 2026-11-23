#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "xpbd/core/storage/component_storage.hpp"

#include <vector>

namespace xpbd
{
void solve_distance_constraints(ParticleStorage &p, DistanceConstraintStorage &dc, float dt);
void solve_distance_constraint(ParticleStorage &p, DistanceConstraint &dc, float dt);
void apply_distance_constraints_damping(ParticleStorage &p, DistanceConstraintStorage &dc, float dt);
void reset_distance_constraints_lambdas(DistanceConstraintStorage &dc);
} // namespace xpbd
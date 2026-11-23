#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "xpbd/core/storage/component_storage.hpp"

#include <vector>

namespace xpbd
{
float polygon_area(const ParticleStorage &storage, const std::vector<ParticleId> &ids);
void solve_volume_constraints(ParticleStorage &p, VolumeConstraintStorage &vc, float dt);
void solve_volume_constraint(ParticleStorage &p, VolumeConstraint &vc, float dt);
void reset_volume_constraints_lambdas(VolumeConstraintStorage &vc);
} // namespace xpbd
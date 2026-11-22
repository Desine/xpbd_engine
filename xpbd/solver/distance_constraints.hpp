#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"

#include <vector>

namespace xpbd
{
    void solve_distance_constraint(std::vector<Particle> &p, DistanceConstraint &dc, float dt);
    void solve_distance_constraints(std::vector<Particle> &p, std::vector<DistanceConstraint> &dc, float dt, bool reversePass);
    void apply_distance_constraints_damping(std::vector<Particle> &p, std::vector<DistanceConstraint> &dc, float dt);
    void reset_distance_constraints_lambdas(std::vector<DistanceConstraint> &dc);
} // namespace xpbd
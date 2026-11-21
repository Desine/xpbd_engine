#pragma once
#include "angle_constraints.hpp"
#include "xpbd/core/particles.hpp"
#include "xpbd/core/constraints.hpp"
#include <vector>

namespace xpbd
{
    void solve_angle_constraint(std::vector<Particle> &p, AngleConstraint &ac, float dt);
    void solve_angle_constraints(std::vector<Particle> &p, std::vector<AngleConstraint> &ac, float dt);
    void reset_angle_constraints_lambdas(std::vector<AngleConstraint> &ac);
} // namespace xpbd
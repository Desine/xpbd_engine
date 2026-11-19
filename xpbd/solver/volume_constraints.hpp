#pragma once
#include "xpbd/core/particles.hpp"
#include "xpbd/core/constraints.hpp"
#include <vector>

namespace xpbd
{
    float compute_polygon_area(const std::vector<glm::vec2> &positions);
    float compute_polygon_area(const Particles &p, const std::vector<size_t> &indices);
    void solve_volume_constraint(Particles &p, VolumeConstraint &vc, float dt);
    void solve_volume_constraints(Particles &p, std::vector<VolumeConstraint> &vc, float dt);
    void reset_volume_constraints_lambdas(std::vector<VolumeConstraint> &vc);
} // namespace xpbd
#pragma once
#include "xpbd/core/particles.hpp"
#include "xpbd/core/colliders.hpp"
#include <vector>

namespace xpbd
{
    void generate_collider_points_aabbs(std::vector<AABB> &out, const Particles &p, const std::vector<PointsCollider> &colliderPoints);

    void add_point_edge_collision_constraints_from_points_polygon_collision(std::vector<PointEdgeCollisionConstraints> &pecc, const Particles &p, const PointsPolygonCollision &collision);
    void solve_point_edge_collision_constraints(Particles &p, std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
    void apply_point_edge_collision_constraints_kinetic_friction(Particles &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
} // namespace xpbd

#pragma once
#include "xpbd/core/particle.hpp"
#include "xpbd/core/colliders.hpp"
#include "xpbd/core/storage/component_storage.hpp"

#include <vector>

namespace xpbd
{
void generate_collider_points_aabbs(std::vector<AABB> &out, const ParticleStorage &p, const PointsColliderStorage &pointsCollider);
void generate_collider_points_aabbs(std::vector<AABB> &out, const ParticleStorage &p, const PolygonColliderStorage &polygonCollider);
void add_point_edge_collision_constraints_from_points_polygon_collision(std::vector<PointEdgeCollisionConstraints> &out, const ParticleStorage &particles, const PointsPolygonCollision &collision);
void solve_point_edge_collision_constraints(ParticleStorage &p, std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
void apply_point_edge_collision_constraints_kinetic_friction(ParticleStorage &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt);
} // namespace xpbd

#pragma once
#include "storage/handles.hpp"
#include "xpbd/solver/spatial_hashing.hpp"

#include "glm/glm.hpp"
#include <vector>

namespace xpbd
{
struct PointsCollider
{
    std::vector<ParticleId> indices;
    float static_friction;
    float kinetic_friction;
    float compliance;
    uint32_t layer;
    uint32_t mask;
};
struct PointEdgeCollision
{
    ParticleId point;
    ParticleId edge1;
    ParticleId edge2;
};
struct PolygonCache
{
    std::vector<glm::vec2> positions;
    size_t positions_size;
    std::vector<glm::vec2> edges;
    std::vector<float> edgesLen2;
};
struct PointsPolygonCollision
{
    std::vector<ParticleId> points;
    std::vector<ParticleId> polygon;
    float static_friction;
    float kinetic_friction;
    float compliance;
    AABB box;
    PolygonCache cache;
};
struct PointEdgeCollisionConstraints
{
    ParticleId point;
    ParticleId edge1;
    ParticleId edge2;
    float static_friction;
    float kinetic_friction;
    float compliance;
    float lambdaN = 0.0f;
    float lambdaT = 0.0f;
};
} // namespace xpbd
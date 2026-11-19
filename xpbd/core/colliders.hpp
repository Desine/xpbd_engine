#pragma once
#include "glm/glm.hpp"
#include <vector>
#include "types.hpp"
namespace xpbd
{
    struct PointsCollider
    {
        std::vector<size_t> indices;
        float staticFriction;
        float kineticFriction;
        float compliance;
        int8_t layer;
        int8_t mask;
    };
    struct PointEdgeCollision
    {
        size_t point;
        size_t edge1;
        size_t edge2;
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
        std::vector<size_t> points;
        std::vector<size_t> polygon;
        float staticFriction;
        float kineticFriction;
        float compliance;
        AABB box;
        PolygonCache cache;
    };
    struct PointEdgeCollisionConstraints
    {
        size_t point;
        size_t edge1;
        size_t edge2;
        float staticFriction;
        float kineticFriction;
        float compliance;
        float lambdaN = 0.0f;
        float lambdaT = 0.0f;
    };
} // namespace xpbd
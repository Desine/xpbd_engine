#include "collision.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/colliders.hpp"
#include "../utils.hpp"

#include <vector>

namespace xpbd
{
    void generate_collider_points_aabbs(
        std::vector<AABB> &out,
        const std::vector<Particle> &p,
        const std::vector<PointsCollider> &colliderPoints)
    {
        out.clear();
        out.reserve(colliderPoints.size());

        for (auto &c : colliderPoints)
        {
            AABB current;
            current.l = std::numeric_limits<float>::max();
            current.r = -std::numeric_limits<float>::max();
            current.b = std::numeric_limits<float>::max();
            current.t = -std::numeric_limits<float>::max();
            current.layer = c.layer;
            current.mask = c.mask;

            for (auto i : c.indices)
            {
                glm::vec2 pos = p[i].pos;
                if (pos.x < current.l)
                    current.l = pos.x;
                if (pos.x > current.r)
                    current.r = pos.x;
                if (pos.y < current.b)
                    current.b = pos.y;
                if (pos.y > current.t)
                    current.t = pos.y;
            }

            out.emplace_back(current);
        }
    }

    bool point_in_polygon(const glm::vec2 &point, const std::vector<glm::vec2> &positions)
    {
        int windingNumber = 0;
        const size_t n = positions.size();

        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &v1 = positions[i];
            const glm::vec2 &v2 = positions[j];
            float cross = Utils::cross(v2 - v1, point - v1);

            if (v1.y <= point.y)
            {
                if (v2.y > point.y && cross > 0)
                    ++windingNumber;
            }
            else
            {
                if (v2.y <= point.y && cross < 0)
                    --windingNumber;
            }
        }
        return windingNumber != 0;
    }
    void get_point_edge_collisions_of_points_inside_polygon(
        std::vector<PointEdgeCollision> &out,
        const std::vector<Particle> &p,
        const std::vector<size_t> &pointIndices,
        const std::vector<size_t> &polygonIndices,
        const PolygonCache &cache)
    {
        out.reserve(pointIndices.size());

        const auto &polygonPositions = cache.positions;
        const auto &edges = cache.edges;
        const auto &edgesLen2 = cache.edgesLen2;
        size_t n = cache.positions_size;

        for (size_t pid : pointIndices)
        {
            const glm::vec2 &point = p[pid].pos;

            if (!point_in_polygon(point, polygonPositions))
                continue;

            float minDist = std::numeric_limits<float>::max();
            size_t nearestEdge = 0;

            for (size_t i = 0; i < n; ++i)
            {
                if (cache.edgesLen2[i] < 1e-12f)
                    continue;

                glm::vec2 toPoint = point - polygonPositions[i];
                float t = glm::dot(toPoint, edges[i]) / edgesLen2[i];
                t = glm::clamp(t, 0.0f, 1.0f);
                glm::vec2 closest = polygonPositions[i] + edges[i] * t;

                float dist = glm::dot(point - closest, point - closest);
                if (dist < minDist)
                {
                    minDist = dist;
                    nearestEdge = i;
                }
            }

            out.emplace_back(pid,
                             polygonIndices[nearestEdge],
                             polygonIndices[(nearestEdge + 1) % n]);
        }
    }
    void add_point_edge_collision_constraints_from_points_polygon_collision(
        std::vector<PointEdgeCollisionConstraints> &pecc,
        const std::vector<Particle> &particles,
        const PointsPolygonCollision &collision)
    {
        std::vector<size_t> filteredPointIndices;
        filteredPointIndices.reserve(collision.points.size());
        for (auto i : collision.points)
            if (collision.box.contains_point(particles[i].pos))
                filteredPointIndices.emplace_back(i);

        if (filteredPointIndices.empty())
            return;

        std::vector<PointEdgeCollision> detections;
        detections.reserve(100);
        get_point_edge_collisions_of_points_inside_polygon(detections, particles, filteredPointIndices, collision.polygon, collision.cache);

        pecc.reserve(pecc.size() + detections.size());
        for (const auto &c : detections)
            pecc.emplace_back(c.point, c.edge1, c.edge2, collision.staticFriction, collision.kineticFriction, collision.compliance, 0.0f, 0.0f);
    }

    void solve_point_edge_collision_constraints(
        std::vector<Particle> &p,
        std::vector<PointEdgeCollisionConstraints> &pecc,
        float dt)
    {
        for (auto &c : pecc)
        {
            size_t ip = c.point;
            size_t ie0 = c.edge1;
            size_t ie1 = c.edge2;

            glm::vec2 &P = p[ip].pos;
            glm::vec2 &Pprev = p[ip].prev;
            glm::vec2 &E0 = p[ie0].pos;
            glm::vec2 &E1 = p[ie1].pos;
            glm::vec2 &E0prev = p[ie0].prev;
            glm::vec2 &E1prev = p[ie1].prev;

            float wp = p[ip].w;
            float we0 = p[ie0].w;
            float we1 = p[ie1].w;

            glm::vec2 edge = E1 - E0;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                continue;

            float t = glm::clamp(
                glm::dot(P - E0, edge) / edgeLengthSq,
                0.0f,
                1.0f);
            glm::vec2 closest = E0 + t * edge;
            glm::vec2 n = P - closest;
            float C = glm::length(n);
            if (C < 1e-6f)
                continue;
            n /= C;

            glm::vec2 gradP = n;
            glm::vec2 gradE0 = -n * (1.0f - t);
            glm::vec2 gradE1 = -n * t;

            float w =
                wp +
                we0 * (1.0f - t) * (1.0f - t) +
                we1 * t * t;
            if (w < 1e-6f)
                continue;

            float alpha = c.compliance / (dt * dt);
            float deltaLambdaN =
                (-C - alpha * c.lambdaN) / (w + alpha);
            c.lambdaN += deltaLambdaN;

            P += wp * deltaLambdaN * gradP;
            E0 += we0 * deltaLambdaN * gradE0;
            E1 += we1 * deltaLambdaN * gradE1;

            // Static friction
            glm::vec2 tangent = Utils::perpCCW2d(n);
            glm::vec2 dp = P - Pprev;
            glm::vec2 de0 = E0 - E0prev;
            glm::vec2 de1 = E1 - E1prev;

            glm::vec2 edgeDisp =
                de0 * (1.0f - t) +
                de1 * t;

            glm::vec2 relDisp = dp - edgeDisp;
            float Ct = glm::dot(relDisp, tangent);
            float deltaLambdaT = -Ct / w;
            float maxLambdaT = c.staticFriction * fabs(c.lambdaN);
            float oldLambdaT = c.lambdaT;

            c.lambdaT = glm::clamp(
                c.lambdaT + deltaLambdaT,
                -maxLambdaT,
                maxLambdaT);
            deltaLambdaT = c.lambdaT - oldLambdaT;

            P += wp * deltaLambdaT * tangent;
            E0 -= we0 * (1.0f - t) * deltaLambdaT * tangent;
            E1 -= we1 * t * deltaLambdaT * tangent;
        }
    }

    void apply_point_edge_collision_constraints_kinetic_friction(std::vector<Particle> &p, const std::vector<PointEdgeCollisionConstraints> &pecc, float dt)
    {
        for (auto &c : pecc)
        {
            size_t ip = c.point;
            size_t ie0 = c.edge1;
            size_t ie1 = c.edge2;

            glm::vec2 &P = p[ip].pos;
            glm::vec2 &E0 = p[ie0].pos;
            glm::vec2 &E1 = p[ie1].pos;

            glm::vec2 &Vp = p[ip].vel;
            glm::vec2 &Ve0 = p[ie0].vel;
            glm::vec2 &Ve1 = p[ie1].vel;

            float wp = p[ip].w;
            float we0 = p[ie0].w;
            float we1 = p[ie1].w;

            // Closest point on edge
            glm::vec2 edge = E1 - E0;
            float edgeLengthSq = glm::dot(edge, edge);
            if (edgeLengthSq < 1e-6f)
                continue;

            float t = glm::clamp(
                glm::dot(P - E0, edge) / edgeLengthSq,
                0.0f,
                1.0f);
            glm::vec2 closest = E0 + t * edge;
            glm::vec2 n = P - closest;
            float nLen = glm::length(n);
            if (nLen < 1e-6f)
                continue;
            n /= nLen;

            float w =
                wp +
                we0 * (1.0f - t) * (1.0f - t) +
                we1 * t * t;
            if (w < 1e-6f)
                continue;

            glm::vec2 edgeVelocity =
                Ve0 * (1.0f - t) +
                Ve1 * t;
            glm::vec2 relVelocity = Vp - edgeVelocity;
            glm::vec2 vt =
                relVelocity - glm::dot(relVelocity, n) * n;
            float vtLength = glm::length(vt);
            if (vtLength < 1e-6f)
                continue;
            glm::vec2 tangentVelocity = vt / vtLength;

            float maxImpulse =
                c.kineticFriction * fabs(c.lambdaN) / dt;

            float impulse =
                glm::min(vtLength / w, maxImpulse);
            glm::vec2 deltaV = -impulse * tangentVelocity;

            Vp += wp * deltaV;
            Ve0 -= we0 * (1.0f - t) * deltaV;
            Ve1 -= we1 * t * deltaV;
        }
    }
} // namespace xpbd

#include "collision.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/colliders.hpp"
#include "xpbd/core/storage/component_storage.hpp"
#include "../utils.hpp"

#include <vector>
#include <limits>
#include <cmath>

namespace xpbd
{
void generate_collider_points_aabbs(
    std::vector<AABB> &out,
    const ParticleStorage &p,
    const PointsColliderStorage &pointsCollider)
{
    out.clear();
    out.reserve(pointsCollider.size());

    for (size_t i = 0; i < pointsCollider.size(); ++i)
    {
        const PointsCollider *c = pointsCollider.at(i);
        if (!c)
            continue;

        AABB current;
        current.l = std::numeric_limits<float>::max();
        current.r = -std::numeric_limits<float>::max();
        current.b = std::numeric_limits<float>::max();
        current.t = -std::numeric_limits<float>::max();
        current.layer = c->layer;
        current.mask = c->mask;

        for (ParticleId id : c->indices)
        {
            const Particle *particle = p.get(id);
            if (!particle)
                continue;

            const glm::vec2 &pos = particle->pos;

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
void generate_collider_points_aabbs(
    std::vector<AABB> &out,
    const ParticleStorage &p,
    const PolygonColliderStorage &polygonCollider)
{
    out.clear();
    out.reserve(polygonCollider.size());

    for (size_t i = 0; i < polygonCollider.size(); ++i)
    {
        const PointsCollider *c = polygonCollider.at(i);
        if (!c)
            continue;

        AABB current;
        current.l = std::numeric_limits<float>::max();
        current.r = -std::numeric_limits<float>::max();
        current.b = std::numeric_limits<float>::max();
        current.t = -std::numeric_limits<float>::max();
        current.layer = c->layer;
        current.mask = c->mask;

        for (ParticleId id : c->indices)
        {
            const Particle *particle = p.get(id);
            if (!particle)
                continue;

            const glm::vec2 &pos = particle->pos;

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

void get_point_edge_collisions_of_points_inside_polygon(
    std::vector<PointEdgeCollision> &out,
    const ParticleStorage &p,
    const std::vector<ParticleId> &pointIds,
    const std::vector<ParticleId> &polygonIds,
    const PolygonCache &cache)
{
    out.reserve(out.size() + pointIds.size());

    const auto &polygonPositions = cache.positions;
    const auto &edges = cache.edges;
    const auto &edgesLen2 = cache.edgesLen2;
    size_t n = cache.positions_size;

    if (n == 0 || polygonIds.size() < n)
        return;

    for (ParticleId pid : pointIds)
    {
        const Particle *particle = p.get(pid);
        if (!particle)
            continue;

        const glm::vec2 &point = particle->pos;

        if (!Utils::point_in_polygon(point, polygonPositions))
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

        out.emplace_back(
            pid,
            polygonIds[nearestEdge],
            polygonIds[(nearestEdge + 1) % n]);
    }
}

void add_point_edge_collision_constraints_from_points_polygon_collision(
    std::vector<PointEdgeCollisionConstraints> &pecc,
    const ParticleStorage &particles,
    const PointsPolygonCollision &collision)
{
    std::vector<ParticleId> filteredPointIds;
    filteredPointIds.reserve(collision.points.size());

    for (ParticleId id : collision.points)
    {
        const Particle *particle = particles.get(id);
        if (particle && collision.box.contains_point(particle->pos))
            filteredPointIds.emplace_back(id);
    }

    if (filteredPointIds.empty())
        return;

    std::vector<PointEdgeCollision> detections;
    detections.reserve(100);

    get_point_edge_collisions_of_points_inside_polygon(
        detections,
        particles,
        filteredPointIds,
        collision.polygon,
        collision.cache);

    pecc.reserve(pecc.size() + detections.size());

    for (const auto &c : detections)
    {
        pecc.emplace_back(
            c.point,
            c.edge1,
            c.edge2,
            collision.staticFriction,
            collision.kineticFriction,
            collision.compliance,
            0.0f,
            0.0f);
    }
}

void solve_point_edge_collision_constraints(
    ParticleStorage &p,
    std::vector<PointEdgeCollisionConstraints> &pecc,
    float dt)
{
    for (auto &c : pecc)
    {
        Particle *particleP = p.get(c.point);
        Particle *particleE0 = p.get(c.edge1);
        Particle *particleE1 = p.get(c.edge2);

        if (!particleP || !particleE0 || !particleE1)
            continue;

        glm::vec2 &P = particleP->pos;
        glm::vec2 &Pprev = particleP->prev;
        glm::vec2 &E0 = particleE0->pos;
        glm::vec2 &E1 = particleE1->pos;
        glm::vec2 &E0prev = particleE0->prev;
        glm::vec2 &E1prev = particleE1->prev;

        float wp = particleP->w;
        float we0 = particleE0->w;
        float we1 = particleE1->w;

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

void apply_point_edge_collision_constraints_kinetic_friction(
    ParticleStorage &p,
    const std::vector<PointEdgeCollisionConstraints> &pecc,
    float dt)
{
    for (const auto &c : pecc)
    {
        Particle *particleP = p.get(c.point);
        Particle *particleE0 = p.get(c.edge1);
        Particle *particleE1 = p.get(c.edge2);

        if (!particleP || !particleE0 || !particleE1)
            continue;

        glm::vec2 &P = particleP->pos;
        glm::vec2 &E0 = particleE0->pos;
        glm::vec2 &E1 = particleE1->pos;

        glm::vec2 &Vp = particleP->vel;
        glm::vec2 &Ve0 = particleE0->vel;
        glm::vec2 &Ve1 = particleE1->vel;

        float wp = particleP->w;
        float we0 = particleE0->w;
        float we1 = particleE1->w;

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
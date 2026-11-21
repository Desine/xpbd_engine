#pragma once
#include "entity.hpp"
#include "particles.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"
#include <vector>
#include <string>

namespace xpbd
{
    struct Cluster
    {
        std::vector<Entity> entities;

        std::vector<std::pair<std::string, size_t>> particleIds;
        std::vector<Particle> particles;
        std::vector<DistanceConstraint> distanceConstraints;
        std::vector<AngleConstraint> angleConstraints;
        std::vector<VolumeConstraint> volumeConstraints;
        std::vector<PointsCollider> polygonColliders;
        std::vector<PointsCollider> pointsColliders;
        std::vector<AngularMotor> angularMotors;
        std::vector<AngularServo> angularServos;
    };

    void merge_clusters(Cluster &into, Cluster &from);
} // namespace xpbd
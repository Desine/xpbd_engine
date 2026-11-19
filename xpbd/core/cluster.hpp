#pragma once
#include <vector>
#include <string>
#include "particles.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"
namespace xpbd
{
    struct Cluster
    {
        std::vector<std::pair<std::string, size_t>> particleIds;
        Particles particles;
        std::vector<DistanceConstraint> distanceConstraints;
        std::vector<AngleConstraint> angleConstraints;
        std::vector<VolumeConstraint> volumeConstraints;
        std::vector<PointsCollider> polygonColliders;
        std::vector<PointsCollider> pointsColliders;
        std::vector<AngularMotor> angularMotors;
        std::vector<AngularServo> angularServos;
    };
} // namespace xpbd
#pragma once
#include "handles.hpp"
#include "dense_storage.hpp"
#include "../constraints.hpp"
#include "../colliders.hpp"
#include "../actuators.hpp"
#include "../particle.hpp"
#include "../entity.hpp"

namespace xpbd
{
    using EntityStorage = DenseStorage<Entity, EntityTag>;

    using ParticleStorage = DenseStorage<Particle, ParticleTag>;

    using DistanceConstraintStorage = DenseStorage<DistanceConstraint, DistanceConstraintTag>;
    using AngleConstraintStorage = DenseStorage<AngleConstraint, AngleConstraintTag>;
    using VolumeConstraintStorage = DenseStorage<VolumeConstraint, VolumeConstraintTag>;

    using PointsColliderStorage = DenseStorage<PointsCollider, PointsColliderTag>;
    using PolygonColliderStorage = DenseStorage<PointsCollider, PolygonColliderTag>;

    using AngularMotorStorage = DenseStorage<AngularMotor, AngularMotorTag>;
    using AngularServoStorage = DenseStorage<AngularServo, AngularServoTag>;

} // namespace xpbd
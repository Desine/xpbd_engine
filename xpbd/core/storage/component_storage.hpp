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
    using DistanceConstraintStorage =
        DenseStorage<DistanceConstraint, DistanceConstraintTag>;

    using AngleConstraintStorage =
        DenseStorage<AngleConstraint, AngleConstraintTag>;

    using VolumeConstraintStorage =
        DenseStorage<VolumeConstraint, VolumeConstraintTag>;

    using PointsColliderStorage =
        DenseStorage<PointsCollider, PointsColliderTag>;

    using PolygonColliderStorage =
        DenseStorage<PointsCollider, PolygonColliderTag>;

    using AngularMotorStorage =
        DenseStorage<AngularMotor, AngularMotorTag>;

    using AngularServoStorage =
        DenseStorage<AngularServo, AngularServoTag>;

    struct ParticleStorage : DenseStorage<Particle, ParticleTag>
    {
        ParticleId create(
            const glm::vec2 &pos,
            float mass,
            const glm::vec2 &vel = {})
        {
            Particle particle{pos, pos, vel, mass > 0.f ? 1.f / mass : 0.f};

            return DenseStorage::create(std::move(particle));
        }
    };

    struct EntityStorage : DenseStorage<Entity, EntityTag>
    {
        EntityId create(Entity entity)
        {
            EntityId id = DenseStorage::create(std::move(entity));

            get(id)->id = id;

            return id;
        }
    };
} // namespace xpbd
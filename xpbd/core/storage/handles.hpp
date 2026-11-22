#pragma once
#include <cstdint>
#include <limits>

namespace xpbd
{
    template <typename Tag>
    struct Handle
    {
        static constexpr uint32_t invalidIndex =
            std::numeric_limits<uint32_t>::max();

        uint32_t index = invalidIndex;
        uint32_t generation = 0;

        bool valid() const noexcept
        {
            return index != invalidIndex;
        }

        bool operator==(const Handle &other) const
        {
            return index == other.index &&
                   generation == other.generation;
        }

        bool operator!=(const Handle &other) const
        {
            return !(*this == other);
        }
    };

    struct EntityTag;
    struct ParticleTag;
    struct DistanceConstraintTag;
    struct AngleConstraintTag;
    struct VolumeConstraintTag;
    struct PointsColliderTag;
    struct PolygonColliderTag;
    struct AngularMotorTag;
    struct AngularServoTag;

    using EntityId = Handle<EntityTag>;
    using ParticleId = Handle<ParticleTag>;
    using DistanceConstraintId = Handle<DistanceConstraintTag>;
    using AngleConstraintId = Handle<AngleConstraintTag>;
    using VolumeConstraintId = Handle<VolumeConstraintTag>;
    using PointsColliderId = Handle<PointsColliderTag>;
    using PolygonColliderId = Handle<PolygonColliderTag>;
    using AngularMotorId = Handle<AngularMotorTag>;
    using AngularServoId = Handle<AngularServoTag>;
} // namespace xpbd
#pragma once
#include <cstdint>
#include <limits>

namespace xpbd
{
    template<typename Tag>
    struct Handle
    {
        static constexpr uint32_t invalidIndex =
            std::numeric_limits<uint32_t>::max();

        uint32_t index = invalidIndex;
        uint32_t generation = 0;

        bool valid() const
        {
            return index != invalidIndex;
        }

        explicit operator bool() const
        {
            return valid();
        }

        bool operator==(const Handle&) const = default;
    };

    struct EntityTag {};
    struct ParticleTag {};
    struct DistanceConstraintTag {};
    struct AngleConstraintTag {};
    struct VolumeConstraintTag {};
    struct PointsColliderTag {};
    struct PolygonColliderTag {};
    struct MotorTag {};
    struct ServoTag {};

    using EntityId = Handle<EntityTag>;
    using ParticleId = Handle<ParticleTag>;
    using DistanceConstraintId = Handle<DistanceConstraintTag>;
    using AngleConstraintId = Handle<AngleConstraintTag>;
    using VolumeConstraintId = Handle<VolumeConstraintTag>;
    using PointsColliderId = Handle<PointsColliderTag>;
    using PolygonColliderId = Handle<PolygonColliderTag>;
    using MotorId = Handle<MotorTag>;
    using ServoId = Handle<ServoTag>;
}
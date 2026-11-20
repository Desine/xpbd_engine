#pragma once
#include "handles.hpp"
#include "object_pool.hpp"
#include <string>
#include <vector>
#include <optional>

namespace xpbd
{
    struct Entity
    {
        std::string name;
        EntityId id;
        EntityId parent;
        std::vector<EntityId> children;

        std::vector<ParticleId> particles;

        std::vector<DistanceConstraintId> distanceConstraints;
        std::vector<AngleConstraintId> angleConstraints;
        std::vector<VolumeConstraintId> volumeConstraints;

        std::vector<PointsColliderId> pointsColliders;
        std::vector<PolygonColliderId> polygonColliders;

        std::vector<MotorId> motors;
        std::vector<ServoId> servos;
    };

    struct EntityStorage
    {
    private:
        std::vector<Entity> entities;
        ObjectPool<EntityTag> handles;

    public:
        EntityId create(Entity entity);
        Entity *get(EntityId id);
        bool remove(EntityId id);
        bool contains(EntityId id)const;
        size_t size()const;
    };
} // namespace xpbd
#include "entity.hpp"
#include "handles.hpp"
#include "object_pool.hpp"
#include <string>
#include <vector>
#include <optional>

namespace xpbd
{
    EntityId EntityStorage::create(Entity entity)
    {
        EntityId id =
            handles.create(static_cast<uint32_t>(entities.size()));
        entity.id = id;

        entities.push_back(std::move(entity));

        return entity.id;
    }

    Entity *EntityStorage::get(EntityId id)
    {
        auto index = handles.get_dense_index(id);
        if (!index)
            return nullptr;

        return &entities[*index];
    }

    bool EntityStorage::remove(EntityId id)
    {
        auto index = handles.get_dense_index(id);
        if (!index)
            return false;

        uint32_t denseIndex = *index;
        uint32_t lastIndex =
            static_cast<uint32_t>(entities.size() - 1);

        if (denseIndex != lastIndex)
        {
            entities[denseIndex] = std::move(entities[lastIndex]);
        }

        entities.pop_back();

        return handles.remove(id);
    }
    bool EntityStorage::contains(EntityId id) const
    {
        return handles.contains(id);
    }

    size_t EntityStorage::size() const
    {
        return entities.size();
    }

} // namespace xpbd
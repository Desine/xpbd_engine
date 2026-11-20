#pragma once
#include "handles.hpp"

#include <optional>
#include <vector>
#include <cstdint>
#include <cassert>

namespace xpbd
{
    template <typename Tag>
    class ObjectPool
    {

    private:
        struct Slot
        {
            uint32_t generation = 0;
            uint32_t denseIndex = 0;
            bool alive = false;
        };

        std::vector<Slot> slots;
        std::vector<uint32_t> freeSlots;
        std::vector<uint32_t> denseToSlot;

    public:
        using Id = Handle<Tag>;

        Id create(uint32_t denseIndex)
        {
            assert(denseIndex == denseToSlot.size());

            uint32_t index;

            if (!freeSlots.empty())
            {
                index = freeSlots.back();
                freeSlots.pop_back();
            }
            else
            {
                index = static_cast<uint32_t>(slots.size());
                slots.emplace_back();
            }

            Slot &slot = slots[index];

            slot.alive = true;
            slot.denseIndex = denseIndex;

            denseToSlot.push_back(index);

            return {
                index,
                slot.generation};
        }

        bool contains(Id id) const
        {
            if (!id.valid() || id.index >= slots.size())
                return false;

            const Slot &slot = slots[id.index];

            return slot.alive &&
                   slot.generation == id.generation;
        }

        std::optional<uint32_t> get_dense_index(Id id) const
        {
            if (contains(id))
                return slots[id.index].denseIndex;

            return std::nullopt;
        }
        std::optional<Id> get_handle(uint32_t denseIndex) const
        {
            if (denseIndex >= denseToSlot.size())
                return std::nullopt;

            uint32_t slotIndex = denseToSlot[denseIndex];

            const Slot &slot = slots[slotIndex];

            if (!slot.alive)
                return std::nullopt;

            return Id{
                slotIndex,
                slot.generation};
        }

        bool update_dense_index(Id id, uint32_t newIndex)
        {
            if (!contains(id))
                return false;

            if (newIndex >= denseToSlot.size())
                return false;

            Slot &slot = slots[id.index];

            slot.denseIndex = newIndex;
            denseToSlot[newIndex] = id.index;

            return true;
        }

        bool remove(Id id)
        {
            if (!contains(id))
                return false;

            Slot &slot = slots[id.index];

            uint32_t denseIndex = slot.denseIndex;
            uint32_t lastDenseIndex =
                static_cast<uint32_t>(denseToSlot.size() - 1);

            if (denseIndex != lastDenseIndex)
            {
                uint32_t movedSlotIndex = denseToSlot[lastDenseIndex];

                denseToSlot[denseIndex] = movedSlotIndex;

                slots[movedSlotIndex].denseIndex = denseIndex;
            }

            denseToSlot.pop_back();

            slot.alive = false;

            if (slot.generation != std::numeric_limits<uint32_t>::max())
            {
                slot.generation++;
                freeSlots.push_back(id.index);
            }

            return true;
        }

        size_t size() const
        {
            return denseToSlot.size();
        }
    };
} // namespace xpbd
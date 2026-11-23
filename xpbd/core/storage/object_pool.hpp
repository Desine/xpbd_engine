#pragma once
#include "handles.hpp"

#include <vector>
#include <cassert>
#include <optional>
#include <stdexcept>


namespace xpbd
{
template <typename Tag>
class ObjectPool
{
public:
    using Id = Handle<Tag>;

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
    Id create();
    bool contains(Id id) const noexcept;
    std::optional<uint32_t> get_dense_index(Id id) const noexcept;
    std::optional<Id> get_handle(uint32_t denseIndex) const noexcept;
    bool remove(Id id);
    size_t size() const noexcept;
    bool empty() const noexcept;
    void reserve(size_t capacity);
    void clear();
};

template <typename Tag>
ObjectPool<Tag>::Id ObjectPool<Tag>::create()
{
    if (denseToSlot.size() >= Id::invalidIndex)
        throw std::length_error("ObjectPool capacity exceeded");

    const uint32_t denseIndex =
        static_cast<uint32_t>(denseToSlot.size());

    const bool reuseSlot = !freeSlots.empty();

    uint32_t slotIndex;

    if (reuseSlot)
    {
        slotIndex = freeSlots.back();
    }
    else
    {
        if (slots.size() >= Id::invalidIndex)
            throw std::length_error("ObjectPool capacity exceeded");

        slotIndex = static_cast<uint32_t>(slots.size());
        slots.emplace_back();
    }

    try
    {
        denseToSlot.push_back(slotIndex);
    }
    catch (...)
    {
        if (!reuseSlot)
            slots.pop_back();

        throw;
    }

    if (reuseSlot)
        freeSlots.pop_back();

    Slot &slot = slots[slotIndex];

    slot.alive = true;
    slot.denseIndex = denseIndex;

    return {slotIndex, slot.generation};
}

template <typename Tag>
bool ObjectPool<Tag>::contains(Id id) const noexcept
{
    if (!id.valid() || id.index >= slots.size())
        return false;

    const Slot &slot = slots[id.index];

    return slot.alive &&
           slot.generation == id.generation;
}

template <typename Tag>
std::optional<uint32_t> ObjectPool<Tag>::get_dense_index(Id id) const noexcept
{
    if (!contains(id))
        return std::nullopt;

    return slots[id.index].denseIndex;
}

template <typename Tag>
std::optional<typename ObjectPool<Tag>::Id> ObjectPool<Tag>::get_handle(uint32_t denseIndex) const noexcept
{
    if (denseIndex >= denseToSlot.size())
        return std::nullopt;

    const uint32_t slotIndex = denseToSlot[denseIndex];

    assert(slotIndex < slots.size());

    const Slot &slot = slots[slotIndex];

    assert(slot.alive);
    assert(slot.denseIndex == denseIndex);

    return Id{
        slotIndex,
        slot.generation};
}

template <typename Tag>
bool ObjectPool<Tag>::remove(Id id)
{
    if (!contains(id))
        return false;

    Slot &slot = slots[id.index];

    const uint32_t denseIndex = slot.denseIndex;
    const uint32_t lastDenseIndex =
        static_cast<uint32_t>(denseToSlot.size() - 1);

    const bool canReuse =
        slot.generation != std::numeric_limits<uint32_t>::max();

    if (canReuse)
        freeSlots.push_back(id.index);

    if (denseIndex != lastDenseIndex)
    {
        const uint32_t movedSlotIndex =
            denseToSlot[lastDenseIndex];

        denseToSlot[denseIndex] = movedSlotIndex;

        slots[movedSlotIndex].denseIndex = denseIndex;
    }

    denseToSlot.pop_back();

    slot.alive = false;

    if (canReuse)
        ++slot.generation;

    return true;
}

template <typename Tag>
size_t ObjectPool<Tag>::size() const noexcept
{
    return denseToSlot.size();
}

template <typename Tag>
bool ObjectPool<Tag>::empty() const noexcept
{
    return denseToSlot.empty();
}

template <typename Tag>
void ObjectPool<Tag>::reserve(size_t capacity)
{
    slots.reserve(capacity);
    denseToSlot.reserve(capacity);
}

template <typename Tag>
void ObjectPool<Tag>::clear()
{
    freeSlots.reserve(slots.size());

    while (!denseToSlot.empty())
    {
        auto id = get_handle(
            static_cast<uint32_t>(denseToSlot.size() - 1));

        assert(id);
        remove(*id);
    }
}
} // namespace xpbd
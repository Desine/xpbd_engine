#pragma once
#include "handles.hpp"
#include "object_pool.hpp"

#include <vector>

namespace xpbd
{
template <typename Data, typename Tag>
class DenseStorage
{
public:
    using Id = Handle<Tag>;

protected:
    ObjectPool<Tag> handles;
    std::vector<Data> data;

public:
    static_assert(
        std::is_nothrow_move_constructible_v<Data>,
        "Data must be nothrow move constructible");

    static_assert(
        std::is_nothrow_move_assignable_v<Data>,
        "Data must be nothrow move assignable");

    Id create(Data value = {});
    Data *get(Id id) noexcept;
    const Data *get(Id id) const noexcept;
    bool remove(Id id);
    bool contains(Id id) const noexcept;
    size_t size() const noexcept;
    bool empty() const noexcept;

    std::optional<Id> handle_at(size_t denseIndex) const noexcept;
    Data *at(size_t denseIndex) noexcept;
    const Data *at(size_t denseIndex) const noexcept;
    void reserve(size_t capacity);
    void clear();
};

template <typename Data, typename Tag>
DenseStorage<Data, Tag>::Id DenseStorage<Data, Tag>::create(Data value)
{
    data.emplace_back(std::move(value));

    try
    {
        return handles.create();
    }
    catch (...)
    {
        data.pop_back();
        throw;
    }
}

template <typename Data, typename Tag>
Data *DenseStorage<Data, Tag>::get(Id id) noexcept
{
    auto index = handles.get_dense_index(id);

    if (!index)
        return nullptr;

    return &data[*index];
}

template <typename Data, typename Tag>
const Data *DenseStorage<Data, Tag>::get(Id id) const noexcept
{
    auto index = handles.get_dense_index(id);

    if (!index)
        return nullptr;

    return &data[*index];
}

template <typename Data, typename Tag>
bool DenseStorage<Data, Tag>::remove(Id id)
{
    auto index = handles.get_dense_index(id);

    if (!index)
        return false;

    const size_t denseIndex = *index;
    const size_t lastIndex = data.size() - 1;

    if (!handles.remove(id))
        return false;

    if (denseIndex != lastIndex)
        data[denseIndex] = std::move(data[lastIndex]);

    data.pop_back();

    return true;
}

template <typename Data, typename Tag>
bool DenseStorage<Data, Tag>::contains(Id id) const noexcept
{
    return handles.contains(id);
}

template <typename Data, typename Tag>
size_t DenseStorage<Data, Tag>::size() const noexcept
{
    return data.size();
}

template <typename Data, typename Tag>
bool DenseStorage<Data, Tag>::empty() const noexcept
{
    return data.empty();
}

template <typename Data, typename Tag>
std::optional<typename DenseStorage<Data, Tag>::Id> DenseStorage<Data, Tag>::handle_at(size_t denseIndex) const noexcept
{
    if (denseIndex >= data.size())
        return std::nullopt;

    return handles.get_handle(
        static_cast<uint32_t>(denseIndex));
}

template <typename Data, typename Tag>
Data *DenseStorage<Data, Tag>::at(size_t denseIndex) noexcept
{
    if (denseIndex >= data.size())
        return nullptr;

    return &data[denseIndex];
}

template <typename Data, typename Tag>
const Data *DenseStorage<Data, Tag>::at(size_t denseIndex) const noexcept
{
    if (denseIndex >= data.size())
        return nullptr;

    return &data[denseIndex];
}

template <typename Data, typename Tag>
void DenseStorage<Data, Tag>::reserve(size_t capacity)
{
    data.reserve(capacity);
    handles.reserve(capacity);
}

template <typename Data, typename Tag>
void DenseStorage<Data, Tag>::clear()
{
    handles.clear();
    data.clear();
}

} // namespace xpbd
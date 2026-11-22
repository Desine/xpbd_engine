#pragma once
#include "handles.hpp"
#include "object_pool.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <type_traits>
#include <utility>
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

        Id create(Data value)
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

        Data *get(Id id) noexcept
        {
            auto index = handles.get_dense_index(id);

            if (!index)
                return nullptr;

            return &data[*index];
        }

        const Data *get(Id id) const noexcept
        {
            auto index = handles.get_dense_index(id);

            if (!index)
                return nullptr;

            return &data[*index];
        }

        bool remove(Id id)
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

        bool contains(Id id) const noexcept
        {
            return handles.contains(id);
        }

        size_t size() const noexcept
        {
            return data.size();
        }

        bool empty() const noexcept
        {
            return data.empty();
        }

        std::optional<Id> handle_at(size_t denseIndex) const noexcept
        {
            if (denseIndex >= data.size())
                return std::nullopt;

            return handles.get_handle(
                static_cast<uint32_t>(denseIndex));
        }

        Data *at(size_t denseIndex) noexcept
        {
            if (denseIndex >= data.size())
                return nullptr;

            return &data[denseIndex];
        }

        const Data *at(size_t denseIndex) const noexcept
        {
            if (denseIndex >= data.size())
                return nullptr;

            return &data[denseIndex];
        }

        void reserve(size_t capacity)
        {
            data.reserve(capacity);
            handles.reserve(capacity);
        }

        void clear()
        {
            handles.clear();
            data.clear();
        }
    };
}
#pragma once
#include "xpbd/core/types.hpp"

#include "glm/glm.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <ankerl/unordered_dense.h>

namespace xpbd
{
    struct CellCoord
    {
        int x;
        int y;

        bool operator==(const CellCoord &other) const noexcept
        {
            return x == other.x && y == other.y;
        }
    };

    struct CellCoordHash
    {
        using is_avalanching = void;
        size_t operator()(const CellCoord &c) const noexcept
        {
            return ankerl::unordered_dense::hash<uint64_t>{}(
                (static_cast<uint64_t>(static_cast<uint32_t>(c.x)) << 32) |
                static_cast<uint32_t>(c.y));
        }
    };

    class SpatialHashAABB
    {
    public:
        size_t cellSize;
        size_t tableReserve;
        size_t bucketReserve;
        size_t outIdsReserve;
        ankerl::unordered_dense::map<CellCoord, std::vector<size_t>, CellCoordHash> idTable;
        ankerl::unordered_dense::map<size_t, AABB> aabbTable;

        explicit SpatialHashAABB(size_t cellSize = 5, size_t tableReserve = 1024, size_t bucketReserve = 32, size_t outIdsReserve = 10)
            : cellSize(cellSize), tableReserve(tableReserve), bucketReserve(bucketReserve), outIdsReserve(outIdsReserve)
        {
            clear();
        }

        void clear()
        {
            idTable.clear();
            idTable.reserve(tableReserve);

            aabbTable.clear();
            aabbTable.reserve(tableReserve);
        }

        void add_aabb(const AABB &aabb, size_t id)
        {
            int startX, endX, startY, endY;
            get_cell_range(aabb, startX, endX, startY, endY);

            for (int cy = startY; cy <= endY; ++cy)
            {
                for (int cx = startX; cx <= endX; ++cx)
                {
                    auto &bucket = idTable[CellCoord{cx, cy}];
                    bucket.reserve(bucketReserve);
                    bucket.emplace_back(id);
                }
            }
            aabbTable.insert({id, aabb});
        }

        void get_overlapping_aabb_ids(std::vector<size_t> &out, const AABB &aabb) const
        {
            out.clear();
            out.reserve(outIdsReserve);
            int startX, endX, startY, endY;
            get_cell_range(aabb, startX, endX, startY, endY);

            for (int y = startY; y <= endY; ++y)
            {
                for (int x = startX; x <= endX; ++x)
                {
                    auto idT = idTable.find(CellCoord{x, y});
                    for (size_t id : idT->second)
                    {
                        auto aabbT = aabbTable.find(id);
                        if (aabbT != aabbTable.end() && layer_mask_collide(aabb, aabbT->second) && aabb.intersects(aabbT->second))
                        {
                            out.emplace_back(id);
                        }
                    }
                }
            }
            dedup_sort(out);
        }

        void get_overlapping_aabb_ids_excludeId(std::vector<size_t> &out, const AABB &aabb, size_t excludeId) const
        {
            out.clear();
            out.reserve(outIdsReserve);
            int startX, endX, startY, endY;
            get_cell_range(aabb, startX, endX, startY, endY);

            for (int y = startY; y <= endY; ++y)
            {
                for (int x = startX; x <= endX; ++x)
                {
                    auto idT = idTable.find(CellCoord{x, y});
                    if (idT != idTable.end())
                    {
                        for (size_t id : idT->second)
                        {
                            if (id == excludeId)
                                continue;
                            auto aabbT = aabbTable.find(id);
                            if (aabbT != aabbTable.end() && layer_mask_collide(aabb, aabbT->second) && aabb.intersects(aabbT->second))
                            {
                                out.emplace_back(id);
                            }
                        }
                    }
                }
            }
            dedup_sort(out);
        }

    private:
        static bool layer_mask_collide(AABB points, AABB polygon)
        {
            return points.mask & polygon.layer || points.layer & polygon.mask;
        }

        inline void get_cell_range(const AABB &aabb,
                                   int &startX, int &endX,
                                   int &startY, int &endY) const noexcept
        {
            startX = static_cast<int>(std::floor(aabb.l / cellSize));
            endX = static_cast<int>(std::floor(aabb.r / cellSize));
            startY = static_cast<int>(std::floor(aabb.b / cellSize));
            endY = static_cast<int>(std::floor(aabb.t / cellSize));
        }

        static inline void dedup_sort(std::vector<size_t> &vec)
        {
            std::sort(vec.begin(), vec.end());
            vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
        }
    };
} // namespace xpbd
#ifndef define_SpatialHashAABB
#define define_SpatialHashAABB
#include <vector>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <ankerl/unordered_dense.h>

#include <stdio.h>
#include <string.h>

struct AABB
{
    float l, r, b, t;

    bool intersects(const AABB &o) const noexcept
    {
        return !(r < o.l || l > o.r || t < o.b || b > o.t);
    }

    AABB get_intersection(const AABB &o) const noexcept
    {
        return {
            std::max(l, o.l),
            std::min(r, o.r),
            std::max(b, o.b),
            std::min(t, o.t)};
    }
};

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
    ankerl::unordered_dense::map<CellCoord, std::vector<size_t>, CellCoordHash> table;

    explicit SpatialHashAABB(size_t cellSize = 300, size_t tableReserve = 1024)
        : cellSize(cellSize), tableReserve(tableReserve)
    {
        clear();
    }

    void clear()
    {
        table.clear();
        table.reserve(tableReserve);
    }

    void add_aabb(const AABB &aabb, size_t id)
    {
        int startX, endX, startY, endY;
        get_cell_range(aabb, startX, endX, startY, endY);

        for (int cy = startY; cy <= endY; ++cy)
        {
            for (int cx = startX; cx <= endX; ++cx)
            {
                auto &bucket = table[CellCoord{cx, cy}];
                bucket.push_back(id);
            }
        }
    }

    std::vector<size_t> get_overlapping_aabb_ids(const AABB &aabb) const
    {
        std::vector<size_t> out;
        int startX, endX, startY, endY;
        get_cell_range(aabb, startX, endX, startY, endY);

        for (int cy = startY; cy <= endY; ++cy)
        {
            for (int cx = startX; cx <= endX; ++cx)
            {
                auto it = table.find(CellCoord{cx, cy});
                if (it != table.end())
                {
                    out.insert(out.end(), it->second.begin(), it->second.end());
                }
            }
        }
        dedup_sort(out);
        return out;
    }

    std::vector<size_t> get_overlapping_aabb_ids_excludeId(const AABB &aabb, size_t excludeId) const
    {
        std::vector<size_t> out;
        int startX, endX, startY, endY;
        get_cell_range(aabb, startX, endX, startY, endY);

        for (int cy = startY; cy <= endY; ++cy)
        {
            for (int cx = startX; cx <= endX; ++cx)
            {
                auto it = table.find(CellCoord{cx, cy});
                if (it != table.end())
                {
                    for (size_t id : it->second)
                    {
                        if (id != excludeId)
                        {
                            out.push_back(id);
                        }
                    }
                }
            }
        }
        dedup_sort(out);
        return out;
    }

private:
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

#endif // define_SpatialHashAABB

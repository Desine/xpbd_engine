#ifndef define_SpatialHashAABB
#define define_SpatialHashAABB
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstddef>
#include <utility>
#include <algorithm>

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

struct PairHash
{
    size_t operator()(const std::pair<int, int> &p) const noexcept
    {
        uint64_t x = static_cast<uint32_t>(p.first);
        uint64_t y = static_cast<uint32_t>(p.second);
        return (x * 73856093ull) ^ (y * 19349663ull);
    }
};

class SpatialHashAABB
{
public:
    size_t cellSize;
    size_t tableReserve;
    std::unordered_map<std::pair<int, int>, std::vector<size_t>, PairHash> table;

    explicit SpatialHashAABB(size_t cellSize = 500, size_t tableReserve = 1024) : cellSize(cellSize), tableReserve(tableReserve)
    {
        clear();
    }

    void clear()
    {
        table.clear();
        
        if (table.max_size() != tableReserve)
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
                auto &bucket = table.try_emplace({cx, cy}).first->second;
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
                auto it = table.find({cx, cy});
                if (it != table.end())
                {
                    for (size_t id : it->second)
                    {
                        out.push_back(id);
                    }
                }
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    std::vector<size_t> get_overlapping_aabb_ids_excludeId(const AABB &aabb, int excludeId) const
    {
        std::vector<size_t> out;
        
        int startX, endX, startY, endY;
        get_cell_range(aabb, startX, endX, startY, endY);

        for (int cy = startY; cy <= endY; ++cy)
        {
            for (int cx = startX; cx <= endX; ++cx)
            {
                auto it = table.find({cx, cy});
                if (it != table.end())
                {
                    for (size_t id : it->second)
                    {
                        if (id != excludeId)
                            out.push_back(id);
                    }
                }
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

private:
    inline void get_cell_range(const AABB &aabb, int &startX, int &endX, int &startY, int &endY) const noexcept
    {
        startX = static_cast<int>(std::floor(aabb.l / cellSize));
        endX = static_cast<int>(std::floor(aabb.r / cellSize));
        startY = static_cast<int>(std::floor(aabb.b / cellSize));
        endY = static_cast<int>(std::floor(aabb.t / cellSize));
    }
};

#endif // define_SpatialHashAABB

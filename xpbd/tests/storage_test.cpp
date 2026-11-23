#include "xpbd/core/storage/object_pool.hpp"
#include "xpbd/core/storage/dense_storage.hpp"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <type_traits>

namespace xpbd
{

int testsPassed = 0;

#define CHECK(expr)                                            \
    do                                                         \
    {                                                          \
        if (!(expr))                                           \
        {                                                      \
            std::cerr << "FAILED: " << #expr                   \
                      << " at " << __FILE__ << ":" << __LINE__ \
                      << '\n';                                 \
            std::exit(EXIT_FAILURE);                           \
        }                                                      \
    } while (false)

struct TestTagA
{
};
struct TestTagB
{
};

using TestIdA = Handle<TestTagA>;
using TestIdB = Handle<TestTagB>;

struct TestData
{
    int value = 0;

    TestData() = default;
    explicit TestData(int v) : value(v)
    {
    }

    TestData(const TestData &) = default;
    TestData &operator=(const TestData &) = default;

    TestData(TestData &&other) noexcept
        : value(other.value)
    {
    }

    TestData &operator=(TestData &&other) noexcept
    {
        value = other.value;
        return *this;
    }
};

void test_object_pool_create_and_contains()
{
    ObjectPool<TestTagA> pool;

    const auto a = pool.create();
    const auto b = pool.create();
    const auto c = pool.create();

    CHECK(pool.size() == 3);
    CHECK(pool.contains(a));
    CHECK(pool.contains(b));
    CHECK(pool.contains(c));

    CHECK(pool.get_dense_index(a).value() == 0);
    CHECK(pool.get_dense_index(b).value() == 1);
    CHECK(pool.get_dense_index(c).value() == 2);

    ++testsPassed;
}

void test_object_pool_swap_and_pop()
{
    ObjectPool<TestTagA> pool;

    const auto a = pool.create();
    const auto b = pool.create();
    const auto c = pool.create();

    CHECK(pool.remove(b));

    CHECK(pool.size() == 2);
    CHECK(pool.contains(a));
    CHECK(!pool.contains(b));
    CHECK(pool.contains(c));

    CHECK(pool.get_dense_index(a).value() == 0);
    CHECK(pool.get_dense_index(c).value() == 1);

    CHECK(pool.get_handle(0).value() == a);
    CHECK(pool.get_handle(1).value() == c);

    ++testsPassed;
}

void test_object_pool_stale_handle()
{
    ObjectPool<TestTagA> pool;

    const auto oldId = pool.create();

    CHECK(pool.remove(oldId));
    CHECK(!pool.contains(oldId));
    CHECK(!pool.remove(oldId));

    const auto newId = pool.create();

    CHECK(pool.contains(newId));
    CHECK(newId.index == oldId.index);
    CHECK(newId.generation != oldId.generation);
    CHECK(!pool.contains(oldId));

    ++testsPassed;
}

void test_object_pool_clear()
{
    ObjectPool<TestTagA> pool;

    const auto a = pool.create();
    const auto b = pool.create();
    const auto c = pool.create();

    pool.clear();

    CHECK(pool.empty());
    CHECK(!pool.contains(a));
    CHECK(!pool.contains(b));
    CHECK(!pool.contains(c));

    const auto d = pool.create();

    CHECK(pool.contains(d));
    CHECK(!pool.contains(a));
    CHECK(!pool.contains(b));
    CHECK(!pool.contains(c));

    ++testsPassed;
}

void test_object_pool_reserve()
{
    ObjectPool<TestTagA> pool;

    const auto a = pool.create();
    const auto b = pool.create();

    pool.reserve(1024);

    CHECK(pool.contains(a));
    CHECK(pool.contains(b));
    CHECK(pool.get_dense_index(a).value() == 0);
    CHECK(pool.get_dense_index(b).value() == 1);

    ++testsPassed;
}

void test_different_handle_types()
{
    static_assert(!std::is_same_v<TestIdA, TestIdB>);

    ObjectPool<TestTagA> poolA;
    ObjectPool<TestTagB> poolB;

    const auto a = poolA.create();
    const auto b = poolB.create();

    CHECK(poolA.contains(a));
    CHECK(poolB.contains(b));

    ++testsPassed;
}

void test_dense_storage_create_and_get()
{
    DenseStorage<TestData, TestTagA> storage;

    const auto a = storage.create(TestData{10});
    const auto b = storage.create(TestData{20});

    CHECK(storage.size() == 2);
    CHECK(storage.contains(a));
    CHECK(storage.contains(b));

    CHECK(storage.get(a) != nullptr);
    CHECK(storage.get(b) != nullptr);

    CHECK(storage.get(a)->value == 10);
    CHECK(storage.get(b)->value == 20);

    ++testsPassed;
}

void test_dense_storage_remove_middle()
{
    DenseStorage<TestData, TestTagA> storage;

    const auto a = storage.create(TestData{10});
    const auto b = storage.create(TestData{20});
    const auto c = storage.create(TestData{30});

    CHECK(storage.remove(b));

    CHECK(storage.size() == 2);
    CHECK(storage.get(a)->value == 10);
    CHECK(storage.get(c)->value == 30);

    CHECK(storage.get(b) == nullptr);
    CHECK(!storage.remove(b));

    CHECK(storage.handle_at(0).value() == a);
    CHECK(storage.handle_at(1).value() == c);

    ++testsPassed;
}

void test_dense_storage_reuse_slot()
{
    DenseStorage<TestData, TestTagA> storage;

    const auto oldId = storage.create(TestData{10});
    CHECK(storage.remove(oldId));

    const auto newId = storage.create(TestData{20});

    CHECK(newId.index == oldId.index);
    CHECK(newId.generation != oldId.generation);

    CHECK(storage.get(oldId) == nullptr);
    CHECK(storage.get(newId) != nullptr);
    CHECK(storage.get(newId)->value == 20);

    ++testsPassed;
}

void test_dense_storage_clear()
{
    DenseStorage<TestData, TestTagA> storage;

    const auto a = storage.create(TestData{10});
    const auto b = storage.create(TestData{20});

    storage.clear();

    CHECK(storage.empty());
    CHECK(storage.get(a) == nullptr);
    CHECK(storage.get(b) == nullptr);

    const auto c = storage.create(TestData{30});

    CHECK(storage.get(c) != nullptr);
    CHECK(storage.get(c)->value == 30);
    CHECK(storage.get(a) == nullptr);
    CHECK(storage.get(b) == nullptr);

    ++testsPassed;
}

} // namespace xpbd

int main()
{
    using namespace xpbd;
    test_object_pool_create_and_contains();
    test_object_pool_swap_and_pop();
    test_object_pool_stale_handle();
    test_object_pool_clear();
    test_object_pool_reserve();
    test_different_handle_types();

    test_dense_storage_create_and_get();
    test_dense_storage_remove_middle();
    test_dense_storage_reuse_slot();
    test_dense_storage_clear();

    std::cout << "All " << testsPassed << " storage tests passed.\n";
    return EXIT_SUCCESS;
}
#include <benchmark/benchmark.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "CHTree.h"
#include <random>
#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

// ============================================================================
// Insertion Benchmarks
// ============================================================================

// Benchmark inserting N points into an empty hull
static void BM_Insert(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    
    for (auto _ : state) {
        CHTree<K> tree;
        for (const auto& p : points) {
            tree.insert(p);
        }
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_Insert)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// Benchmark single insertion into a hull of size N
static void BM_InsertSingle(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    
    // Pre-build a hull of size n
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    
    CHTree<K> tree;
    for (const auto& p : points) {
        tree.insert(p);
    }
    
    // Generate new points to insert
    std::vector<Point_2> new_points;
    for (int i = 0; i < 1000; i++) {
        new_points.emplace_back(dist(rng), dist(rng));
    }
    
    int idx = 0;
    for (auto _ : state) {
        tree.insert(new_points[idx % new_points.size()]);
        idx++;
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_InsertSingle)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// ============================================================================
// Removal Benchmarks
// ============================================================================

// Benchmark removing points from a hull of size N
static void BM_Remove(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    
    for (auto _ : state) {
        state.PauseTiming();
        CHTree<K> tree;
        for (const auto& p : points) {
            tree.insert(p);
        }
        std::vector<Point_2> to_remove = points;
        std::shuffle(to_remove.begin(), to_remove.end(), rng);
        state.ResumeTiming();
        
        for (const auto& p : to_remove) {
            tree.remove(p);
        }
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_Remove)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// Benchmark single removal from a hull of size N
static void BM_RemoveSingle(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    
    std::vector<Point_2> points;
    for (int i = 0; i < n * 2; i++) {  // Extra points to remove and re-add
        points.emplace_back(dist(rng), dist(rng));
    }
    
    CHTree<K> tree;
    for (int i = 0; i < n; i++) {
        tree.insert(points[i]);
    }
    
    int idx = n;
    for (auto _ : state) {
        // Remove a random existing point
        size_t remove_idx = rng() % n;
        tree.remove(points[remove_idx]);
        // Re-add a new point to maintain size
        tree.insert(points[idx % points.size()]);
        idx++;
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_RemoveSingle)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// ============================================================================
// Query Benchmarks
// ============================================================================

// Benchmark covers() query on a hull of size N
static void BM_Covers(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    
    CHTree<K> tree;
    for (const auto& p : points) {
        tree.insert(p);
    }
    
    // Generate query points
    std::vector<Point_2> queries;
    for (int i = 0; i < 1000; i++) {
        queries.emplace_back(dist(rng), dist(rng));
    }
    
    int idx = 0;
    for (auto _ : state) {
        bool result = tree.covers(queries[idx % queries.size()]);
        benchmark::DoNotOptimize(result);
        idx++;
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_Covers)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// Benchmark upperHullPoints() retrieval
static void BM_UpperHullPoints(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    
    CHTree<K> tree;
    for (const auto& p : points) {
        tree.insert(p);
    }
    
    for (auto _ : state) {
        auto hull = tree.upperHullPoints();
        benchmark::DoNotOptimize(hull);
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_UpperHullPoints)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// ============================================================================
// Mixed Operations Benchmark
// ============================================================================

static void BM_MixedOperations(benchmark::State& state) {
    const int n = state.range(0);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    std::uniform_int_distribution<int> action(0, 4);  // 0=remove, 1-4=insert
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    
    CHTree<K> tree;
    for (const auto& p : points) {
        tree.insert(p);
    }
    
    for (auto _ : state) {
        if (!points.empty() && action(rng) == 0) {
            // Remove
            size_t idx = rng() % points.size();
            tree.remove(points[idx]);
            points.erase(points.begin() + idx);
        } else {
            // Insert
            Point_2 p(dist(rng), dist(rng));
            tree.insert(p);
            points.push_back(p);
        }
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_MixedOperations)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// ============================================================================
// O(n) Build Construction Benchmarks
// ============================================================================

// Helper to generate sorted points
static std::vector<Point_2> generateSortedPoints(int n, uint64_t seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    std::sort(points.begin(), points.end(), [](const Point_2& a, const Point_2& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    return points;
}

// Benchmark O(n) build from sorted points
static void BM_Build(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    
    for (auto _ : state) {
        CHTree<K> tree;
        tree.build(points);
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_Build)->RangeMultiplier(2)->Range(64, 16384)->Complexity();

// Benchmark incremental insertion for comparison
static void BM_IncrementalInsert(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    
    for (auto _ : state) {
        CHTree<K> tree;
        for (const auto& p : points) {
            tree.insert(p);
        }
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_IncrementalInsert)->RangeMultiplier(2)->Range(64, 16384)->Complexity();

// Benchmark build vs insert speedup (larger sizes to see the difference)
static void BM_BuildLarge(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    
    for (auto _ : state) {
        CHTree<K> tree;
        tree.build(points);
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_BuildLarge)->RangeMultiplier(2)->Range(1024, 65536)->Complexity();

static void BM_IncrementalInsertLarge(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    
    for (auto _ : state) {
        CHTree<K> tree;
        for (const auto& p : points) {
            tree.insert(p);
        }
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_IncrementalInsertLarge)->RangeMultiplier(2)->Range(1024, 65536)->Complexity();

// Benchmark operations after build (to verify tree validity)
static void BM_InsertAfterBuild(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    
    CHTree<K> tree;
    tree.build(points);
    
    std::mt19937 rng(123);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    std::vector<Point_2> new_points;
    for (int i = 0; i < 1000; i++) {
        new_points.emplace_back(dist(rng), dist(rng));
    }
    
    int idx = 0;
    for (auto _ : state) {
        tree.insert(new_points[idx % new_points.size()]);
        idx++;
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_InsertAfterBuild)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// Benchmark covers() after build
static void BM_CoversAfterBuild(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    
    CHTree<K> tree;
    tree.build(points);
    
    std::mt19937 rng(123);
    std::uniform_real_distribution<double> dist(-1000, 1000);
    std::vector<Point_2> queries;
    for (int i = 0; i < 1000; i++) {
        queries.emplace_back(dist(rng), dist(rng));
    }
    
    int idx = 0;
    for (auto _ : state) {
        bool result = tree.covers(queries[idx % queries.size()]);
        benchmark::DoNotOptimize(result);
        idx++;
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_CoversAfterBuild)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

// Benchmark rebuild (to test clear + build)
static void BM_Rebuild(benchmark::State& state) {
    const int n = state.range(0);
    auto points1 = generateSortedPoints(n, 42);
    auto points2 = generateSortedPoints(n, 123);
    
    CHTree<K> tree;
    tree.build(points1);
    
    bool use_first = true;
    for (auto _ : state) {
        tree.build(use_first ? points2 : points1);
        use_first = !use_first;
        benchmark::ClobberMemory();
    }
    state.SetComplexityN(n);
}
BENCHMARK(BM_Rebuild)->RangeMultiplier(2)->Range(64, 8192)->Complexity();

BENCHMARK_MAIN();

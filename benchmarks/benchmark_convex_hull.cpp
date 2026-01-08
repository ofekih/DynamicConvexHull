#include <benchmark/benchmark.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "CHTree.h"
#include <random>
#include <vector>

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

BENCHMARK_MAIN();

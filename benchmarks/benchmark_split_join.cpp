#include <benchmark/benchmark.h>
#include "inexact.h"
#include "CHTree.h"
#include <random>
#include <vector>
#include <algorithm>

using namespace dch;

using K = Inexact_kernel<double>;
using Point_2 = K::Point_2;

// Helper to generate sorted points (ensures consistent setup)
static std::vector<Point_2> generateSortedPoints(int n, uint64_t seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-100000, 100000); 
    std::vector<Point_2> points;
    points.reserve(n);
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    std::sort(points.begin(), points.end(), [](const Point_2& a, const Point_2& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    return points;
}

// ============================================================================
// Ping-Pong Benchmark (Split then Join)
// ============================================================================

static void BM_SplitJoin_PingPong(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateSortedPoints(n, 42);
    double splitX = points[n/2].x(); // Median split
    
    // Setup: Build the tree ONCE
    CHTree<K> tree;
    tree.Build(points);
    
    for (auto _ : state) {
        // Op 1: Split
        // Returns a new tree with the right half. 'tree' keeps left half.
        CHTree<K> rightHull = tree.Split(splitX);
        
        // Op 2: Join
        // Merges 'rightHull' back into 'tree'.
        // 'rightHull' becomes empty.
        tree.Join(rightHull);
        
        // Destruction of 'rightHull' happens here, but it's empty, so O(1).
    }
    
    state.SetComplexityN(n);
}

// Expected Complexity: O(log^2 N)
BENCHMARK(BM_SplitJoin_PingPong)->RangeMultiplier(2)->Range(1024, 262144)->Complexity();

BENCHMARK_MAIN();

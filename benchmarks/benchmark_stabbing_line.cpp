// Benchmark for StabbingLineStructure - verify O(log² n) complexity
#include <benchmark/benchmark.h>
#include <vector>
#include <random>
#include "StabbingLineStructure.h"
#include "inexact.h"

using namespace dch;

using K = Inexact_kernel<double>;
using Point_2 = K::Point_2;
typedef StabbingLineStructure<K> SLS;

// Generate points near y = x with small noise
std::vector<Point_2> generateNearLinePoints(int n, int seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> noise(-0.1, 0.1);
    
    std::vector<Point_2> points;
    points.reserve(n);
    for (int i = 0; i < n; ++i) {
        double x = i * 0.1;
        double y = x + noise(rng);
        points.emplace_back(x, y);
    }
    
    // Already sorted by x
    return points;
}

// Benchmark: findStabbingLine query
static void BM_FindStabbingLine(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    SLS sls(0.5);
    sls.Build(points);
    
    for (auto _ : state) {
        auto line = sls.FindStabbingLine();
        benchmark::DoNotOptimize(line);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_FindStabbingLine)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

// Benchmark: hasStabbingLine query
static void BM_HasStabbingLine(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    SLS sls(0.5);
    sls.Build(points);
    
    for (auto _ : state) {
        bool has = sls.HasStabbingLine();
        benchmark::DoNotOptimize(has);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_HasStabbingLine)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

// Benchmark: build operation
static void BM_Build(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    for (auto _ : state) {
        SLS sls(0.5);
        sls.Build(points);
        benchmark::DoNotOptimize(sls);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_Build)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

// Benchmark: amortized insert (incremental build)
static void BM_InsertAmortized(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    for (auto _ : state) {
        SLS sls(0.5);
        for (const auto& p : points) {
            sls.Insert(p);
        }
        benchmark::DoNotOptimize(sls);
    }
    
    state.SetComplexityN(n);
    // Metrics: Time per insert = Time / N
    state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_InsertAmortized)
    ->RangeMultiplier(2)
    ->Range(64, 4096)
    ->Complexity();

// Benchmark: getMinimumGap
static void BM_GetMinimumGap(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    SLS sls(0.5);
    sls.Build(points);
    
    for (auto _ : state) {
        double gap = sls.GetMinimumGap();
        benchmark::DoNotOptimize(gap);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_GetMinimumGap)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

// Benchmark: split operation
static void BM_Split(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    for (auto _ : state) {
        state.PauseTiming();
        SLS sls(0.5);
        sls.Build(points);
        state.ResumeTiming();
        
        auto right = sls.Split(n * 0.05 / 2);  // Split at midpoint
        benchmark::DoNotOptimize(right);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_Split)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

// Benchmark: join operation
static void BM_Join(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    // Build and split once to get two halves
    SLS original(0.5);
    original.Build(points);
    
    for (auto _ : state) {
        state.PauseTiming();
        SLS left(0.5);
        left.Build(std::vector<Point_2>(points.begin(), points.begin() + n/2));
        SLS right(0.5);
        right.Build(std::vector<Point_2>(points.begin() + n/2, points.end()));
        state.ResumeTiming();
        
        left.Join(right);
        benchmark::DoNotOptimize(left);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_Join)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

BENCHMARK_MAIN();

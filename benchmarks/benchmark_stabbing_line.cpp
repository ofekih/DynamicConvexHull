// Benchmark for StabbingLineStructure - verify O(log² n) complexity
#include <benchmark/benchmark.h>
#include <vector>
#include <random>
#include "StabbingLineStructure.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
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
    sls.build(points);
    
    for (auto _ : state) {
        auto line = sls.findStabbingLine();
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
    sls.build(points);
    
    for (auto _ : state) {
        bool has = sls.hasStabbingLine();
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
        sls.build(points);
        benchmark::DoNotOptimize(sls);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_Build)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

// Benchmark: single insert
static void BM_Insert(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    for (auto _ : state) {
        state.PauseTiming();
        SLS sls(0.5);
        for (int i = 0; i < n - 1; ++i) {
            sls.insert(points[i]);
        }
        state.ResumeTiming();
        
        // Time just the last insert
        sls.insert(points.back());
        benchmark::DoNotOptimize(sls);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_Insert)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 14)
    ->Complexity();

// Benchmark: getMinimumGap
static void BM_GetMinimumGap(benchmark::State& state) {
    const int n = state.range(0);
    auto points = generateNearLinePoints(n);
    
    SLS sls(0.5);
    sls.build(points);
    
    for (auto _ : state) {
        double gap = sls.getMinimumGap();
        benchmark::DoNotOptimize(gap);
    }
    
    state.SetComplexityN(n);
}
BENCHMARK(BM_GetMinimumGap)
    ->RangeMultiplier(2)
    ->Range(64, 1 << 18)
    ->Complexity();

BENCHMARK_MAIN();

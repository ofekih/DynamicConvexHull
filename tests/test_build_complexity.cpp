// Test to analyze build complexity
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <chrono>
#include "CHTree.h"
#include "inexact.h"
#include <iomanip>

using namespace dch;

using K = Inexact_kernel<double>;
using Point_2 = K::Point_2;

// Generate sorted random points
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

// Measure time for different sizes and compute complexity
TEST(BuildComplexity, AnalyzeScaling) {
    std::vector<int> sizes = {100, 200, 400, 800, 1600, 3200, 6400, 12800};
    
    std::cout << "\n=== Build Complexity Analysis ===" << std::endl;
    std::cout << "Size\t\tTime (μs)\tTime/N\t\tTime/(N*logN)" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    
    for (int n : sizes) {
        auto points = generateSortedPoints(n, 42);
        
        // Warm up
        {
            CHTree<K> tree;
            tree.Build(points);
        }
        
        // Measure
        const int iterations = 10;
        double total_time = 0;
        
        for (int i = 0; i < iterations; i++) {
            auto start = std::chrono::high_resolution_clock::now();
            CHTree<K> tree;
            tree.Build(points);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }
        
        double avg_time = total_time / iterations;
        double log_n = std::log2(n);
        double time_per_n = avg_time / n;
        double time_per_nlogn = avg_time / (n * log_n);
        
        std::cout << n << "\t\t" 
                  << std::fixed << std::setprecision(1) << avg_time << "\t\t"
                  << std::setprecision(4) << time_per_n << "\t\t"
                  << std::setprecision(4) << time_per_nlogn << std::endl;
    }
    
    std::cout << "\nIf complexity is O(N), Time/N should be roughly constant." << std::endl;
    std::cout << "If complexity is O(N log N), Time/(N*logN) should be roughly constant." << std::endl;
}

// Compare build time vs incremental insert time
TEST(BuildComplexity, CompareWithInsert) {
    std::vector<int> sizes = {500, 1000, 2000, 4000};
    
    std::cout << "\n=== Build vs Insert Comparison ===" << std::endl;
    std::cout << "Size\tBuild (μs)\tInsert (μs)\tSpeedup" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    
    for (int n : sizes) {
        auto points = generateSortedPoints(n, 42);
        
        // Measure build
        auto start = std::chrono::high_resolution_clock::now();
        CHTree<K> build_tree;
        build_tree.Build(points);
        auto end = std::chrono::high_resolution_clock::now();
        double build_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        // Measure insert
        start = std::chrono::high_resolution_clock::now();
        CHTree<K> insert_tree;
        for (const auto& p : points) {
            insert_tree.Insert(p);
        }
        end = std::chrono::high_resolution_clock::now();
        double insert_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        std::cout << n << "\t" 
                  << std::fixed << std::setprecision(0) << build_time << "\t\t"
                  << insert_time << "\t\t"
                  << std::setprecision(1) << (insert_time / build_time) << "x" << std::endl;
    }
}

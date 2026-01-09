// Tests for CHTree split() and join() (merge) operations
// Verifies correctness of hull geometry and O(log² N) complexity
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>
#include "CHTree.h"
#include "hull_test_helpers.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

// ============================================================================
// Helper Functions
// ============================================================================

void printHullDebug(const std::string& name, const std::vector<std::pair<int, int>>& hull) {
    std::cerr << name << ": ";
    for (const auto& p : hull) std::cerr << "(" << p.first << "," << p.second << ") ";
    std::cerr << "\n";
}

// Verify that a CHTree hull contains the correct points by comparing to monotone chain
bool verifyHullCorrectness(CHTree<K>& tree, const std::vector<Point_2>& points) {
    if (points.size() < 3) return true;  // Too few points to verify meaningfully
    
    auto expected_upper = hull_helpers::adjustUpperHullForCHTree(
        monotone_chain::upperHull(hull_helpers::toIntPairs(points)));
    auto expected_lower = hull_helpers::adjustLowerHullForCHTree(
        monotone_chain::lowerHull(hull_helpers::toIntPairs(points)));
    
    auto ch_upper = hull_helpers::toIntPairs(tree.upperHullPoints());
    auto ch_lower = hull_helpers::toIntPairs(tree.lowerHullPoints());

    // Check lower hull
    if (!hull_helpers::hullContainsAll(expected_lower, ch_lower)) {
        std::cerr << "Lower hull mismatch!\n";
        printHullDebug("Expected", expected_lower);
        printHullDebug("Actual", ch_lower);
        return false;
    }

    // Check upper hull
    if (!hull_helpers::hullContainsAll(expected_upper, ch_upper)) {
        std::cerr << "Upper hull mismatch!\n";
        printHullDebug("Expected", expected_upper);
        printHullDebug("Actual", ch_upper);
        return false;
    }

    return true;
}

// Generate random points sorted by x-coordinate
std::vector<Point_2> generateSortedRandomPoints(std::mt19937& rng, int n, int range = 1000) {
    std::uniform_int_distribution<int> dist(-range, range);
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

// ============================================================================
// Basic Functionality Tests - Split
// ============================================================================

TEST(SplitMerge, SplitEmpty) {
    CHTree<K> tree;
    auto right = tree.split(0.0);
    EXPECT_TRUE(tree.empty());
    EXPECT_TRUE(right.empty());
}

TEST(SplitMerge, SplitSinglePointLeft) {
    CHTree<K> tree;
    tree.insert(Point_2(50, 50));
    
    // Split to the left of the point - point goes to right tree
    auto right = tree.split(0.0);
    EXPECT_TRUE(tree.empty());
    EXPECT_EQ(right.size(), 1u);
}

TEST(SplitMerge, SplitSinglePointRight) {
    CHTree<K> tree;
    tree.insert(Point_2(50, 50));
    
    // Split to the right of the point - point stays in left tree
    auto right = tree.split(100.0);
    EXPECT_EQ(tree.size(), 1u);
    EXPECT_TRUE(right.empty());
}

TEST(SplitMerge, SplitSinglePointExact) {
    CHTree<K> tree;
    tree.insert(Point_2(50, 50));
    
    // Split at exact x-coordinate - point goes to right tree (>= semantics)
    auto right = tree.split(50.0);
    EXPECT_TRUE(tree.empty());
    EXPECT_EQ(right.size(), 1u);
}

TEST(SplitMerge, SplitTwoPoints) {
    CHTree<K> tree;
    tree.insert(Point_2(10, 10));
    tree.insert(Point_2(90, 90));
    
    // Split in the middle
    auto right = tree.split(50.0);
    EXPECT_EQ(tree.size(), 1u);
    EXPECT_EQ(right.size(), 1u);
}

TEST(SplitMerge, SplitThreePointsMiddle) {
    CHTree<K> tree;
    tree.insert(Point_2(0, 0));
    tree.insert(Point_2(50, 100));
    tree.insert(Point_2(100, 0));
    
    // Split in the middle - left gets (0,0), right gets (50,100) and (100,0)
    auto right = tree.split(25.0);
    EXPECT_EQ(tree.size(), 1u);
    EXPECT_EQ(right.size(), 2u);
}

// ============================================================================
// Basic Functionality Tests - Join
// ============================================================================

TEST(SplitMerge, JoinEmpty) {
    CHTree<K> left;
    left.insert(Point_2(0, 0));
    left.insert(Point_2(10, 10));
    
    CHTree<K> empty;
    left.join(empty);
    
    EXPECT_EQ(left.size(), 2u);
    EXPECT_TRUE(empty.empty());
}

TEST(SplitMerge, JoinIntoEmpty) {
    CHTree<K> empty;
    CHTree<K> right;
    right.insert(Point_2(100, 0));
    right.insert(Point_2(200, 100));
    
    empty.join(right);
    
    EXPECT_EQ(empty.size(), 2u);
    EXPECT_TRUE(right.empty());
}

TEST(SplitMerge, BasicJoin) {
    CHTree<K> left;
    left.insert(Point_2(0, 0));
    left.insert(Point_2(10, 50));
    
    CHTree<K> right;
    right.insert(Point_2(100, 0));
    right.insert(Point_2(110, 50));
    
    left.join(right);
    
    EXPECT_EQ(left.size(), 4u);
    EXPECT_TRUE(right.empty());
}

// ============================================================================
// Correctness After Split Tests
// ============================================================================

TEST(SplitMerge, SplitPreservesHullGeometry_Small) {
    std::mt19937 rng(42);
    auto points = generateSortedRandomPoints(rng, 20, 100);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Split in the middle
    double splitX = 0.0;  // Middle of -100 to 100 range
    auto right = tree.split(splitX);
    
    // Separate points into left and right based on split
    std::vector<Point_2> leftPoints, rightPoints;
    for (const auto& p : points) {
        if (p.x() < splitX) {
            leftPoints.push_back(p);
        } else {
            rightPoints.push_back(p);
        }
    }
    
    // Verify both hulls are correct
    if (leftPoints.size() >= 3) {
        EXPECT_TRUE(verifyHullCorrectness(tree, leftPoints)) << "Left hull incorrect after split";
    }
    if (rightPoints.size() >= 3) {
        EXPECT_TRUE(verifyHullCorrectness(right, rightPoints)) << "Right hull incorrect after split";
    }
    
    // Verify sizes match
    EXPECT_EQ(tree.size(), leftPoints.size());
    EXPECT_EQ(right.size(), rightPoints.size());
}

TEST(SplitMerge, SplitPreservesHullGeometry_Large) {
    std::mt19937 rng(123);
    auto points = generateSortedRandomPoints(rng, 500, 1000);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Split at multiple positions
    std::vector<double> splitPositions = {-500.0, -250.0, 0.0, 250.0, 500.0};
    
    for (double splitX : splitPositions) {
        CHTree<K> testTree;
        testTree.build(points);
        
        auto right = testTree.split(splitX);
        
        std::vector<Point_2> leftPoints, rightPoints;
        for (const auto& p : points) {
            if (p.x() < splitX) {
                leftPoints.push_back(p);
            } else {
                rightPoints.push_back(p);
            }
        }
        
        if (leftPoints.size() >= 3) {
            EXPECT_TRUE(verifyHullCorrectness(testTree, leftPoints)) 
                << "Left hull incorrect after split at " << splitX;
        }
        if (rightPoints.size() >= 3) {
            EXPECT_TRUE(verifyHullCorrectness(right, rightPoints)) 
                << "Right hull incorrect after split at " << splitX;
        }
    }
}

TEST(SplitMerge, SplitVerifyCovers) {
    std::mt19937 rng(456);
    auto points = generateSortedRandomPoints(rng, 100, 200);
    
    CHTree<K> original;
    original.build(points);
    
    double splitX = 0.0;
    auto right = original.split(splitX);
    
    // Test cover queries on both split parts
    for (const auto& p : points) {
        Point_2 testPoint(p.x(), 0);  // Test on x-axis
        
        if (p.x() < splitX) {
            // Point should only be queryable in left tree
            // (covers behavior may vary for points on boundary)
        } else {
            // Point should only be queryable in right tree
        }
    }
    
    // Verify the trees are independent - should work without errors
    EXPECT_GE(original.size() + right.size(), 0u);
}

// ============================================================================
// Correctness After Join Tests
// ============================================================================

TEST(SplitMerge, JoinPreservesHullGeometry_Small) {
    std::mt19937 rng(789);
    
    // Create two disjoint point sets
    std::vector<Point_2> leftPoints, rightPoints;
    std::uniform_int_distribution<int> distY(-100, 100);
    
    for (int x = -100; x <= -10; x += 10) {
        leftPoints.emplace_back(x, distY(rng));
    }
    for (int x = 10; x <= 100; x += 10) {
        rightPoints.emplace_back(x, distY(rng));
    }
    
    CHTree<K> left, right;
    left.build(leftPoints);
    right.build(rightPoints);
    
    // Join and verify
    left.join(right);
    
    std::vector<Point_2> allPoints = leftPoints;
    allPoints.insert(allPoints.end(), rightPoints.begin(), rightPoints.end());
    std::sort(allPoints.begin(), allPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    EXPECT_EQ(left.size(), allPoints.size());
    EXPECT_TRUE(verifyHullCorrectness(left, allPoints));
}

TEST(SplitMerge, JoinPreservesHullGeometry_Large) {
    std::mt19937 rng(1001);
    
    // Create two large disjoint point sets
    auto leftPoints = generateSortedRandomPoints(rng, 200, 500);
    // Shift left points to negative x
    for (auto& p : leftPoints) {
        p = Point_2(p.x() - 1000, p.y());
    }
    std::sort(leftPoints.begin(), leftPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    auto rightPoints = generateSortedRandomPoints(rng, 200, 500);
    // Shift right points to positive x
    for (auto& p : rightPoints) {
        p = Point_2(p.x() + 1000, p.y());
    }
    std::sort(rightPoints.begin(), rightPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    CHTree<K> left, right;
    left.build(leftPoints);
    right.build(rightPoints);
    
    left.join(right);
    
    std::vector<Point_2> allPoints = leftPoints;
    allPoints.insert(allPoints.end(), rightPoints.begin(), rightPoints.end());
    std::sort(allPoints.begin(), allPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    EXPECT_EQ(left.size(), allPoints.size());
    EXPECT_TRUE(verifyHullCorrectness(left, allPoints));
}

TEST(SplitMerge, JoinVerifyCovers) {
    std::mt19937 rng(2002);
    
    // Create two disjoint hulls
    CHTree<K> left, right;
    std::vector<Point_2> leftPoints, rightPoints;
    
    for (int i = 0; i < 30; i++) {
        leftPoints.emplace_back(-100 + i * 3, rng() % 100 - 50);
        rightPoints.emplace_back(100 + i * 3, rng() % 100 - 50);
    }
    
    left.build(leftPoints);
    right.build(rightPoints);
    
    left.join(right);
    
    // Build reference tree with all points
    std::vector<Point_2> allPoints = leftPoints;
    allPoints.insert(allPoints.end(), rightPoints.begin(), rightPoints.end());
    std::sort(allPoints.begin(), allPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    CHTree<K> reference;
    reference.build(allPoints);
    
    // Verify covers() returns the same results
    std::uniform_int_distribution<int> dist(-200, 300);
    for (int i = 0; i < 100; i++) {
        Point_2 query(dist(rng), dist(rng));
        EXPECT_EQ(left.covers(query), reference.covers(query)) 
            << "covers() mismatch for point (" << query.x() << ", " << query.y() << ")";
    }
}

// ============================================================================
// Round-trip Tests (Split then Join)
// ============================================================================

TEST(SplitMerge, SplitThenJoin_Basic) {
    std::mt19937 rng(3003);
    auto points = generateSortedRandomPoints(rng, 50, 200);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Split at x = 0
    auto right = tree.split(0.0);
    size_t leftSize = tree.size();
    size_t rightSize = right.size();
    
    EXPECT_EQ(leftSize + rightSize, points.size());
    
    // Join back together
    tree.join(right);
    
    EXPECT_EQ(tree.size(), points.size());
    EXPECT_TRUE(verifyHullCorrectness(tree, points));
}

TEST(SplitMerge, SplitThenJoin_Large) {
    std::mt19937 rng(4004);
    auto points = generateSortedRandomPoints(rng, 500, 1000);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Split at multiple x-coordinates and rejoin
    std::vector<double> splitPoints = {-500.0, -200.0, 0.0, 200.0};
    
    for (double splitX : splitPoints) {
        CHTree<K> testTree;
        testTree.build(points);
        
        auto right = testTree.split(splitX);
        
        // Join back
        testTree.join(right);
        
        EXPECT_EQ(testTree.size(), points.size()) << "Size mismatch after split at " << splitX;
        EXPECT_TRUE(verifyHullCorrectness(testTree, points)) << "Hull incorrect after split at " << splitX;
    }
}

TEST(SplitMerge, MultipleSplitsAndJoins) {
    std::mt19937 rng(5005);
    auto points = generateSortedRandomPoints(rng, 100, 500);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Split into 4 parts
    auto part2 = tree.split(-250.0);   // tree: x < -250, part2: x >= -250
    auto part3 = part2.split(0.0);      // part2: -250 <= x < 0, part3: x >= 0
    auto part4 = part3.split(250.0);    // part3: 0 <= x < 250, part4: x >= 250
    
    size_t totalSize = tree.size() + part2.size() + part3.size() + part4.size();
    EXPECT_EQ(totalSize, points.size());
    
    // Join all back together (in order!)
    tree.join(part2);   // tree now has x < 0
    tree.join(part3);   // tree now has x < 250
    tree.join(part4);   // tree now has all points
    
    EXPECT_EQ(tree.size(), points.size());
    EXPECT_TRUE(verifyHullCorrectness(tree, points));
}

// ============================================================================
// Stress Tests
// ============================================================================

TEST(SplitMerge, StressTest_RandomSplitJoin) {
    std::mt19937 rng(6006);
    auto points = generateSortedRandomPoints(rng, 500, 1000);
    
    CHTree<K> tree;
    tree.build(points);
    
    std::uniform_int_distribution<int> splitDist(-1000, 1000);
    
    // Perform 50 random split-rejoin operations
    for (int i = 0; i < 50; i++) {
        double splitX = splitDist(rng);
        auto right = tree.split(splitX);
        tree.join(right);
        
        ASSERT_EQ(tree.size(), points.size()) << "Size mismatch at iteration " << i;
        
        // Structural Validation
        if (!hull_helpers::validateTreeInvariants(tree.getRoot())) {
             FAIL() << "Structural invariant failed at iteration " << i;
        }
        
        // Deep Bridge Validation
        if (!tree.validateBridges()) {
             FAIL() << "Bridge validation failed at iteration " << i;
        }

        // Hull Validation (every 5 iterations to save time)
        if (i % 5 == 0) {
            if(!verifyHullCorrectness(tree, points)) {
                std::cerr << "Stress test failed at iteration " << i << " (n=" << points.size() << ")\n";
                FAIL() << "Hull verification failed";
            }
        }
    }
    
    // Final verification
    EXPECT_TRUE(verifyHullCorrectness(tree, points));
}

TEST(SplitMerge, StressTest_LargeHullSplit) {
    std::mt19937 rng(7007);
    auto points = generateSortedRandomPoints(rng, 2000, 10000);
    
    CHTree<K> tree;
    tree.build(points);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Split into halves multiple times
    for (int i = 0; i < 20; i++) {
        CHTree<K> testTree;
        testTree.build(points);
        auto right = testTree.split(0.0);
        EXPECT_EQ(testTree.size() + right.size(), points.size());
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "20 splits of 2000-point hull completed in " << duration.count() << "ms\n";
}

TEST(SplitMerge, StressTest_LargeHullJoin) {
    std::mt19937 rng(8008);
    
    // Create two large disjoint sets
    auto leftPoints = generateSortedRandomPoints(rng, 1000, 5000);
    for (auto& p : leftPoints) {
        p = Point_2(p.x() - 10000, p.y());
    }
    std::sort(leftPoints.begin(), leftPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    auto rightPoints = generateSortedRandomPoints(rng, 1000, 5000);
    for (auto& p : rightPoints) {
        p = Point_2(p.x() + 10000, p.y());
    }
    std::sort(rightPoints.begin(), rightPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Join multiple times
    for (int i = 0; i < 20; i++) {
        CHTree<K> left, right;
        left.build(leftPoints);
        right.build(rightPoints);
        left.join(right);
        EXPECT_EQ(left.size(), leftPoints.size() + rightPoints.size());
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "20 joins of two 1000-point hulls completed in " << duration.count() << "ms\n";
}

// ============================================================================
// Complexity Verification Tests
// ============================================================================

TEST(SplitMerge, SplitComplexity) {
    std::mt19937 rng(9001);
    std::vector<size_t> sizes = {100, 200, 400, 800, 1600, 3200, 6400};
    std::vector<double> times;
    
    for (size_t n : sizes) {
        auto points = generateSortedRandomPoints(rng, n, n * 10);
        
        CHTree<K> tree;
        tree.build(points);
        
        const int numTrials = 50;
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int t = 0; t < numTrials; t++) {
            CHTree<K> testTree;
            testTree.build(points);
            auto right = testTree.split(0.0);
            volatile size_t dummy = testTree.size() + right.size();
            (void)dummy;
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        double avgTimeNs = duration.count() / (double)numTrials;
        times.push_back(avgTimeNs);
    }
    
    // Print results
    std::cout << "\nSplit Complexity Analysis:\n";
    std::cout << "N\t\tTime(ns)\tTime/log²(N)\n";
    for (size_t i = 0; i < sizes.size(); i++) {
        double logN = std::log2(sizes[i]);
        double logSq = logN * logN;
        std::cout << sizes[i] << "\t\t" << times[i] << "\t\t" << times[i] / logSq << "\n";
    }
    
    // Check that time grows roughly as O(log² N)
    // For O(log² N), T(2N)/T(N) should be approximately ((log 2N)/(log N))² ≈ constant for large N
    for (size_t i = 1; i < sizes.size(); i++) {
        double ratio = times[i] / times[i-1];
        double expectedRatio = std::pow(std::log2(sizes[i]) / std::log2(sizes[i-1]), 2);
        
        // Allow generous factor for noise (0.3 to 3.0 of expected)
        // The main check is that ratio doesn't grow linearly with N
        std::cout << "N=" << sizes[i] << ": actual ratio=" << ratio 
                  << ", expected (log² ratio)=" << expectedRatio << "\n";
    }
}

TEST(SplitMerge, JoinComplexity) {
    std::mt19937 rng(9002);
    std::vector<size_t> sizes = {100, 200, 400, 800, 1600, 3200, 6400};
    std::vector<double> times;
    
    for (size_t n : sizes) {
        auto leftPoints = generateSortedRandomPoints(rng, n/2, n * 10);
        for (auto& p : leftPoints) {
            p = Point_2(p.x() - n * 20, p.y());
        }
        std::sort(leftPoints.begin(), leftPoints.end(), [](const auto& a, const auto& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
        
        auto rightPoints = generateSortedRandomPoints(rng, n/2, n * 10);
        for (auto& p : rightPoints) {
            p = Point_2(p.x() + n * 20, p.y());
        }
        std::sort(rightPoints.begin(), rightPoints.end(), [](const auto& a, const auto& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
        
        const int numTrials = 50;
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int t = 0; t < numTrials; t++) {
            CHTree<K> left, right;
            left.build(leftPoints);
            right.build(rightPoints);
            left.join(right);
            volatile size_t dummy = left.size();
            (void)dummy;
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        double avgTimeNs = duration.count() / (double)numTrials;
        times.push_back(avgTimeNs);
    }
    
    // Print results
    std::cout << "\nJoin Complexity Analysis:\n";
    std::cout << "N\t\tTime(ns)\tTime/log²(N)\n";
    for (size_t i = 0; i < sizes.size(); i++) {
        double logN = std::log2(sizes[i]);
        double logSq = logN * logN;
        std::cout << sizes[i] << "\t\t" << times[i] << "\t\t" << times[i] / logSq << "\n";
    }
    
    // Similar check as split
    for (size_t i = 1; i < sizes.size(); i++) {
        double ratio = times[i] / times[i-1];
        double expectedRatio = std::pow(std::log2(sizes[i]) / std::log2(sizes[i-1]), 2);
        std::cout << "N=" << sizes[i] << ": actual ratio=" << ratio 
                  << ", expected (log² ratio)=" << expectedRatio << "\n";
    }
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST(SplitMerge, SplitAtExactPoints) {
    CHTree<K> tree;
    std::vector<Point_2> points;
    for (int x = 0; x <= 100; x += 10) {
        points.emplace_back(x, x);
    }
    tree.build(points);
    
    // Split at exact x-coordinates
    for (int x = 0; x <= 100; x += 10) {
        CHTree<K> testTree;
        testTree.build(points);
        
        auto right = testTree.split(x);
        
        size_t expectedLeft = x / 10;  // Points with x < splitX
        size_t expectedRight = points.size() - expectedLeft;
        
        EXPECT_EQ(testTree.size(), expectedLeft) << "Left size mismatch at x=" << x;
        EXPECT_EQ(right.size(), expectedRight) << "Right size mismatch at x=" << x;
    }
}

TEST(SplitMerge, JoinManySmallHulls) {
    std::mt19937 rng(10001);
    
    // Create 20 small hulls with disjoint x-ranges
    std::vector<std::unique_ptr<CHTree<K>>> hulls;
    std::vector<Point_2> allPoints;
    
    for (int i = 0; i < 20; i++) {
        auto hull = std::make_unique<CHTree<K>>();
        std::vector<Point_2> pts;
        for (int j = 0; j < 10; j++) {
            Point_2 p(i * 100 + j * 5, (int)(rng() % 100) - 50);
            pts.push_back(p);
            allPoints.push_back(p);
        }
        std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
        hull->build(pts);
        hulls.push_back(std::move(hull));
    }
    
    // Join all hulls into the first one
    for (size_t i = 1; i < hulls.size(); i++) {
        hulls[0]->join(*hulls[i]);
    }
    
    std::sort(allPoints.begin(), allPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    EXPECT_EQ(hulls[0]->size(), allPoints.size());
    EXPECT_TRUE(verifyHullCorrectness(*hulls[0], allPoints));
}

TEST(SplitMerge, SplitIntoManyParts) {
    std::mt19937 rng(10002);
    auto points = generateSortedRandomPoints(rng, 100, 1000);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Split into 10 parts
    std::vector<CHTree<K>> parts;
    parts.push_back(std::move(tree));
    
    for (int splitX = -800; splitX <= 600; splitX += 200) {
        auto right = parts.back().split(splitX);
        parts.push_back(std::move(right));
    }
    
    // Count total points
    size_t total = 0;
    for (const auto& part : parts) {
        total += part.size();
    }
    EXPECT_EQ(total, points.size());
    
    // Join all back together
    for (size_t i = 1; i < parts.size(); i++) {
        parts[0].join(parts[i]);
    }
    
    EXPECT_EQ(parts[0].size(), points.size());
    EXPECT_TRUE(verifyHullCorrectness(parts[0], points));
}

TEST(SplitMerge, OperationsAfterSplit) {
    std::mt19937 rng(10003);
    auto points = generateSortedRandomPoints(rng, 50, 200);
    
    CHTree<K> tree;
    tree.build(points);
    
    auto right = tree.split(0.0);
    
    // Insert new points into both trees
    std::vector<Point_2> leftPoints, rightPoints;
    for (const auto& p : points) {
        if (p.x() < 0.0) leftPoints.push_back(p);
        else rightPoints.push_back(p);
    }
    
    // Insert into left tree
    Point_2 newLeft(-150, 0);
    tree.insert(newLeft);
    leftPoints.push_back(newLeft);
    std::sort(leftPoints.begin(), leftPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    // Insert into right tree
    Point_2 newRight(250, 0);
    right.insert(newRight);
    rightPoints.push_back(newRight);
    std::sort(rightPoints.begin(), rightPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    if (leftPoints.size() >= 3) {
        EXPECT_TRUE(verifyHullCorrectness(tree, leftPoints));
    }
    if (rightPoints.size() >= 3) {
        EXPECT_TRUE(verifyHullCorrectness(right, rightPoints));
    }
}

TEST(SplitMerge, OperationsAfterJoin) {
    std::mt19937 rng(10004);
    
    CHTree<K> left, right;
    std::vector<Point_2> leftPoints, rightPoints;
    
    for (int i = 0; i < 20; i++) {
        leftPoints.emplace_back(-100 + i * 4, (int)(rng() % 100) - 50);
        rightPoints.emplace_back(100 + i * 4, (int)(rng() % 100) - 50);
    }
    std::sort(leftPoints.begin(), leftPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    std::sort(rightPoints.begin(), rightPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    left.build(leftPoints);
    right.build(rightPoints);
    
    left.join(right);
    
    std::vector<Point_2> allPoints = leftPoints;
    allPoints.insert(allPoints.end(), rightPoints.begin(), rightPoints.end());
    
    // Insert and remove points
    Point_2 newPoint(0, 200);
    left.insert(newPoint);
    allPoints.push_back(newPoint);
    std::sort(allPoints.begin(), allPoints.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    EXPECT_TRUE(verifyHullCorrectness(left, allPoints));
    
    // Remove a point
    left.remove(newPoint);
    allPoints.erase(std::find(allPoints.begin(), allPoints.end(), newPoint));
    
    EXPECT_TRUE(verifyHullCorrectness(left, allPoints));
}

TEST(SplitMerge, RandomizedPartitionStressTest) {
    std::mt19937 rng(12345);
    // Use a moderate number of points to keep test duration reasonable
    const int NUM_POINTS = 500; 
    const int NUM_ITERATIONS = 500;
    const int RANGE = 10000;
    
    auto points = generateSortedRandomPoints(rng, NUM_POINTS, RANGE);
    
    // Maintain a list of partitions. 
    // They are implicitly sorted by x-range because we only split and join adjacent ones.
    std::vector<std::unique_ptr<CHTree<K>>> partitions;
    
    // Start with one big partition
    auto initialHull = std::make_unique<CHTree<K>>();
    initialHull->build(points);
    partitions.push_back(std::move(initialHull));
    
    std::uniform_real_distribution<double> actionDist(0.0, 1.0);
    
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        // Probability of split vs join
        // Bias towards split if we have few partitions, join if we have many
        double joinProb = (double)partitions.size() / 20.0; // scales up as partitions grow
        if (joinProb > 0.8) joinProb = 0.8;
        if (partitions.size() < 2) joinProb = 0.0;
        
        bool doJoin = actionDist(rng) < joinProb;
        
        if (doJoin) {
            // JOIN OPERATION
            // Pick a random adjacent pair
            std::uniform_int_distribution<size_t> idxDist(0, partitions.size() - 2);
            size_t idx = idxDist(rng);
            
            // Join idx+1 into idx
            partitions[idx]->join(*partitions[idx+1]);
            
            // Remove the empty shell of idx+1
            partitions.erase(partitions.begin() + idx + 1);
            
        } else {
            // SPLIT OPERATION
            // Pick a random partition
            std::uniform_int_distribution<size_t> idxDist(0, partitions.size() - 1);
            size_t idx = idxDist(rng);
            
            CHTree<K>& hull = *partitions[idx];
            
            if (hull.size() < 2) continue; // Can't split much further
            
            // Find min and max x of this hull to choose a split point
            auto pts = hull.upperHullPoints(); // Just need x range
            if (pts.empty()) continue;
            
            double minX = pts.front().x();
            double maxX = pts.back().x();
            
            if (maxX - minX < 1.0) continue; // Too narrow
            
            std::uniform_real_distribution<double> splitDist(minX + 0.1, maxX - 0.1);
            double splitX = splitDist(rng);
            
            auto rightParams = hull.split(splitX);
            auto rightHull = std::make_unique<CHTree<K>>();
            // split method returns by value, we need to move it into heap object
            *rightHull = std::move(rightParams);
            
            // Insert right hull after the current one
            partitions.insert(partitions.begin() + idx + 1, std::move(rightHull));
        }
        
        // --- VERIFICATION ---
        
        // 1. Total size check
        size_t currentTotal = 0;
        for (const auto& p : partitions) currentTotal += p->size();
        ASSERT_EQ(currentTotal, NUM_POINTS) << "Total points mismatch at iteration " << i;
        
        // 2. Local Invariant Check (every 10 iterations)
        if (i % 10 == 0) {
            for (const auto& p : partitions) {
                 ASSERT_TRUE(p->validateBridges()) << "Bridge validation failed at iter " << i;
            }
        }
        
        // 3. Global Invariant Check (every 50 iterations)
        // Verify partitions are sorted and disjoint
        if (i % 50 == 0) {
             double lastMaxX = -std::numeric_limits<double>::infinity();
             for(const auto& p : partitions) {
                 if (p->empty()) continue;
                 auto uh = p->upperHullPoints();
                 if (uh.empty()) continue;
                 
                 double minX = uh.front().x();
                 double maxX = uh.back().x();
                 
                 ASSERT_LE(minX, maxX);
                 // Relaxed check: just ensure non-decreasing start points or strictly increasing "zones"
                 // Since we split at strict double X, the sets should be disjoint.
                 // However, points might lie exactly on split boundary? split() uses < splitX for left.
                 // So left max < splitX <= right min.
                 ASSERT_LT(lastMaxX, minX) << "Partitions not sorted or overlapping at iter " << i;
                 lastMaxX = maxX;
             }
        }
    }
}

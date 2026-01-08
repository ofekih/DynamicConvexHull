// Minimal test case to debug split/join
// Tests various sizes to find threshold where it breaks

#include <gtest/gtest.h>
#include <CGAL/Simple_cartesian.h>
#include "CHTree.h"
#include "hull_test_helpers.hpp"
#include <random>
#include <iostream>

using K = CGAL::Simple_cartesian<double>;
using Point_2 = K::Point_2;

// Generate n sorted random points with seed
std::vector<Point_2> genPoints(int n, int seed, int range = 1000) {
    std::mt19937 rng(seed);
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

// Verify hull correctness
bool verifyHull(CHTree<K>& tree, const std::vector<Point_2>& points) {
    if (points.size() < 3) return true;
    
    auto expected_upper = hull_helpers::adjustUpperHullForCHTree(
        monotone_chain::upperHull(hull_helpers::toIntPairs(points)));
    auto expected_lower = hull_helpers::adjustLowerHullForCHTree(
        monotone_chain::lowerHull(hull_helpers::toIntPairs(points)));
    
    auto ch_upper = hull_helpers::toIntPairs(tree.upperHullPoints());
    auto ch_lower = hull_helpers::toIntPairs(tree.lowerHullPoints());
    
    return hull_helpers::hullContainsAll(expected_upper, ch_upper) &&
           hull_helpers::hullContainsAll(expected_lower, ch_lower);
}

// Print hull points for debugging
void printHulls(CHTree<K>& tree, const std::vector<Point_2>& points) {
    std::cerr << "Expected upper hull: ";
    auto exp_upper = monotone_chain::upperHull(hull_helpers::toIntPairs(points));
    for (auto& p : exp_upper) std::cerr << "(" << p.first << "," << p.second << ") ";
    std::cerr << "\nActual upper hull: ";
    auto act_upper = tree.upperHullPoints();
    for (auto& p : act_upper) std::cerr << "(" << p.x() << "," << p.y() << ") ";
    std::cerr << "\n";
    
    std::cerr << "Expected lower hull: ";
    auto exp_lower = monotone_chain::lowerHull(hull_helpers::toIntPairs(points));
    for (auto& p : exp_lower) std::cerr << "(" << p.first << "," << p.second << ") ";
    std::cerr << "\nActual lower hull: ";
    auto act_lower = tree.lowerHullPoints();
    for (auto& p : act_lower) std::cerr << "(" << p.x() << "," << p.y() << ") ";
    std::cerr << "\n";
}

// Test split at x=0 for various sizes
TEST(Debug, FindMinimalFailingSize) {
    int seed = 123;
    double splitX = 0.0;
    
    for (int n = 10; n <= 100; n += 5) {
        auto points = genPoints(n, seed);
        CHTree<K> tree;
        tree.build(points);
        
        auto right = tree.split(splitX);
        
        // Separate points
        std::vector<Point_2> leftPoints, rightPoints;
        for (const auto& p : points) {
            if (p.x() < splitX) leftPoints.push_back(p);
            else rightPoints.push_back(p);
        }
        
        bool leftOK = leftPoints.size() < 3 || verifyHull(tree, leftPoints);
        bool rightOK = rightPoints.size() < 3 || verifyHull(right, rightPoints);
        
        if (!leftOK || !rightOK) {
            std::cerr << "FAILED at n=" << n << " (left=" << leftPoints.size() 
                      << ", right=" << rightPoints.size() << ")\n";
            
            if (!leftOK) {
                std::cerr << "Left tree failure:\n";
                printHulls(tree, leftPoints);
            }
            if (!rightOK) {
                std::cerr << "Right tree failure:\n";
                printHulls(right, rightPoints);
            }
            FAIL() << "Failed at n=" << n;
        } else {
            std::cerr << "OK at n=" << n << "\n";
        }
    }
}

// Test split with specific small failing case
TEST(Debug, SmallestFailingCase) {
    int n = 25;
    int seed = 123;
    double splitX = 0.0;
    
    auto points = genPoints(n, seed);
    
    std::cerr << "=== Points (n=" << n << ", splitX=" << splitX << ") ===\n";
    for (size_t i = 0; i < points.size(); i++) {
        std::cerr << i << ": (" << points[i].x() << ", " << points[i].y() << ")";
        if (points[i].x() < splitX) std::cerr << " -> LEFT\n";
        else std::cerr << " -> RIGHT\n";
    }
    
    CHTree<K> tree;
    tree.build(points);
    
    // Separate expected points
    std::vector<Point_2> leftPoints, rightPoints;
    for (const auto& p : points) {
        if (p.x() < splitX) leftPoints.push_back(p);
        else rightPoints.push_back(p);
    }
    
    std::cerr << "\n=== Before split ===\n";
    std::cerr << "Tree size: " << tree.size() << "\n";
    std::cerr << "Expected left: " << leftPoints.size() << ", right: " << rightPoints.size() << "\n";
    
    auto right = tree.split(splitX);
    
    std::cerr << "\n=== After split ===\n";
    std::cerr << "Left tree size: " << tree.size() << " (expected " << leftPoints.size() << ")\n";
    std::cerr << "Right tree size: " << right.size() << " (expected " << rightPoints.size() << ")\n";
    

    
    // Check sizes first
    ASSERT_EQ(tree.size(), leftPoints.size()) << "Left tree has wrong number of points!";
    ASSERT_EQ(right.size(), rightPoints.size()) << "Right tree has wrong number of points!";
    
    // Structural Validation
    std::cerr << "\n=== Structural Validation ===\n";
    if (!hull_helpers::validateTreeInvariants(tree.getRoot())) {
        FAIL() << "Left tree (tree) validation failed!";
    } else {
        std::cerr << "Left tree structure valid.\n";
    }
    
    if (!hull_helpers::validateTreeInvariants(right.getRoot())) {
        FAIL() << "Right tree (right) validation failed!";
    } else {
        std::cerr << "Right tree structure valid.\n";
    }
    
    if (leftPoints.size() >= 3) {
        std::cerr << "\n=== Left hull verification ===\n";
        printHulls(tree, leftPoints);
        EXPECT_TRUE(verifyHull(tree, leftPoints)) << "Left hull incorrect";
    }
    
    if (rightPoints.size() >= 3) {
        std::cerr << "\n=== Right hull verification ===\n";
        printHulls(right, rightPoints);
        EXPECT_TRUE(verifyHull(right, rightPoints)) << "Right hull incorrect";
    }
}

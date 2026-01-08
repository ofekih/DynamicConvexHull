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

// Helper to print tree structure recursively
template<typename Node>
void printTreeStructure(Node* node, int depth = 0, const std::string& prefix = "") {
    if (!node) return;
    
    std::string indent(depth * 2, ' ');
    bool is_leaf_node = !node->left && !node->right;
    
    if (is_leaf_node) {
        // Leaf: print point
        auto p = node->val[0].min();
        std::cerr << indent << prefix << "LEAF: (" << p.x() << ", " << p.y() << ")\n";
    } else {
        // Internal: print bridge info
        auto ub = node->val[0]; // val[0] = upper bridge (findBridge<false>)
        auto lb = node->val[1]; // val[1] = lower bridge (findBridge<true>)
        std::cerr << indent << prefix << "INTERNAL [" << node->size << " pts, h=" << node->height << "]\n";
        std::cerr << indent << "  UpperBridge: (" << ub.min().x() << "," << ub.min().y() << ") -> (" 
                  << ub.max().x() << "," << ub.max().y() << ")\n";
        std::cerr << indent << "  LowerBridge: (" << lb.min().x() << "," << lb.min().y() << ") -> (" 
                  << lb.max().x() << "," << lb.max().y() << ")\n";
        std::cerr << indent << "  Range: [" << node->min[0].min().x() << " .. " << node->max[0].max().x() << "]\n";
        
        printTreeStructure(node->left, depth + 1, "L: ");
        printTreeStructure(node->right, depth + 1, "R: ");
    }
}

// Find minimum failing stress test case
TEST(Debug, FindMinimalStressFailure) {
    std::cerr << "\n=== Searching for minimum failing stress test case ===\n";
    
    // Try different seeds and sizes
    for (int seed = 42; seed <= 52; seed++) {
        for (int n = 10; n <= 100; n += 10) {
            std::mt19937 rng(seed);
            std::uniform_int_distribution<int> dist(-1000, 1000);
            
            // Generate points
            std::vector<Point_2> points;
            for (int i = 0; i < n; i++) {
                points.emplace_back(dist(rng), dist(rng));
            }
            std::sort(points.begin(), points.end(), [](const Point_2& a, const Point_2& b) {
                if (a.x() != b.x()) return a.x() < b.x();
                return a.y() < b.y();
            });
            
            CHTree<K> tree;
            tree.build(points);
            
            // Check initial build
            if (!tree.validateBridges()) {
                std::cerr << "FOUND: Build failed at seed=" << seed << ", n=" << n << "\n";
                FAIL() << "Initial build has invalid bridges";
            }
            
            std::uniform_int_distribution<int> splitDist(-1000, 1000);
            
            // Perform split/join operations
            for (int iter = 0; iter < 20; iter++) {
                double splitX = splitDist(rng);
                
                // Validate before split
                bool validBefore = tree.validateBridges();
                bool correctBefore = verifyHull(tree, points);
                
                auto right = tree.split(splitX);
                
                // Validate after split
                bool validAfterLeft = tree.validateBridges();
                bool validAfterRight = right.validateBridges();
                
                tree.join(right);
                
                // Validate after join
                bool validAfterJoin = tree.validateBridges();
                bool correctAfterJoin = verifyHull(tree, points);
                
                // Check for first failure
                if (!validAfterJoin || !correctAfterJoin) {
                    std::cerr << "\n=== FOUND MINIMAL FAILURE ===\n";
                    std::cerr << "seed=" << seed << ", n=" << n << ", iter=" << iter << ", splitX=" << splitX << "\n";
                    std::cerr << "validBefore=" << validBefore << ", correctBefore=" << correctBefore << "\n";
                    std::cerr << "validAfterLeft=" << validAfterLeft << ", validAfterRight=" << validAfterRight << "\n";
                    std::cerr << "validAfterJoin=" << validAfterJoin << ", correctAfterJoin=" << correctAfterJoin << "\n";
                    
                    if (!validBefore || !correctBefore) {
                        std::cerr << "\n*** Issue existed BEFORE this operation! ***\n";
                    } else {
                        std::cerr << "\n*** This operation caused the issue ***\n";
                    }
                    
                    // Print hull comparison
                    std::cerr << "\nHull comparison:\n";
                    printHulls(tree, points);
                    
                    FAIL() << "Found minimal failure at seed=" << seed << ", n=" << n << ", iter=" << iter;
                }
            }
        }
    }
    std::cerr << "All combinations passed!\n";
}

// Test specific stress case with detailed state printing
TEST(Debug, DetailedStressDebug) {
    // Use the same parameters as StressTest_RandomSplitJoin
    int seed = 5005;
    int n = 500;
    
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(-1000, 1000);
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    std::sort(points.begin(), points.end(), [](const Point_2& a, const Point_2& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    CHTree<K> tree;
    tree.build(points);
    
    std::uniform_int_distribution<int> splitDist(-1000, 1000);
    
    std::cerr << "\n=== Starting detailed stress debug ===\n";
    std::cerr << "seed=" << seed << ", n=" << n << "\n";
    
    for (int iter = 0; iter < 50; iter++) {
        double splitX = splitDist(rng);
        
        // Check state BEFORE operation
        bool bridgesValidBefore = tree.validateBridges();
        bool hullCorrectBefore = verifyHull(tree, points);
        
        std::cerr << "Iter " << iter << ": splitX=" << splitX 
                  << ", bridgesOK=" << bridgesValidBefore 
                  << ", hullOK=" << hullCorrectBefore << "\n";
        
        if (!bridgesValidBefore) {
            std::cerr << "*** Bridges already invalid BEFORE iter " << iter << "! ***\n";
            std::cerr << "Tree structure:\n";
            printTreeStructure(tree.getRoot());
            FAIL() << "Bridges invalid before iter " << iter;
        }
        
        // Perform split
        auto right = tree.split(splitX);
        
        // Check after split
        bool leftBridgesOK = tree.validateBridges();
        bool rightBridgesOK = right.validateBridges();
        
        if (!leftBridgesOK || !rightBridgesOK) {
            std::cerr << "*** Split at iter " << iter << " corrupted bridges! ***\n";
            std::cerr << "leftBridgesOK=" << leftBridgesOK << ", rightBridgesOK=" << rightBridgesOK << "\n";
            if (!leftBridgesOK) {
                std::cerr << "Left tree structure:\n";
                printTreeStructure(tree.getRoot());
            }
            if (!rightBridgesOK) {
                std::cerr << "Right tree structure:\n";
                printTreeStructure(right.getRoot());
            }
        }
        
        // Perform join
        tree.join(right);
        
        // Check after join
        bool bridgesValidAfter = tree.validateBridges();
        bool hullCorrectAfter = verifyHull(tree, points);
        
        if (!bridgesValidAfter || !hullCorrectAfter) {
            std::cerr << "*** Join at iter " << iter << " failed! ***\n";
            std::cerr << "bridgesValidAfter=" << bridgesValidAfter << ", hullCorrectAfter=" << hullCorrectAfter << "\n";
            
            if (!hullCorrectAfter) {
                std::cerr << "\nHull mismatch:\n";
                printHulls(tree, points);
            }
            
            if (!bridgesValidAfter) {
                std::cerr << "Tree structure after failed join:\n";
                printTreeStructure(tree.getRoot());
            }
            
            FAIL() << "Failed at iter " << iter;
        }
    }
    
    std::cerr << "All 50 iterations passed!\n";
}

// Targeted test for the minimal failing case found
TEST(Debug, TargetedMinimalCase) {
    // Exact failing parameters
    int seed = 46;
    int n = 50;
    int numIterations = 3; // Fail happens at iter=2
    
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(-1000, 1000);
    
    // Generate points
    std::vector<Point_2> points;
    for (int i = 0; i < n; i++) {
        points.emplace_back(dist(rng), dist(rng));
    }
    std::sort(points.begin(), points.end(), [](const Point_2& a, const Point_2& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    
    std::cerr << "\n=== TARGETED MINIMAL CASE ===\n";
    std::cerr << "seed=" << seed << ", n=" << n << "\n";
    std::cerr << "\nPoints:\n";
    for (size_t i = 0; i < points.size(); i++) {
        std::cerr << i << ": (" << points[i].x() << ", " << points[i].y() << ")\n";
    }
    
    CHTree<K> tree;
    tree.build(points);
    
    std::cerr << "\n=== Initial tree ===\n";
    printTreeStructure(tree.getRoot());
    
    std::uniform_int_distribution<int> splitDist(-1000, 1000);
    
    for (int iter = 0; iter < numIterations; iter++) {
        double splitX = splitDist(rng);
        
        std::cerr << "\n==============================\n";
        std::cerr << "=== ITERATION " << iter << ": splitX=" << splitX << " ===\n";
        std::cerr << "==============================\n";
        
        // State before split
        std::cerr << "\n--- BEFORE SPLIT ---\n";
        std::cerr << "Tree size: " << tree.size() << "\n";
        std::cerr << "validateBridges: " << tree.validateBridges() << "\n";
        std::cerr << "verifyHull: " << verifyHull(tree, points) << "\n";
        
        // Perform split
        auto right = tree.split(splitX);
        
        std::cerr << "\n--- AFTER SPLIT ---\n";
        std::cerr << "Left size: " << tree.size() << ", Right size: " << right.size() << "\n";
        std::cerr << "Left validateBridges: " << tree.validateBridges() << "\n";
        std::cerr << "Right validateBridges: " << right.validateBridges() << "\n";
        
        // State before join
        std::cerr << "\n--- BEFORE JOIN ---\n";
        std::cerr << "Left tree:\n";
        printTreeStructure(tree.getRoot());
        std::cerr << "\nRight tree:\n";
        printTreeStructure(right.getRoot());
        
        // Perform join
        tree.join(right);
        
        std::cerr << "\n--- AFTER JOIN ---\n";
        std::cerr << "Tree size: " << tree.size() << "\n";
        std::cerr << "validateBridges: " << tree.validateBridges() << "\n";
        
        bool correct = verifyHull(tree, points);
        std::cerr << "verifyHull: " << correct << "\n";
        
        if (!correct) {
            std::cerr << "\n*** HULL FAILED AT ITER " << iter << " ***\n";
            std::cerr << "\nTree structure after failed join:\n";
            printTreeStructure(tree.getRoot());
            
            std::cerr << "\nHull comparison:\n";
            printHulls(tree, points);
            
            FAIL() << "Hull verification failed at iter " << iter;
        }
    }
    
    std::cerr << "\nAll iterations passed!\n";
}



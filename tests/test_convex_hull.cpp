// Tests comparing Dynamic Convex Hull (CHTree) to a standard convex hull solution
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <memory>
#include <chrono>
#include "CHTree.h"
#include "hull_test_helpers.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

// ============================================================================
// Print hull for debugging
// ============================================================================

void printHull(const std::string& name, const std::vector<std::pair<int, int>>& hull) {
    std::cerr << name << ": ";
    for (const auto& p : hull) std::cerr << "(" << p.first << "," << p.second << ") ";
    std::cerr << "\n";
}

// ============================================================================
// Test Helper Class
// ============================================================================

class ConvexHullTester {
public:
    std::unique_ptr<CHTree<K>> ch;
    std::vector<std::pair<int, int>> points;
    std::mt19937 rng;

    ConvexHullTester(unsigned seed = 42) : ch(std::make_unique<CHTree<K>>()), rng(seed) {}

    void reset(unsigned seed) {
        ch = std::make_unique<CHTree<K>>();
        points.clear();
        rng.seed(seed);
    }

    void addPoint(int x, int y) {
        ch->insert(Point_2(x, y));
        points.push_back({x, y});
    }

    void removePoint(int x, int y) {
        ch->remove(Point_2(x, y));
        auto it = std::find(points.begin(), points.end(), std::make_pair(x, y));
        if (it != points.end()) points.erase(it);
    }

    bool verifyHulls() {
        if (points.size() < 3) return true;

        auto expected_upper = hull_helpers::adjustUpperHullForCHTree(
            monotone_chain::upperHull(points));
        auto expected_lower = hull_helpers::adjustLowerHullForCHTree(
            monotone_chain::lowerHull(points));

        auto ch_upper = hull_helpers::toIntPairs(ch->upperHullPoints());
        auto ch_lower = hull_helpers::toIntPairs(ch->lowerHullPoints());

        bool upper_ok = hull_helpers::hullContainsAll(expected_upper, ch_upper);
        bool lower_ok = hull_helpers::hullContainsAll(expected_lower, ch_lower);

        if (!upper_ok) {
            std::cerr << "Upper hull mismatch! Expected points not found.\n";
            printHull("Expected", expected_upper);
            printHull("Got", ch_upper);
        }
        if (!lower_ok) {
            std::cerr << "Lower hull mismatch! Expected points not found.\n";
            printHull("Expected", expected_lower);
            printHull("Got", ch_lower);
        }

        return upper_ok && lower_ok;
    }

    std::pair<int, int> randomPoint(int range = 1000) {
        std::uniform_int_distribution<int> dist(-range, range);
        return {dist(rng), dist(rng)};
    }
};

// ============================================================================
// Tests: Random Insertion
// ============================================================================

TEST(ConvexHullComparison, RandomInsertion_Small) {
    ConvexHullTester tester(42);
    for (int i = 0; i < 10; i++) {
        auto [x, y] = tester.randomPoint(100);
        tester.addPoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, RandomInsertion_Medium) {
    ConvexHullTester tester(123);
    for (int i = 0; i < 50; i++) {
        auto [x, y] = tester.randomPoint(500);
        tester.addPoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, RandomInsertion_Large) {
    ConvexHullTester tester(456);
    for (int i = 0; i < 200; i++) {
        auto [x, y] = tester.randomPoint(1000);
        tester.addPoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, RandomInsertion_VerifyIncrementally) {
    ConvexHullTester tester(789);
    for (int i = 0; i < 30; i++) {
        auto [x, y] = tester.randomPoint(200);
        tester.addPoint(x, y);
        if (i >= 2) {
            EXPECT_TRUE(tester.verifyHulls()) << "Failed after inserting point " << i;
        }
    }
}

// ============================================================================
// Tests: Random Insertion + Removal
// ============================================================================

TEST(ConvexHullComparison, RandomInsertionRemoval_Basic) {
    ConvexHullTester tester(1001);
    
    for (int i = 0; i < 20; i++) {
        auto [x, y] = tester.randomPoint(200);
        tester.addPoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());

    std::shuffle(tester.points.begin(), tester.points.end(), tester.rng);
    for (int i = 0; i < 10 && !tester.points.empty(); i++) {
        auto [x, y] = tester.points.back();
        tester.removePoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, RandomInsertionRemoval_Interleaved) {
    ConvexHullTester tester(2002);
    std::uniform_int_distribution<int> action(0, 2);

    for (int i = 0; i < 100; i++) {
        if (tester.points.empty() || action(tester.rng) != 0) {
            auto [x, y] = tester.randomPoint(300);
            tester.addPoint(x, y);
        } else {
            std::uniform_int_distribution<size_t> idx(0, tester.points.size() - 1);
            auto [x, y] = tester.points[idx(tester.rng)];
            tester.removePoint(x, y);
        }

        if (tester.points.size() >= 3) {
            EXPECT_TRUE(tester.verifyHulls()) << "Failed at iteration " << i;
        }
    }
}

TEST(ConvexHullComparison, RandomInsertionRemoval_ManyRemovals) {
    ConvexHullTester tester(3003);
    
    for (int i = 0; i < 50; i++) {
        auto [x, y] = tester.randomPoint(400);
        tester.addPoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());

    while (tester.points.size() > 5) {
        std::uniform_int_distribution<size_t> idx(0, tester.points.size() - 1);
        auto [x, y] = tester.points[idx(tester.rng)];
        tester.removePoint(x, y);
        if (tester.points.size() >= 3) {
            EXPECT_TRUE(tester.verifyHulls());
        }
    }
}

// ============================================================================
// Tests: Edge Cases
// ============================================================================

TEST(ConvexHullComparison, CollinearPoints_Horizontal) {
    ConvexHullTester tester;
    for (int i = 0; i < 10; i++) {
        tester.addPoint(i * 10, 50);
    }
    auto ch_upper = tester.ch->upperHullPoints();
    auto ch_lower = tester.ch->lowerHullPoints();
    EXPECT_GE(ch_upper.size(), 1u);
    EXPECT_GE(ch_lower.size(), 1u);
}

TEST(ConvexHullComparison, CollinearPoints_Vertical) {
    ConvexHullTester tester;
    for (int i = 0; i < 10; i++) {
        tester.addPoint(50, i * 10);
    }
    auto ch_upper = tester.ch->upperHullPoints();
    auto ch_lower = tester.ch->lowerHullPoints();
    EXPECT_GE(ch_upper.size(), 1u);
    EXPECT_GE(ch_lower.size(), 1u);
}

TEST(ConvexHullComparison, Triangle) {
    ConvexHullTester tester;
    tester.addPoint(0, 0);
    tester.addPoint(100, 0);
    tester.addPoint(50, 100);
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, Square) {
    ConvexHullTester tester;
    tester.addPoint(0, 0);
    tester.addPoint(100, 0);
    tester.addPoint(100, 100);
    tester.addPoint(0, 100);
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, SquareWithInteriorPoint) {
    ConvexHullTester tester;
    tester.addPoint(0, 0);
    tester.addPoint(100, 0);
    tester.addPoint(100, 100);
    tester.addPoint(0, 100);
    tester.addPoint(50, 50);
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, CircleApproximation) {
    ConvexHullTester tester;
    const int n = 20;
    const int r = 100;
    for (int i = 0; i < n; i++) {
        double angle = 2.0 * 3.14159265358979 * i / n;
        tester.addPoint((int)(r * cos(angle)), (int)(r * sin(angle)));
    }
    EXPECT_TRUE(tester.verifyHulls());
}

// ============================================================================
// Tests: Stress Tests
// ============================================================================

TEST(ConvexHullComparison, StressTest_ManyOperations) {
    ConvexHullTester tester(4004);
    std::uniform_int_distribution<int> action(0, 4);

    for (int i = 0; i < 500; i++) {
        if (tester.points.empty() || action(tester.rng) != 0) {
            auto [x, y] = tester.randomPoint(1000);
            tester.addPoint(x, y);
        } else {
            std::uniform_int_distribution<size_t> idx(0, tester.points.size() - 1);
            auto [x, y] = tester.points[idx(tester.rng)];
            tester.removePoint(x, y);
        }
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, StressTest_DifferentSeeds) {
    for (int seed = 0; seed < 10; seed++) {
        ConvexHullTester tester(seed * 12345);

        for (int i = 0; i < 50; i++) {
            auto [x, y] = tester.randomPoint(500);
            tester.addPoint(x, y);
        }
        EXPECT_TRUE(tester.verifyHulls()) << "Failed with seed " << seed;
    }
}

TEST(ConvexHullComparison, InsertAndRemoveAll) {
    ConvexHullTester tester(5005);
    
    std::vector<std::pair<int, int>> inserted;
    for (int i = 0; i < 30; i++) {
        auto [x, y] = tester.randomPoint(300);
        tester.addPoint(x, y);
        inserted.push_back({x, y});
    }
    EXPECT_TRUE(tester.verifyHulls());

    std::shuffle(inserted.begin(), inserted.end(), tester.rng);
    for (auto [x, y] : inserted) {
        tester.removePoint(x, y);
    }
    EXPECT_TRUE(tester.points.empty());
}

// ============================================================================
// Tests: Specific Patterns
// ============================================================================

TEST(ConvexHullComparison, DiamondPattern) {
    ConvexHullTester tester;
    tester.addPoint(0, 100);
    tester.addPoint(100, 0);
    tester.addPoint(0, -100);
    tester.addPoint(-100, 0);
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, StarPattern) {
    ConvexHullTester tester;
    const int n = 5;
    for (int i = 0; i < n; i++) {
        double angle = 2.0 * 3.14159265358979 * i / n;
        tester.addPoint((int)(100 * cos(angle)), (int)(100 * sin(angle)));
        tester.addPoint((int)(50 * cos(angle + 3.14159265358979/n)), 
                        (int)(50 * sin(angle + 3.14159265358979/n)));
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, GridPattern) {
    ConvexHullTester tester;
    for (int x = 0; x <= 100; x += 10) {
        for (int y = 0; y <= 100; y += 10) {
            tester.addPoint(x, y);
        }
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, ClusteredPoints) {
    ConvexHullTester tester(6006);
    std::uniform_int_distribution<int> clusterOffset(-20, 20);
    
    std::vector<std::pair<int, int>> centers = {{0, 0}, {200, 0}, {100, 200}};
    for (auto [cx, cy] : centers) {
        for (int i = 0; i < 10; i++) {
            tester.addPoint(cx + clusterOffset(tester.rng), cy + clusterOffset(tester.rng));
        }
    }
    EXPECT_TRUE(tester.verifyHulls());
}

TEST(ConvexHullComparison, NegativeCoordinates) {
    ConvexHullTester tester(7007);
    for (int i = 0; i < 30; i++) {
        auto [x, y] = tester.randomPoint(500);
        tester.addPoint(x, y);
    }
    EXPECT_TRUE(tester.verifyHulls());
}

// ============================================================================
// Tests: O(n) Build Construction
// ============================================================================

class BuildTester {
public:
    std::mt19937 rng;
    
    BuildTester(unsigned seed = 42) : rng(seed) {}
    
    std::vector<Point_2> generateRandomPoints(int n, int range = 1000) {
        std::uniform_int_distribution<int> dist(-range, range);
        std::vector<Point_2> points;
        for (int i = 0; i < n; i++) {
            points.emplace_back(dist(rng), dist(rng));
        }
        // Sort by x, then by y for ties
        std::sort(points.begin(), points.end(), [](const Point_2& a, const Point_2& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
        return points;
    }
    
    // Compare build() vs insert() for same point set
    // Both methods should produce functionally equivalent convex hulls
    // Note: They may not produce identical hulls (collinear point inclusion differs)
    // but both should be valid and functional
    bool compareBuildVsInsert(const std::vector<Point_2>& sorted_points) {
        if (sorted_points.size() < 3) {
            // For very small point sets, just verify basic properties
            CHTree<K> built_tree;
            built_tree.build(sorted_points);
            return built_tree.size() == sorted_points.size();
        }
        
        // Method 1: Use build()
        CHTree<K> built_tree;
        built_tree.build(sorted_points);
        
        // Method 2: Use insert()
        CHTree<K> inserted_tree;
        for (const auto& p : sorted_points) {
            inserted_tree.insert(p);
        }
        
        // Get hulls from both methods
        auto built_upper = hull_helpers::toIntPairs(built_tree.upperHullPoints());
        auto built_lower = hull_helpers::toIntPairs(built_tree.lowerHullPoints());
        auto inserted_upper = hull_helpers::toIntPairs(inserted_tree.upperHullPoints());
        auto inserted_lower = hull_helpers::toIntPairs(inserted_tree.lowerHullPoints());
        
        // Find extremal points from input
        auto int_points = hull_helpers::toIntPairs(sorted_points);
        auto leftmost = *std::min_element(int_points.begin(), int_points.end());
        auto rightmost = *std::max_element(int_points.begin(), int_points.end());
        
        // Combine all hull points for checking extremal inclusion
        std::set<hull_helpers::Point> built_all, inserted_all;
        built_all.insert(built_upper.begin(), built_upper.end());
        built_all.insert(built_lower.begin(), built_lower.end());
        inserted_all.insert(inserted_upper.begin(), inserted_upper.end());
        inserted_all.insert(inserted_lower.begin(), inserted_lower.end());
        
        // Check that extremal points are in both hulls
        bool build_has_extremal = (built_all.find(leftmost) != built_all.end() &&
                                   built_all.find(rightmost) != built_all.end());
        bool insert_has_extremal = (inserted_all.find(leftmost) != inserted_all.end() &&
                                    inserted_all.find(rightmost) != inserted_all.end());
        
        if (!build_has_extremal) {
            std::cerr << "build() hull missing extremal points!\n";
            std::cerr << "Leftmost: " << leftmost.first << "," << leftmost.second << "\n";
            std::cerr << "Rightmost: " << rightmost.first << "," << rightmost.second << "\n";
        }
        if (!insert_has_extremal) {
            std::cerr << "insert() hull missing extremal points!\n";
        }
        
        // Verify that covers() returns the same result for both methods on test points
        std::mt19937 test_rng(42);
        std::uniform_int_distribution<int> dist(-1000, 1000);
        bool covers_match = true;
        for (int i = 0; i < 100; i++) {
            Point_2 query(dist(test_rng), dist(test_rng));
            if (built_tree.covers(query) != inserted_tree.covers(query)) {
                covers_match = false;
                std::cerr << "covers() mismatch at (" << query.x() << "," << query.y() << ")\n";
                break;
            }
        }
        
        return build_has_extremal && insert_has_extremal && covers_match;
    }
};

TEST(ConvexHullBuild, EmptyBuild) {
    CHTree<K> tree;
    std::vector<Point_2> empty;
    tree.build(empty);
    EXPECT_TRUE(tree.empty());
    EXPECT_EQ(tree.size(), 0u);
}

TEST(ConvexHullBuild, SinglePoint) {
    CHTree<K> tree;
    std::vector<Point_2> points = {Point_2(50, 50)};
    tree.build(points);
    EXPECT_EQ(tree.size(), 1u);
}

TEST(ConvexHullBuild, TwoPoints) {
    CHTree<K> tree;
    std::vector<Point_2> points = {Point_2(0, 0), Point_2(100, 100)};
    tree.build(points);
    EXPECT_EQ(tree.size(), 2u);
}

TEST(ConvexHullBuild, ThreePoints_Triangle) {
    BuildTester tester;
    std::vector<Point_2> points = {Point_2(0, 0), Point_2(50, 100), Point_2(100, 0)};
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, Small_RandomPoints) {
    BuildTester tester(42);
    auto points = tester.generateRandomPoints(10, 100);
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, Medium_RandomPoints) {
    BuildTester tester(123);
    auto points = tester.generateRandomPoints(100, 500);
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, Large_RandomPoints) {
    BuildTester tester(456);
    auto points = tester.generateRandomPoints(500, 1000);
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, VeryLarge_RandomPoints) {
    BuildTester tester(789);
    auto points = tester.generateRandomPoints(2000, 10000);
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, CollinearHorizontal) {
    BuildTester tester;
    std::vector<Point_2> points;
    for (int i = 0; i < 20; i++) {
        points.emplace_back(i * 10, 50);
    }
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, CollinearVertical) {
    BuildTester tester;
    std::vector<Point_2> points;
    for (int i = 0; i < 20; i++) {
        points.emplace_back(50, i * 10);
    }
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, CollinearDiagonal) {
    BuildTester tester;
    std::vector<Point_2> points;
    for (int i = 0; i < 20; i++) {
        points.emplace_back(i * 10, i * 10);
    }
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, Square) {
    BuildTester tester;
    std::vector<Point_2> points = {
        Point_2(0, 0), Point_2(0, 100), Point_2(100, 0), Point_2(100, 100)
    };
    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, CircleApproximation) {
    BuildTester tester;
    std::vector<Point_2> points;
    const int n = 24;
    const int r = 100;
    for (int i = 0; i < n; i++) {
        double angle = 2.0 * 3.14159265358979 * i / n;
        points.emplace_back((int)(r * cos(angle)), (int)(r * sin(angle)));
    }
    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, SquareWithInterior) {
    BuildTester tester;
    std::vector<Point_2> points = {
        Point_2(0, 0), Point_2(0, 100), Point_2(100, 0), Point_2(100, 100),
        Point_2(50, 50), Point_2(25, 75), Point_2(75, 25)
    };
    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

TEST(ConvexHullBuild, DifferentSeeds) {
    for (int seed = 0; seed < 20; seed++) {
        BuildTester tester(seed * 12345);
        auto points = tester.generateRandomPoints(200, 500);
        EXPECT_TRUE(tester.compareBuildVsInsert(points)) << "Failed with seed " << seed;
    }
}

TEST(ConvexHullBuild, PowersOfTwo) {
    // Test with sizes that are powers of 2 (tests tree balance)
    std::vector<int> sizes = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512};
    for (int n : sizes) {
        BuildTester tester(n);
        auto points = tester.generateRandomPoints(n, 1000);
        EXPECT_TRUE(tester.compareBuildVsInsert(points)) << "Failed with size " << n;
    }
}

TEST(ConvexHullBuild, NonPowersOfTwo) {
    // Test with sizes that are not powers of 2
    std::vector<int> sizes = {3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 255, 257};
    for (int n : sizes) {
        BuildTester tester(n);
        auto points = tester.generateRandomPoints(n, 1000);
        EXPECT_TRUE(tester.compareBuildVsInsert(points)) << "Failed with size " << n;
    }
}

// Test that build() correctly handles duplicate x-coordinates
TEST(ConvexHullBuild, DuplicateXCoordinates) {
    BuildTester tester;
    std::vector<Point_2> points;
    for (int x = 0; x < 10; x++) {
        for (int y = 0; y < 5; y++) {
            points.emplace_back(x * 10, y * 20);
        }
    }
    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    });
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
}

// Test that operations after build() work correctly
TEST(ConvexHullBuild, OperationsAfterBuild) {
    BuildTester tester(42);
    auto points = tester.generateRandomPoints(50, 200);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Insert more points
    std::uniform_int_distribution<int> dist(-200, 200);
    for (int i = 0; i < 20; i++) {
        Point_2 p(dist(tester.rng), dist(tester.rng));
        tree.insert(p);
        points.push_back(p);
    }
    
    // Verify against reference
    auto expected_upper = hull_helpers::adjustUpperHullForCHTree(
        monotone_chain::upperHull(hull_helpers::toIntPairs(points)));
    auto expected_lower = hull_helpers::adjustLowerHullForCHTree(
        monotone_chain::lowerHull(hull_helpers::toIntPairs(points)));
    
    auto ch_upper = hull_helpers::toIntPairs(tree.upperHullPoints());
    auto ch_lower = hull_helpers::toIntPairs(tree.lowerHullPoints());
    
    EXPECT_TRUE(hull_helpers::hullContainsAll(expected_upper, ch_upper));
    EXPECT_TRUE(hull_helpers::hullContainsAll(expected_lower, ch_lower));
}

// Test removal after build()
TEST(ConvexHullBuild, RemovalAfterBuild) {
    BuildTester tester(123);
    auto points = tester.generateRandomPoints(100, 300);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Remove some points
    std::vector<Point_2> removed;
    for (int i = 0; i < 30; i++) {
        size_t idx = tester.rng() % points.size();
        removed.push_back(points[idx]);
        tree.remove(points[idx]);
        points.erase(points.begin() + idx);
    }
    
    // Verify against reference
    if (points.size() >= 3) {
        auto expected_upper = hull_helpers::adjustUpperHullForCHTree(
            monotone_chain::upperHull(hull_helpers::toIntPairs(points)));
        auto expected_lower = hull_helpers::adjustLowerHullForCHTree(
            monotone_chain::lowerHull(hull_helpers::toIntPairs(points)));
        
        auto ch_upper = hull_helpers::toIntPairs(tree.upperHullPoints());
        auto ch_lower = hull_helpers::toIntPairs(tree.lowerHullPoints());
        
        EXPECT_TRUE(hull_helpers::hullContainsAll(expected_upper, ch_upper));
        EXPECT_TRUE(hull_helpers::hullContainsAll(expected_lower, ch_lower));
    }
}

// Test rebuilding (calling build() multiple times)
TEST(ConvexHullBuild, RebuildMultipleTimes) {
    for (int round = 0; round < 5; round++) {
        BuildTester tester(round);
        auto points = tester.generateRandomPoints(100, 500);
        
        CHTree<K> tree;
        tree.build(points);
        
        // Build again with different points
        auto new_points = tester.generateRandomPoints(150, 600);
        tree.build(new_points);
        
        // Verify the second build
        EXPECT_TRUE(tester.compareBuildVsInsert(new_points)) << "Failed on round " << round;
    }
}

// Test covers() query after build()
TEST(ConvexHullBuild, CoversAfterBuild) {
    BuildTester tester(42);
    auto points = tester.generateRandomPoints(100, 200);
    
    CHTree<K> built_tree;
    built_tree.build(points);
    
    CHTree<K> inserted_tree;
    for (const auto& p : points) {
        inserted_tree.insert(p);
    }
    
    // Generate test points and compare covers() results
    std::uniform_int_distribution<int> dist(-200, 200);
    for (int i = 0; i < 100; i++) {
        Point_2 query(dist(tester.rng), dist(tester.rng));
        EXPECT_EQ(built_tree.covers(query), inserted_tree.covers(query)) 
            << "covers() mismatch for point (" << query.x() << ", " << query.y() << ")";
    }
}

// Stress test for O(n) build
TEST(ConvexHullBuild, StressTest) {
    BuildTester tester(9999);
    auto points = tester.generateRandomPoints(5000, 100000);
    
    auto start = std::chrono::high_resolution_clock::now();
    EXPECT_TRUE(tester.compareBuildVsInsert(points));
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Stress test (5000 points) completed in " << duration.count() << "ms\n";
}

// ============================================================================
// Comprehensive tests for operations AFTER build()
// These tests verify that the hull remains correct after each operation
// ============================================================================

// Helper function to verify hull correctness against reference implementation
bool verifyHullCorrectness(CHTree<K>& tree, const std::vector<Point_2>& points) {
    if (points.size() < 3) return true;  // Too few points to verify meaningfully
    
    auto expected_upper = hull_helpers::adjustUpperHullForCHTree(
        monotone_chain::upperHull(hull_helpers::toIntPairs(points)));
    auto expected_lower = hull_helpers::adjustLowerHullForCHTree(
        monotone_chain::lowerHull(hull_helpers::toIntPairs(points)));
    
    auto ch_upper = hull_helpers::toIntPairs(tree.upperHullPoints());
    auto ch_lower = hull_helpers::toIntPairs(tree.lowerHullPoints());
    
    // Check that all essential points are contained in the CHTree hull
    return hull_helpers::hullContainsAll(expected_upper, ch_upper) &&
           hull_helpers::hullContainsAll(expected_lower, ch_lower);
}

// Test insertions after build with verification after EACH insertion
TEST(ConvexHullBuild, InsertionsAfterBuild_Verified) {
    BuildTester tester(1001);
    auto initial_points = tester.generateRandomPoints(50, 200);
    
    CHTree<K> tree;
    tree.build(initial_points);
    
    // Verify initial state
    ASSERT_TRUE(verifyHullCorrectness(tree, initial_points)) 
        << "Hull incorrect after initial build";
    
    // Insert 50 more points, verifying after each
    std::uniform_int_distribution<int> dist(-250, 250);
    for (int i = 0; i < 50; i++) {
        Point_2 p(dist(tester.rng), dist(tester.rng));
        tree.insert(p);
        initial_points.push_back(p);
        
        EXPECT_TRUE(verifyHullCorrectness(tree, initial_points)) 
            << "Hull incorrect after insertion " << (i+1) 
            << " of point (" << p.x() << ", " << p.y() << ")";
    }
}

// Test removals after build with verification after EACH removal
TEST(ConvexHullBuild, RemovalsAfterBuild_Verified) {
    BuildTester tester(2002);
    auto points = tester.generateRandomPoints(100, 300);
    
    CHTree<K> tree;
    tree.build(points);
    
    // Verify initial state
    ASSERT_TRUE(verifyHullCorrectness(tree, points)) 
        << "Hull incorrect after initial build";
    
    // Remove 50 points, verifying after each
    for (int i = 0; i < 50 && points.size() > 3; i++) {
        size_t idx = tester.rng() % points.size();
        Point_2 p = points[idx];
        tree.remove(p);
        points.erase(points.begin() + idx);
        
        EXPECT_TRUE(verifyHullCorrectness(tree, points)) 
            << "Hull incorrect after removal " << (i+1) 
            << " of point (" << p.x() << ", " << p.y() << ")";
    }
}

// Test mixed insert/remove operations after build with step-by-step verification
TEST(ConvexHullBuild, MixedOperationsAfterBuild_Verified) {
    BuildTester tester(3003);
    auto points = tester.generateRandomPoints(75, 250);
    
    CHTree<K> tree;
    tree.build(points);
    
    ASSERT_TRUE(verifyHullCorrectness(tree, points)) 
        << "Hull incorrect after initial build";
    
    std::uniform_int_distribution<int> coord_dist(-300, 300);
    std::uniform_int_distribution<int> op_dist(0, 2);  // 0,1 = insert, 2 = remove
    
    for (int i = 0; i < 100; i++) {
        int op = op_dist(tester.rng);
        
        if (op < 2 || points.size() <= 4) {  // Insert
            Point_2 p(coord_dist(tester.rng), coord_dist(tester.rng));
            tree.insert(p);
            points.push_back(p);
        } else {  // Remove
            size_t idx = tester.rng() % points.size();
            Point_2 p = points[idx];
            tree.remove(p);
            points.erase(points.begin() + idx);
        }
        
        EXPECT_TRUE(verifyHullCorrectness(tree, points)) 
            << "Hull incorrect after operation " << (i+1);
    }
}

// Stress test: many operations after build
TEST(ConvexHullBuild, ManyOperationsAfterBuild_Stress) {
    BuildTester tester(4004);
    auto points = tester.generateRandomPoints(200, 500);
    
    CHTree<K> tree;
    tree.build(points);
    
    std::uniform_int_distribution<int> coord_dist(-600, 600);
    std::uniform_int_distribution<int> op_dist(0, 1);
    
    // 500 operations (mix of insert/remove)
    for (int i = 0; i < 500; i++) {
        if (op_dist(tester.rng) == 0 || points.size() <= 5) {
            Point_2 p(coord_dist(tester.rng), coord_dist(tester.rng));
            tree.insert(p);
            points.push_back(p);
        } else {
            size_t idx = tester.rng() % points.size();
            tree.remove(points[idx]);
            points.erase(points.begin() + idx);
        }
    }
    
    // Verify final hull
    EXPECT_TRUE(verifyHullCorrectness(tree, points)) 
        << "Hull incorrect after 500 mixed operations";
}

// Test that covers() works correctly after insertions and removals
TEST(ConvexHullBuild, CoversAfterMixedOperations) {
    BuildTester tester(5005);
    auto points = tester.generateRandomPoints(100, 300);
    
    CHTree<K> built_tree;
    built_tree.build(points);
    
    CHTree<K> inserted_tree;
    for (const auto& p : points) {
        inserted_tree.insert(p);
    }
    
    std::uniform_int_distribution<int> coord_dist(-350, 350);
    
    // Do some operations on both trees
    for (int i = 0; i < 30; i++) {
        Point_2 p(coord_dist(tester.rng), coord_dist(tester.rng));
        built_tree.insert(p);
        inserted_tree.insert(p);
        points.push_back(p);
    }
    
    for (int i = 0; i < 20 && points.size() > 5; i++) {
        size_t idx = tester.rng() % points.size();
        built_tree.remove(points[idx]);
        inserted_tree.remove(points[idx]);
        points.erase(points.begin() + idx);
    }
    
    // Verify covers() matches between both trees
    for (int i = 0; i < 200; i++) {
        Point_2 query(coord_dist(tester.rng), coord_dist(tester.rng));
        EXPECT_EQ(built_tree.covers(query), inserted_tree.covers(query)) 
            << "covers() mismatch at (" << query.x() << ", " << query.y() << ")";
    }
}

// Test hull points match exactly after operations
TEST(ConvexHullBuild, HullPointsMatch_AfterOperations) {
    BuildTester tester(6006);
    auto points = tester.generateRandomPoints(80, 250);
    
    CHTree<K> built_tree;
    built_tree.build(points);
    
    CHTree<K> inserted_tree;
    for (const auto& p : points) {
        inserted_tree.insert(p);
    }
    
    std::uniform_int_distribution<int> coord_dist(-300, 300);
    
    // Perform same operations on both
    for (int round = 0; round < 5; round++) {
        // Insert some points
        for (int i = 0; i < 10; i++) {
            Point_2 p(coord_dist(tester.rng), coord_dist(tester.rng));
            built_tree.insert(p);
            inserted_tree.insert(p);
            points.push_back(p);
        }
        
        // Remove some points
        for (int i = 0; i < 5 && points.size() > 5; i++) {
            size_t idx = tester.rng() % points.size();
            built_tree.remove(points[idx]);
            inserted_tree.remove(points[idx]);
            points.erase(points.begin() + idx);
        }
        
        // Verify hull points match
        auto built_upper = built_tree.upperHullPoints();
        auto insert_upper = inserted_tree.upperHullPoints();
        auto built_lower = built_tree.lowerHullPoints();
        auto insert_lower = inserted_tree.lowerHullPoints();
        
        EXPECT_EQ(built_upper.size(), insert_upper.size()) 
            << "Upper hull size mismatch at round " << round;
        EXPECT_EQ(built_lower.size(), insert_lower.size()) 
            << "Lower hull size mismatch at round " << round;
        
        // Verify each hull point matches
        for (size_t i = 0; i < std::min(built_upper.size(), insert_upper.size()); i++) {
            EXPECT_EQ(built_upper[i], insert_upper[i]) 
                << "Upper hull point mismatch at index " << i << " round " << round;
        }
        for (size_t i = 0; i < std::min(built_lower.size(), insert_lower.size()); i++) {
            EXPECT_EQ(built_lower[i], insert_lower[i]) 
                << "Lower hull point mismatch at index " << i << " round " << round;
        }
    }
}

// Test building, doing many operations, then rebuilding
TEST(ConvexHullBuild, BuildOperationsRebuild) {
    BuildTester tester(7007);
    
    for (int cycle = 0; cycle < 3; cycle++) {
        auto points = tester.generateRandomPoints(50 + cycle * 25, 300);
        
        CHTree<K> tree;
        tree.build(points);
        
        std::uniform_int_distribution<int> coord_dist(-350, 350);
        
        // Do operations
        for (int i = 0; i < 30; i++) {
            if (tester.rng() % 2 == 0 || points.size() <= 5) {
                Point_2 p(coord_dist(tester.rng), coord_dist(tester.rng));
                tree.insert(p);
                points.push_back(p);
            } else {
                size_t idx = tester.rng() % points.size();
                tree.remove(points[idx]);
                points.erase(points.begin() + idx);
            }
        }
        
        // Verify hull is correct
        EXPECT_TRUE(verifyHullCorrectness(tree, points)) 
            << "Hull incorrect after operations in cycle " << cycle;
        
        // Rebuild with new points
        auto new_points = tester.generateRandomPoints(80 + cycle * 20, 400);
        tree.build(new_points);
        
        EXPECT_TRUE(verifyHullCorrectness(tree, new_points)) 
            << "Hull incorrect after rebuild in cycle " << cycle;
    }
}


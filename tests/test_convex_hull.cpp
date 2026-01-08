// Tests comparing Dynamic Convex Hull (CHTree) to a standard convex hull solution
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <memory>
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

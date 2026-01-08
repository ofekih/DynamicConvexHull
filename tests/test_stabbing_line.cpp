// Tests for StabbingLineStructure - epsilon-stabbing line query
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include "StabbingLineStructure.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef StabbingLineStructure<K> SLS;
typedef StabbingLine<K> Line;

// ============================================================================
// Brute-force reference implementation for verification
// ============================================================================

// Check if a line with given slope and intercept covers all points within epsilon
bool bruteForceCheck(const std::vector<Point_2>& points, double slope, double intercept, double epsilon) {
    for (const auto& p : points) {
        double lineY = slope * p.x() + intercept;
        if (std::abs(p.y() - lineY) > epsilon + 1e-9) {
            return false;
        }
    }
    return true;
}

// Try to find any stabbing line using brute force (for small point sets)
// Returns true if a stabbing line exists
bool bruteForceFindStabbingLine(const std::vector<Point_2>& points, double epsilon) {
    if (points.empty()) return true;
    if (points.size() == 1) return true;
    
    // For each pair of points, try the line through their epsilon boundaries
    // This is O(n² * n) but fine for testing small sets
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i; j < points.size(); ++j) {
            // Try 4 combinations of upper/lower epsilon bounds
            std::vector<Point_2> corners = {
                Point_2(points[i].x(), points[i].y() + epsilon),
                Point_2(points[i].x(), points[i].y() - epsilon),
                Point_2(points[j].x(), points[j].y() + epsilon),
                Point_2(points[j].x(), points[j].y() - epsilon)
            };
            
            for (size_t a = 0; a < 2; ++a) {
                for (size_t b = 2; b < 4; ++b) {
                    Point_2 p1 = corners[a];
                    Point_2 p2 = corners[b];
                    
                    double dx = p2.x() - p1.x();
                    double slope, intercept;
                    if (std::abs(dx) < 1e-12) {
                        // Vertical - check horizontal line through midpoint
                        slope = 0;
                        intercept = (p1.y() + p2.y()) / 2.0;
                    } else {
                        slope = (p2.y() - p1.y()) / dx;
                        intercept = p1.y() - slope * p1.x();
                    }
                    
                    if (bruteForceCheck(points, slope, intercept, epsilon)) {
                        return true;
                    }
                }
            }
        }
    }
    
    // Also try horizontal lines at various y-values
    for (const auto& p : points) {
        if (bruteForceCheck(points, 0, p.y(), epsilon)) {
            return true;
        }
    }
    
    return false;
}

// ============================================================================
// Basic Edge Cases
// ============================================================================

TEST(StabbingLine, EmptySet) {
    SLS sls(1.0);
    EXPECT_TRUE(sls.hasStabbingLine());
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
}

TEST(StabbingLine, SinglePoint) {
    SLS sls(1.0);
    sls.insert(Point_2(5, 10));
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    
    // The line should pass within epsilon of the point
    double lineY = line->at(5);
    EXPECT_LE(std::abs(lineY - 10), 1.0 + 1e-9);
}

TEST(StabbingLine, TwoPoints_Horizontal) {
    SLS sls(1.0);
    sls.insert(Point_2(0, 10));
    sls.insert(Point_2(10, 10));
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    EXPECT_TRUE(line->coversPoint(Point_2(0, 10), 1.0));
    EXPECT_TRUE(line->coversPoint(Point_2(10, 10), 1.0));
}

TEST(StabbingLine, TwoPoints_Diagonal) {
    SLS sls(1.0);
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(10, 10));
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    EXPECT_TRUE(line->coversPoint(Point_2(0, 0), 1.0));
    EXPECT_TRUE(line->coversPoint(Point_2(10, 10), 1.0));
}

TEST(StabbingLine, TwoPoints_TooFarApart) {
    // Two points at the SAME x-coordinate, far apart vertically
    // No non-vertical line can pass within epsilon of both
    SLS sls(0.1);
    sls.insert(Point_2(5, 0));
    sls.insert(Point_2(5, 10));  // Same x, far apart in y
    EXPECT_FALSE(sls.hasStabbingLine());
}

// ============================================================================
// Collinear Points
// ============================================================================

TEST(StabbingLine, CollinearHorizontal) {
    SLS sls(1.0);
    for (int x = 0; x < 10; ++x) {
        sls.insert(Point_2(x, 5));
    }
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    // Should be approximately y = 5
    EXPECT_NEAR(line->at(5), 5.0, 1.0);
}

TEST(StabbingLine, CollinearDiagonal) {
    SLS sls(1.0);
    for (int i = 0; i < 10; ++i) {
        sls.insert(Point_2(i, i * 2));  // y = 2x
    }
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    // All points should be covered
    for (int i = 0; i < 10; ++i) {
        EXPECT_TRUE(line->coversPoint(Point_2(i, i * 2), 1.0));
    }
}

// ============================================================================
// Valid Stabbing Lines (Should Exist)
// ============================================================================

TEST(StabbingLine, PointsNearLine) {
    SLS sls(1.0);
    // Points close to y = x
    sls.insert(Point_2(0, 0.5));
    sls.insert(Point_2(5, 5.3));
    sls.insert(Point_2(10, 9.8));
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    EXPECT_TRUE(line->coversPoint(Point_2(0, 0.5), 1.0));
    EXPECT_TRUE(line->coversPoint(Point_2(5, 5.3), 1.0));
    EXPECT_TRUE(line->coversPoint(Point_2(10, 9.8), 1.0));
}

TEST(StabbingLine, LargeEpsilon) {
    SLS sls(100.0);
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(10, 50));
    sls.insert(Point_2(20, -30));
    EXPECT_TRUE(sls.hasStabbingLine());
}

// ============================================================================
// Invalid Stabbing Lines (Should Not Exist)
// ============================================================================

TEST(StabbingLine, OutlierPoint) {
    SLS sls(1.0);
    // Three collinear points and one outlier
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(5, 5));
    sls.insert(Point_2(10, 10));
    sls.insert(Point_2(5, 20));  // Outlier - far from y=x
    EXPECT_FALSE(sls.hasStabbingLine());
}

TEST(StabbingLine, ZigZagPattern) {
    SLS sls(0.5);
    // Points that zig-zag too much for epsilon=0.5
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(1, 5));
    sls.insert(Point_2(2, 0));
    sls.insert(Point_2(3, 5));
    // No line can fit through this with epsilon=0.5
    EXPECT_FALSE(sls.hasStabbingLine());
}

// ============================================================================
// Dynamic Operations
// ============================================================================

TEST(StabbingLine, InsertMakesInvalid) {
    SLS sls(1.0);
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(10, 10));
    EXPECT_TRUE(sls.hasStabbingLine());
    
    // Add an outlier
    sls.insert(Point_2(5, 100));
    EXPECT_FALSE(sls.hasStabbingLine());
}

TEST(StabbingLine, RemoveMakesValid) {
    SLS sls(1.0);
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(10, 10));
    sls.insert(Point_2(5, 100));  // Outlier
    EXPECT_FALSE(sls.hasStabbingLine());
    
    // Remove the outlier
    sls.remove(Point_2(5, 100));
    EXPECT_TRUE(sls.hasStabbingLine());
}

TEST(StabbingLine, BuildFromSorted) {
    SLS sls(1.0);
    std::vector<Point_2> points = {
        Point_2(0, 0),
        Point_2(2, 2.1),
        Point_2(5, 4.9),
        Point_2(8, 8),
        Point_2(10, 10.2)
    };
    sls.build(points);
    EXPECT_TRUE(sls.hasStabbingLine());
    EXPECT_EQ(sls.size(), 5u);
}

// ============================================================================
// Split and Join
// ============================================================================

TEST(StabbingLine, SplitBasic) {
    SLS sls(1.0);
    for (int x = 0; x < 10; ++x) {
        sls.insert(Point_2(x, x));
    }
    EXPECT_TRUE(sls.hasStabbingLine());
    
    auto right = sls.split(5);
    
    // Both parts should have stabbing lines
    EXPECT_TRUE(sls.hasStabbingLine());
    EXPECT_TRUE(right.hasStabbingLine());
}

TEST(StabbingLine, JoinBasic) {
    SLS left(1.0);
    SLS right(1.0);
    
    for (int x = 0; x < 5; ++x) left.insert(Point_2(x, x));
    for (int x = 5; x < 10; ++x) right.insert(Point_2(x, x));
    
    left.join(right);
    
    EXPECT_TRUE(left.hasStabbingLine());
    EXPECT_EQ(right.size(), 0u);
}

TEST(StabbingLine, JoinMakesInvalid) {
    SLS left(1.0);
    SLS right(1.0);
    
    // Two sets that individually have stabbing lines
    // but together don't
    for (int x = 0; x < 5; ++x) left.insert(Point_2(x, 0));
    for (int x = 10; x < 15; ++x) right.insert(Point_2(x, 100));
    
    EXPECT_TRUE(left.hasStabbingLine());
    EXPECT_TRUE(right.hasStabbingLine());
    
    left.join(right);
    EXPECT_FALSE(left.hasStabbingLine());
}

// ============================================================================
// Randomized Tests with Brute-Force Verification
// ============================================================================

TEST(StabbingLine, RandomSmall_Valid) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos(-10, 10);
    std::uniform_real_distribution<double> noise(-0.5, 0.5);
    
    for (int trial = 0; trial < 20; ++trial) {
        SLS sls(1.0);
        std::vector<Point_2> points;
        
        // Generate points close to y = mx + b
        double slope = pos(rng) / 5.0;
        double intercept = pos(rng);
        
        for (int i = 0; i < 10; ++i) {
            double x = i * 2.0;
            double y = slope * x + intercept + noise(rng);
            points.emplace_back(x, y);
            sls.insert(Point_2(x, y));
        }
        
        // Should have a valid stabbing line
        bool expected = bruteForceFindStabbingLine(points, 1.0);
        bool actual = sls.hasStabbingLine();
        
        EXPECT_EQ(expected, actual) 
            << "Trial " << trial << ": expected=" << expected << ", actual=" << actual;
    }
}

TEST(StabbingLine, RandomSmall_Invalid) {
    std::mt19937 rng(123);
    std::uniform_real_distribution<double> pos(-100, 100);
    
    for (int trial = 0; trial < 20; ++trial) {
        SLS sls(0.5);  // Small epsilon
        std::vector<Point_2> points;
        
        // Generate random scattered points
        for (int i = 0; i < 8; ++i) {
            double x = i * 5.0;
            double y = pos(rng);
            points.emplace_back(x, y);
            sls.insert(Point_2(x, y));
        }
        
        // Compare with brute force
        bool expected = bruteForceFindStabbingLine(points, 0.5);
        bool actual = sls.hasStabbingLine();
        
        EXPECT_EQ(expected, actual) 
            << "Trial " << trial << ": expected=" << expected << ", actual=" << actual;
    }
}

TEST(StabbingLine, RandomMedium) {
    std::mt19937 rng(456);
    std::uniform_real_distribution<double> noise(-0.3, 0.3);
    
    for (int trial = 0; trial < 10; ++trial) {
        SLS sls(0.5);
        std::vector<Point_2> points;
        
        // Generate points with small noise around y = x
        for (int i = 0; i < 50; ++i) {
            double x = i;
            double y = x + noise(rng);
            points.emplace_back(x, y);
        }
        
        // Sort and build
        std::sort(points.begin(), points.end(), 
            [](const Point_2& a, const Point_2& b) {
                if (a.x() != b.x()) return a.x() < b.x();
                return a.y() < b.y();
            });
        sls.build(points);
        
        // Should have a stabbing line
        EXPECT_TRUE(sls.hasStabbingLine()) << "Trial " << trial;
        
        // Verify the returned line
        auto line = sls.findStabbingLine();
        ASSERT_TRUE(line.has_value());
        for (const auto& p : points) {
            EXPECT_TRUE(line->coversPoint(p, 0.5)) 
                << "Line doesn't cover point (" << p.x() << "," << p.y() << ")";
        }
    }
}

// ============================================================================
// Stress Test
// ============================================================================

TEST(StabbingLine, StressTest_Large) {
    std::mt19937 rng(789);
    std::uniform_real_distribution<double> noise(-0.1, 0.1);
    
    SLS sls(0.5);
    std::vector<Point_2> points;
    
    // Generate 1000 points close to y = 2x + 3
    for (int i = 0; i < 1000; ++i) {
        double x = i * 0.1;
        double y = 2.0 * x + 3.0 + noise(rng);
        points.emplace_back(x, y);
    }
    
    // Sort and build
    std::sort(points.begin(), points.end(), 
        [](const Point_2& a, const Point_2& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
    sls.build(points);
    
    EXPECT_TRUE(sls.hasStabbingLine());
    EXPECT_EQ(sls.size(), 1000u);
    
    // Verify the returned line covers all points
    auto line = sls.findStabbingLine();
    ASSERT_TRUE(line.has_value());
    for (const auto& p : points) {
        EXPECT_TRUE(line->coversPoint(p, 0.5)) 
            << "Line doesn't cover point (" << p.x() << "," << p.y() << ")";
    }
}

TEST(StabbingLine, StressTest_DynamicOperations) {
    std::mt19937 rng(1001);
    std::uniform_real_distribution<double> noise(-0.2, 0.2);
    std::uniform_int_distribution<int> action(0, 4);
    
    SLS sls(0.5);
    std::vector<Point_2> points;
    
    // Start with some points near y = x
    for (int i = 0; i < 20; ++i) {
        double x = i;
        double y = x + noise(rng);
        points.emplace_back(x, y);
        sls.insert(points.back());
    }
    
    // Perform random insertions and removals
    for (int iter = 0; iter < 100; ++iter) {
        if (points.empty() || action(rng) != 0) {
            // Insert
            double x = 50.0 + iter * 0.1;
            double y = x + noise(rng);
            points.emplace_back(x, y);
            sls.insert(points.back());
        } else {
            // Remove
            size_t idx = rng() % points.size();
            sls.remove(points[idx]);
            points.erase(points.begin() + idx);
        }
    }
    
    // Should still have a stabbing line (all points near y=x)
    EXPECT_TRUE(sls.hasStabbingLine());
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// Tests for StabbingLineStructure - epsilon-stabbing line query
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>
#include <sstream>
#include "StabbingLineStructure.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef StabbingLineStructure<K> SLS;
typedef StabbingLine<K> Line;

// ============================================================================
// Comprehensive Verification Framework
// ============================================================================

// O(n) verification: Check if a line covers ALL points within epsilon
// This is the definitive correctness check for any returned line
bool verifyLineCoversAllPoints(const std::vector<Point_2>& points, 
                                double slope, double intercept, 
                                double epsilon) {
    for (const auto& p : points) {
        double lineY = slope * p.x() + intercept;
        double error = std::abs(p.y() - lineY);
        if (error > epsilon + 1e-9) {
            return false;
        }
    }
    return true;
}

// Wrapper for Line object
bool verifyLineCoversAllPoints(const std::vector<Point_2>& points,
                                const Line& line,
                                double epsilon) {
    return verifyLineCoversAllPoints(points, line.slope, line.intercept, epsilon);
}

// ============================================================================
// Independent Reference Implementation using Half-Plane Intersection
// ============================================================================
//
// For each point (xi, yi), the constraint |yi - (slope * xi + intercept)| <= epsilon
// gives two linear constraints on (slope, intercept):
//   slope * xi + intercept >= yi - epsilon   (lower bound)
//   slope * xi + intercept <= yi + epsilon   (upper bound)
//
// Rewriting in terms of intercept = b and slope = m:
//   b >= (yi - epsilon) - m * xi   [lower constraint]
//   b <= (yi + epsilon) - m * xi   [upper constraint]
//
// A valid line exists iff the feasible region (intersection of half-planes) is non-empty.
//
// For a given slope m, the feasible range of intercepts is:
//   b_min(m) = max_i { yi - epsilon - m * xi }
//   b_max(m) = min_i { yi + epsilon - m * xi }
// A valid line with slope m exists iff b_min(m) <= b_max(m)
//
// The overall problem is feasible iff there exists m such that b_min(m) <= b_max(m)
//
// b_min(m) is a piecewise-linear convex function (max of lines)
// b_max(m) is a piecewise-linear concave function (min of lines)
// They intersect where the problem becomes infeasible.

struct HalfPlaneIntersection {
    std::vector<Point_2> points;
    double epsilon;
    
    HalfPlaneIntersection(const std::vector<Point_2>& pts, double eps) 
        : points(pts), epsilon(eps) {}
    
    // Get b_min for a given slope m
    double getBMin(double m) const {
        double bMin = -std::numeric_limits<double>::infinity();
        for (const auto& p : points) {
            double constraint = (p.y() - epsilon) - m * p.x();
            bMin = std::max(bMin, constraint);
        }
        return bMin;
    }
    
    // Get b_max for a given slope m
    double getBMax(double m) const {
        double bMax = std::numeric_limits<double>::infinity();
        for (const auto& p : points) {
            double constraint = (p.y() + epsilon) - m * p.x();
            bMax = std::min(bMax, constraint);
        }
        return bMax;
    }
    
    // Check if slope m is feasible
    bool isFeasible(double m) const {
        return getBMin(m) <= getBMax(m) + 1e-9;
    }
    
    // Find a valid line using ternary search on slope
    // Returns {hasLine, slope, intercept}
    std::tuple<bool, double, double> findLine() const {
        if (points.empty()) {
            return {true, 0, 0};
        }
        if (points.size() == 1) {
            return {true, 0, points[0].y()};
        }
        
        // Get slope bounds from points
        double minSlope = -1e9, maxSlope = 1e9;
        
        // Binary search to find if there's any feasible slope
        // First, sample a range of slopes to find one that works
        std::vector<double> slopesToTry;
        
        // Add slopes defined by pairs of points' boundaries
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                double dx = points[j].x() - points[i].x();
                if (std::abs(dx) < 1e-12) continue;
                
                // Slopes connecting boundary corners
                double dy1 = (points[j].y() - epsilon) - (points[i].y() + epsilon);
                double dy2 = (points[j].y() + epsilon) - (points[i].y() - epsilon);
                double dy3 = (points[j].y() - epsilon) - (points[i].y() - epsilon);
                double dy4 = (points[j].y() + epsilon) - (points[i].y() + epsilon);
                
                slopesToTry.push_back(dy1 / dx);
                slopesToTry.push_back(dy2 / dx);
                slopesToTry.push_back(dy3 / dx);
                slopesToTry.push_back(dy4 / dx);
            }
        }
        
        // Also try horizontal line
        slopesToTry.push_back(0.0);
        
        // Try each candidate slope
        for (double m : slopesToTry) {
            if (isFeasible(m)) {
                double bMin = getBMin(m);
                double bMax = getBMax(m);
                double b = (bMin + bMax) / 2.0;
                return {true, m, b};
            }
        }
        
        // If no candidate worked, do a denser grid search
        for (double m = -1000; m <= 1000; m += 0.1) {
            if (isFeasible(m)) {
                double bMin = getBMin(m);
                double bMax = getBMax(m);
                double b = (bMin + bMax) / 2.0;
                return {true, m, b};
            }
        }
        
        return {false, 0, 0};
    }
    
    bool hasStabbingLine() const {
        return std::get<0>(findLine());
    }
};

// ============================================================================
// Simple O(n² * n) brute force for small sets
// ============================================================================

bool bruteForceCheck(const std::vector<Point_2>& points, double slope, double intercept, double epsilon) {
    return verifyLineCoversAllPoints(points, slope, intercept, epsilon);
}

bool bruteForceFindStabbingLine(const std::vector<Point_2>& points, double epsilon) {
    if (points.empty()) return true;
    if (points.size() == 1) return true;
    
    // Use the more sophisticated half-plane intersection method
    HalfPlaneIntersection hpi(points, epsilon);
    return hpi.hasStabbingLine();
}

// ============================================================================
// Verification helper: check that algorithm agrees with reference
// ============================================================================

struct VerificationResult {
    bool algorithmResult;
    bool referenceResult;
    bool algorithmLineValid;  // If algorithm returns line, is it actually valid?
    std::string details;
    
    bool passed() const {
        if (algorithmResult != referenceResult) return false;
        if (algorithmResult && !algorithmLineValid) return false;
        return true;
    }
};

VerificationResult verifyStabbingLineQuery(const std::vector<Point_2>& points, 
                                            double epsilon,
                                            bool debug = false) {
    VerificationResult result;
    
    // Build the structure
    SLS sls(epsilon);
    std::vector<Point_2> sortedPoints = points;
    std::sort(sortedPoints.begin(), sortedPoints.end(), 
        [](const Point_2& a, const Point_2& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
    if (!sortedPoints.empty()) {
        sls.build(sortedPoints);
    }
    
    if (debug) {
        std::cout << "DEBUG: Input points (" << points.size() << "): ";
        for (const auto& p : points) std::cout << "(" << p.x() << "," << p.y() << ") ";
        std::cout << "\n";
        std::cout << "DEBUG: Sorted points (" << sortedPoints.size() << "): ";
        for (const auto& p : sortedPoints) std::cout << "(" << p.x() << "," << p.y() << ") ";
        std::cout << "\n";
        std::cout << "DEBUG: SLS points (" << sls.getOriginalPoints().size() << "): ";
        for (const auto& p : sls.getOriginalPoints()) std::cout << "(" << p.x() << "," << p.y() << ") ";
        std::cout << "\n";
        std::cout << "DEBUG: epsilon = " << epsilon << "\n";
    }
    
    // Get algorithm result
    result.algorithmResult = sls.hasStabbingLine();
    
    // Get reference result
    HalfPlaneIntersection hpi(points, epsilon);
    result.referenceResult = hpi.hasStabbingLine();
    
    if (debug && result.algorithmResult != result.referenceResult) {
        // Try to find what slope the reference finds
        auto [ok, m, b] = hpi.findLine();
        if (ok) {
            std::cout << "DEBUG: Reference found line: y = " << m << "x + " << b << "\n";
            // Check if algorithm's checkSlope would find this valid
            double bMin = -1e308, bMax = 1e308;
            for (const auto& p : sls.getOriginalPoints()) {
                bMin = std::max(bMin, p.y() - epsilon - m * p.x());
                bMax = std::min(bMax, p.y() + epsilon - m * p.x());
            }
            std::cout << "DEBUG: Algorithm bMin=" << bMin << ", bMax=" << bMax 
                      << ", gap=" << (bMax - bMin) << "\n";
        }
    }
    
    // If algorithm says there's a line, verify it
    result.algorithmLineValid = true;
    if (result.algorithmResult) {
        auto line = sls.findStabbingLine();
        if (line.has_value()) {
            result.algorithmLineValid = verifyLineCoversAllPoints(points, *line, epsilon);
            if (!result.algorithmLineValid) {
                std::ostringstream oss;
                oss << "Line (slope=" << line->slope << ", intercept=" << line->intercept 
                    << ") doesn't cover all points!";
                // Find the first uncovered point
                for (const auto& p : points) {
                    double lineY = line->at(p.x());
                    double error = std::abs(p.y() - lineY);
                    if (error > epsilon + 1e-9) {
                        oss << " Point (" << p.x() << "," << p.y() << ") has error " << error;
                        break;
                    }
                }
                result.details = oss.str();
            }
        } else {
            result.algorithmLineValid = false;
            result.details = "hasStabbingLine() returned true but findStabbingLine() returned nullopt";
        }
    }
    
    // Build details if there's a mismatch
    if (result.algorithmResult != result.referenceResult) {
        std::ostringstream oss;
        oss << "Algorithm says " << (result.algorithmResult ? "YES" : "NO")
            << ", Reference says " << (result.referenceResult ? "YES" : "NO");
        result.details = oss.str();
    }
    
    return result;
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

TEST(StabbingLine, TrianglePattern) {
    SLS sls(0.5);
    // Triangle pattern - no line fits with small epsilon
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(5, 10));
    sls.insert(Point_2(10, 0));
    EXPECT_FALSE(sls.hasStabbingLine());
}

TEST(StabbingLine, WideSpread) {
    SLS sls(1.0);
    // Points spread across a wide y-range at same x
    sls.insert(Point_2(0, -100));
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(0, 100));
    // Gap = 200, epsilon = 1, so need gap <= 2 = impossible
    EXPECT_FALSE(sls.hasStabbingLine());
}

TEST(StabbingLine, ParabolicSpread) {
    SLS sls(0.5);
    // Points that follow a parabola - no line fits
    sls.insert(Point_2(-2, 4));
    sls.insert(Point_2(-1, 1));
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(1, 1));
    sls.insert(Point_2(2, 4));
    EXPECT_FALSE(sls.hasStabbingLine());
}

TEST(StabbingLine, BoxPattern) {
    SLS sls(0.1);
    // Four corners of a box - no line fits with small epsilon
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(0, 10));
    sls.insert(Point_2(10, 0));
    sls.insert(Point_2(10, 10));
    EXPECT_FALSE(sls.hasStabbingLine());
}

TEST(StabbingLine, SteepVsShallow) {
    SLS sls(0.1);
    // One set of points requires steep line, another requires shallow
    sls.insert(Point_2(0, 0));
    sls.insert(Point_2(1, 10));  // Suggests slope ~10
    sls.insert(Point_2(10, 1));  // Suggests slope ~0.1
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

// ============================================================================
// RIGOROUS VERIFICATION TESTS
// Compare algorithm against independent reference and verify all returned lines
// ============================================================================

// Test: verify that returned lines ALWAYS cover all points
TEST(StabbingLineRigorous, ReturnedLineAlwaysValid) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> noise(-0.8, 0.8);
    
    for (int trial = 0; trial < 100; ++trial) {
        std::vector<Point_2> points;
        double epsilon = 1.0;
        
        // Generate points near y = mx + b
        double slope = (rng() % 100 - 50) / 10.0;
        double intercept = (rng() % 100 - 50);
        int n = 5 + (rng() % 20);
        
        for (int i = 0; i < n; ++i) {
            double x = i * 2.0;
            double y = slope * x + intercept + noise(rng);
            points.emplace_back(x, y);
        }
        
        VerificationResult result = verifyStabbingLineQuery(points, epsilon);
        
        EXPECT_TRUE(result.algorithmLineValid) 
            << "Trial " << trial << ": " << result.details;
    }
}

// Test: algorithm matches reference on random VALID instances
TEST(StabbingLineRigorous, MatchesReferenceOnValid) {
    std::mt19937 rng(100);
    std::uniform_real_distribution<double> noise(-0.4, 0.4);
    
    for (int trial = 0; trial < 50; ++trial) {
        std::vector<Point_2> points;
        double epsilon = 0.5;
        
        // Generate points guaranteed to be near a line
        double slope = (rng() % 20 - 10) / 5.0;
        double intercept = rng() % 50;
        int n = 5 + (rng() % 15);
        
        for (int i = 0; i < n; ++i) {
            double x = i * 3.0;
            double y = slope * x + intercept + noise(rng);
            points.emplace_back(x, y);
        }
        
        VerificationResult result = verifyStabbingLineQuery(points, epsilon);
        
        EXPECT_TRUE(result.passed()) 
            << "Trial " << trial << ": " << result.details
            << " (algo=" << result.algorithmResult 
            << ", ref=" << result.referenceResult 
            << ", lineValid=" << result.algorithmLineValid << ")";
    }
}

// Test: algorithm matches reference on random INVALID instances
TEST(StabbingLineRigorous, MatchesReferenceOnInvalid) {
    std::mt19937 rng(200);
    std::uniform_real_distribution<double> yDist(-100, 100);
    
    for (int trial = 0; trial < 50; ++trial) {
        std::vector<Point_2> points;
        double epsilon = 0.5;  // Small epsilon
        
        // Generate scattered points that likely don't fit a line
        int n = 5 + (rng() % 10);
        for (int i = 0; i < n; ++i) {
            double x = i * 2.0;
            double y = yDist(rng);
            points.emplace_back(x, y);
        }
        
        VerificationResult result = verifyStabbingLineQuery(points, epsilon);
        
        EXPECT_TRUE(result.passed()) 
            << "Trial " << trial << ": " << result.details;
    }
}

// Test: rigorous verification across mixed valid/invalid instances
TEST(StabbingLineRigorous, MixedRandomInstances) {
    std::mt19937 rng(300);
    std::uniform_real_distribution<double> noise(-0.5, 0.5);
    std::uniform_real_distribution<double> wild(-50, 50);
    std::uniform_real_distribution<double> epsDist(0.5, 1.5);
    std::uniform_int_distribution<int> nDist(3, 14);
    std::uniform_int_distribution<int> slopeDist(-50, 49);
    std::uniform_int_distribution<int> interceptDist(-50, 49);
    
    int passCount = 0;
    int failCount = 0;
    
    for (int trial = 0; trial < 200; ++trial) {
        std::vector<Point_2> points;
        double epsilon = epsDist(rng);
        
        int n = nDist(rng);
        bool addNoise = rng() % 2 == 0;
        double slope = slopeDist(rng) / 10.0;
        double intercept = interceptDist(rng);
        
        for (int i = 0; i < n; ++i) {
            double x = i * 2.0;
            double y;
            if (addNoise) {
                y = slope * x + intercept + noise(rng);
            } else {
                y = wild(rng);
            }
            points.emplace_back(x, y);
        }
        
        VerificationResult result = verifyStabbingLineQuery(points, epsilon);
        
        EXPECT_TRUE(result.passed()) 
            << "Trial " << trial << ": " << result.details;
        
        if (result.algorithmResult) passCount++;
        else failCount++;
    }
    
    // Should have some of both valid and invalid
    std::cout << "Rigorous test: " << passCount << " valid, " << failCount << " invalid\n";
}

// Test: verify line validation is working correctly
TEST(StabbingLineRigorous, LineValidationWorks) {
    std::vector<Point_2> points = {
        Point_2(0, 0),
        Point_2(5, 5),
        Point_2(10, 10)
    };
    
    // Line y = x should cover all within epsilon=1
    EXPECT_TRUE(verifyLineCoversAllPoints(points, 1.0, 0.0, 1.0));
    
    // Line y = x + 5 should NOT cover (0,0) within epsilon=1
    EXPECT_FALSE(verifyLineCoversAllPoints(points, 1.0, 5.0, 1.0));
    
    // Line y = 0 should NOT cover (5,5) within epsilon=1
    EXPECT_FALSE(verifyLineCoversAllPoints(points, 0.0, 0.0, 1.0));
}

// Test: edge case - two points with same x
TEST(StabbingLineRigorous, SameXCoordinate) {
    for (double eps : {0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) {
        std::vector<Point_2> points = {
            Point_2(5, 0),
            Point_2(5, 3)
        };
        
        VerificationResult result = verifyStabbingLineQuery(points, eps);
        
        // Reference: no line if gap > 2*epsilon
        double gap = 3.0;
        bool expected = (gap <= 2.0 * eps);
        
        EXPECT_EQ(result.algorithmResult, expected) 
            << "eps=" << eps << ": " << result.details;
        EXPECT_TRUE(result.passed()) << "eps=" << eps;
    }
}

// Test: stress test with many random instances and full verification
TEST(StabbingLineRigorous, StressTestWithFullVerification) {
    std::mt19937 rng(999);
    std::uniform_real_distribution<double> noise(-0.3, 0.3);
    
    int totalTests = 0;
    int lineTests = 0;
    
    for (int trial = 0; trial < 100; ++trial) {
        std::vector<Point_2> points;
        double epsilon = 0.5;
        
        // Generate points near y = x
        int n = 10 + (rng() % 50);
        for (int i = 0; i < n; ++i) {
            double x = i;
            double y = x + noise(rng);
            points.emplace_back(x, y);
        }
        
        totalTests++;
        VerificationResult result = verifyStabbingLineQuery(points, epsilon);
        
        EXPECT_TRUE(result.passed()) 
            << "Trial " << trial << " (n=" << n << "): " << result.details;
        
        if (result.algorithmResult && result.algorithmLineValid) {
            lineTests++;
        }
    }
    
    std::cout << "Stress test: " << totalTests << " tests, " << lineTests << " valid lines\n";
}

// Test: verify dynamic operations maintain correctness
TEST(StabbingLineRigorous, DynamicOperationsVerified) {
    std::mt19937 rng(555);
    std::uniform_real_distribution<double> noise(-0.2, 0.2);
    std::uniform_int_distribution<int> action(0, 2);
    
    double epsilon = 0.5;
    SLS sls(epsilon);
    std::vector<Point_2> points;
    
    // Build initial set
    for (int i = 0; i < 10; ++i) {
        double x = i;
        double y = x + noise(rng);
        points.emplace_back(x, y);
        sls.insert(points.back());
    }
    
    // Perform operations and verify after each
    for (int iter = 0; iter < 50; ++iter) {
        int act = action(rng);
        
        if (points.empty() || act == 0) {
            // Insert
            double x = points.empty() ? 0 : points.back().x() + 1;
            double y = x + noise(rng);
            points.emplace_back(x, y);
            sls.insert(points.back());
        } else if (act == 1 && points.size() > 2) {
            // Remove
            size_t idx = rng() % points.size();
            sls.remove(points[idx]);
            points.erase(points.begin() + idx);
        }
        // act == 2: just verify
        
        // Verify current state
        bool algResult = sls.hasStabbingLine();
        bool refResult = bruteForceFindStabbingLine(points, epsilon);
        
        EXPECT_EQ(algResult, refResult) 
            << "Iter " << iter << ": algo=" << algResult << ", ref=" << refResult;
        
        // If there's a line, verify it covers all points
        if (algResult) {
            auto line = sls.findStabbingLine();
            ASSERT_TRUE(line.has_value());
            EXPECT_TRUE(verifyLineCoversAllPoints(points, *line, epsilon))
                << "Iter " << iter << ": Line doesn't cover all points";
        }
    }
}

// Test: verify half-plane intersection reference is correct
TEST(StabbingLineRigorous, ReferenceImplementationIsCorrect) {
    // Known case: three collinear points
    {
        std::vector<Point_2> pts = {Point_2(0,0), Point_2(5,5), Point_2(10,10)};
        HalfPlaneIntersection hpi(pts, 1.0);
        EXPECT_TRUE(hpi.hasStabbingLine());
        auto [ok, m, b] = hpi.findLine();
        EXPECT_TRUE(ok);
        EXPECT_TRUE(verifyLineCoversAllPoints(pts, m, b, 1.0));
    }
    
    // Known case: zigzag (no line with small epsilon)
    {
        std::vector<Point_2> pts = {Point_2(0,0), Point_2(1,10), Point_2(2,0), Point_2(3,10)};
        HalfPlaneIntersection hpi(pts, 1.0);
        EXPECT_FALSE(hpi.hasStabbingLine());
    }
    
    // Known case: two points same x, far apart
    {
        std::vector<Point_2> pts = {Point_2(5,0), Point_2(5,10)};
        HalfPlaneIntersection hpi(pts, 1.0);
        EXPECT_FALSE(hpi.hasStabbingLine());
    }
    
    // Known case: two points same x, close together
    {
        std::vector<Point_2> pts = {Point_2(5,0), Point_2(5,1)};
        HalfPlaneIntersection hpi(pts, 1.0);
        EXPECT_TRUE(hpi.hasStabbingLine());
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


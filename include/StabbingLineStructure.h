//
// StabbingLineStructure - Query for epsilon-stabbing lines in O(log² n)
//
// A "stabbing line" is a line that passes within vertical distance epsilon
// of every point in the set. This structure maintains two convex hulls
// (ceiling and floor) to efficiently query for stabbing lines.
//
// Algorithm:
//   1. Tree-based search: Find critical vertex where slope difference changes sign
//   2. Subgradient intersection: Compute valid slope range at critical vertex
//   3. Line construction: Use slope from intersection, verify at endpoints
//   
// Time Complexity:
//   - insert/remove: O(log² n)
//   - hasStabbingLine/findStabbingLine: O(log² n)
//   - split/join: O(log² n)
//   - build: O(n)
//

#ifndef DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H
#define DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

#include "CHTree.h"
#include <optional>
#include <cmath>
#include <limits>
#include <iostream>

// Debug flag - set to true to enable verbose output
#ifndef STABBING_DEBUG
#define STABBING_DEBUG false
#endif

template<class Traits>
struct StabbingLine {
    using Point = typename Traits::Point_2;
    double slope;
    double intercept;  // y = slope * x + intercept
    
    StabbingLine() : slope(0), intercept(0) {}
    StabbingLine(double s, double i) : slope(s), intercept(i) {}
    
    double at(double x) const { return slope * x + intercept; }
    
    bool coversPoint(const Point& p, double epsilon) const {
        double lineY = at(p.x());
        return std::abs(p.y() - lineY) <= epsilon + 1e-9;
    }
};

template<class Traits>
class StabbingLineStructure {
public:
    using Point = typename Traits::Point_2;
    using Line = StabbingLine<Traits>;
    
private:
    CHTree<Traits> ceiling;  // stores (x, y+ε), query lower hull
    CHTree<Traits> floor;    // stores (x, y-ε), query upper hull
    double epsilon;
    size_t pointCount = 0;
    
    Point shiftUp(const Point& p) const {
        return Point(p.x(), p.y() + epsilon);
    }
    
    Point shiftDown(const Point& p) const {
        return Point(p.x(), p.y() - epsilon);
    }
    
public:
    explicit StabbingLineStructure(double eps = 1.0) : epsilon(eps) {}
    
    double getEpsilon() const { return epsilon; }
    size_t size() const { return pointCount; }
    bool empty() const { return pointCount == 0; }
    
    void insert(const Point& p) {
        ceiling.insert(shiftUp(p));
        floor.insert(shiftDown(p));
        ++pointCount;
    }
    
    void remove(const Point& p) {
        ceiling.remove(shiftUp(p));
        floor.remove(shiftDown(p));
        if (pointCount > 0) --pointCount;
    }
    
    void build(const std::vector<Point>& sortedPoints) {
        std::vector<Point> ceilingPoints, floorPoints;
        ceilingPoints.reserve(sortedPoints.size());
        floorPoints.reserve(sortedPoints.size());
        
        for (const auto& p : sortedPoints) {
            ceilingPoints.push_back(shiftUp(p));
            floorPoints.push_back(shiftDown(p));
        }
        
        ceiling.build(ceilingPoints);
        floor.build(floorPoints);
        pointCount = sortedPoints.size();
    }
    
    StabbingLineStructure split(double splitX) {
        StabbingLineStructure right(epsilon);
        auto ceilingRight = ceiling.split(splitX);
        auto floorRight = floor.split(splitX);
        right.ceiling = std::move(ceilingRight);
        right.floor = std::move(floorRight);
        right.pointCount = right.ceiling.size();
        pointCount = ceiling.size();
        return right;
    }
    
    void join(StabbingLineStructure& other) {
        ceiling.join(other.ceiling);
        floor.join(other.floor);
        pointCount = ceiling.size();
        other.pointCount = 0;
    }
    
    // ========================================================================
    // Core Query: Find a stabbing line in PURE O(log² n)
    // ========================================================================
    //
    // Expert's Algorithm:
    // 1. Traverse ceiling tree to find candidate vertex (using floor slopes)
    // 2. Traverse floor tree to find candidate vertex (using ceiling slopes)
    // 3. Evaluate gap at both candidates, pick the one with minimum gap
    // 4. If gap >= 0, construct line and return; else return nullopt
    //
    // Each tree traversal is O(log n) with O(log n) slope query per step = O(log² n)
    // ========================================================================
    
    std::optional<Line> findStabbingLine() const {
        if constexpr (STABBING_DEBUG) {
            std::cout << "=== findStabbingLine() ===" << std::endl;
            std::cout << "pointCount = " << pointCount << ", epsilon = " << epsilon << std::endl;
        }
        
        // Edge case: empty set
        if (pointCount == 0) {
            return Line(0.0, 0.0);
        }
        
        // Get x-range
        auto [minX, maxX] = ceiling.getXRange();
        if constexpr (STABBING_DEBUG) {
            std::cout << "x-range: [" << minX << ", " << maxX << "]" << std::endl;
        }
        
        if (minX > maxX) {
            return Line(0.0, 0.0);
        }
        
        // Edge case: single point or degenerate range
        if (std::abs(maxX - minX) < 1e-12) {
            double yCeil = ceiling.evaluateLowerHullAt(minX);
            double yFloor = floor.evaluateUpperHullAt(minX);
            double gap = yCeil - yFloor;
            if (gap >= -1e-9) {
                return Line(0.0, (yCeil + yFloor) / 2.0);
            }
            return std::nullopt;
        }
        
        // Phase A: Find candidate vertex by traversing ceiling tree
        auto [xCeil, slopeCeil] = ceiling.findCriticalVertexFromCeiling(floor);
        
        // Phase B: Find candidate vertex by traversing floor tree
        auto [xFloor, slopeFloor] = floor.findCriticalVertexFromFloor(ceiling);
        
        if constexpr (STABBING_DEBUG) {
            std::cout << "Ceiling candidate: x=" << xCeil << ", slope=" << slopeCeil << std::endl;
            std::cout << "Floor candidate: x=" << xFloor << ", slope=" << slopeFloor << std::endl;
        }
        
        // Phase C: Evaluate gap at both candidates and endpoints
        double gapCeil = computeGap(xCeil);
        double gapFloor = computeGap(xFloor);
        double gapMin = computeGap(minX);
        double gapMax = computeGap(maxX);
        
        if constexpr (STABBING_DEBUG) {
            std::cout << "Gap at ceiling candidate: " << gapCeil << std::endl;
            std::cout << "Gap at floor candidate: " << gapFloor << std::endl;
            std::cout << "Gap at minX: " << gapMin << std::endl;
            std::cout << "Gap at maxX: " << gapMax << std::endl;
        }
        
        // Find the minimum gap and its location
        double bestGap = gapCeil;
        double bestX = xCeil;
        
        if (gapFloor < bestGap) { bestGap = gapFloor; bestX = xFloor; }
        if (gapMin < bestGap) { bestGap = gapMin; bestX = minX; }
        if (gapMax < bestGap) { bestGap = gapMax; bestX = maxX; }
        
        if constexpr (STABBING_DEBUG) {
            std::cout << "Best: x=" << bestX << ", gap=" << bestGap << std::endl;
        }
        
        // Check if stabbing line exists
        if (bestGap < -1e-9) {
            if constexpr (STABBING_DEBUG) {
                std::cout << "No stabbing line exists (gap < 0)" << std::endl;
            }
            return std::nullopt;
        }
        
        // Construct line at the critical point using SUBGRADIENT INTERSECTION
        // At a vertex, each hull has a RANGE of valid tangent slopes
        // We need to find a slope in the intersection of both ranges
        
        double yCeil = ceiling.evaluateLowerHullAt(bestX);
        double yFloor = floor.evaluateUpperHullAt(bestX);
        double yMid = (yCeil + yFloor) / 2.0;
        
        // Get slope ranges at the critical point
        // Ceiling (lower hull): returns [leftSlope, rightSlope] where valid tangent is in this range
        // Floor (upper hull): returns [leftSlope, rightSlope] - for upper hull, valid tangent is in [right, left]
        auto [ceilLeft, ceilRight] = ceiling.getLowerHullSlopeRange(bestX);
        auto [floorLeft, floorRight] = floor.getUpperHullSlopeRange(bestX);
        
        if constexpr (STABBING_DEBUG) {
            std::cout << "Ceiling slope range: [" << ceilLeft << ", " << ceilRight << "]" << std::endl;
            std::cout << "Floor slope range: [" << floorLeft << ", " << floorRight << "]" << std::endl;
        }
        
        // Both hulls return ranges in [min, max] format
        // Ceiling (lower hull): valid slopes in [ceilLeft, ceilRight]
        // Floor (upper hull): valid slopes in [floorLeft, floorRight]
        // Intersect by taking max of left bounds and min of right bounds
        double minValid, maxValid;
        
        // Handle infinity cases
        if (std::isinf(ceilLeft) && ceilLeft < 0) ceilLeft = -1e15;
        if (std::isinf(ceilRight) && ceilRight > 0) ceilRight = 1e15;
        if (std::isinf(floorLeft) && floorLeft < 0) floorLeft = -1e15;
        if (std::isinf(floorRight) && floorRight > 0) floorRight = 1e15;
        
        // Intersect: valid slope must be in BOTH ranges
        minValid = std::max(ceilLeft, floorLeft);
        maxValid = std::min(ceilRight, floorRight);
        
        if constexpr (STABBING_DEBUG) {
            std::cout << "Valid slope intersection: [" << minValid << ", " << maxValid << "]" << std::endl;
        }
        
        // Pick slope from the valid range
        double slope;
        if (minValid <= maxValid + 1e-9) {
            slope = (minValid + maxValid) / 2.0;  // Midpoint of valid range
        } else {
            // Ranges don't overlap - no valid stabbing line exists
            if constexpr (STABBING_DEBUG) {
                std::cout << "No slope intersection! Returning nullopt" << std::endl;
            }
            return std::nullopt;
        }
        
        // Compute intercept using expert's formula:
        // b = (y_ceil - gap/2) - m * x
        double gap = yCeil - yFloor;
        double intercept = (yCeil - gap / 2.0) - slope * bestX;
        
        // CRITICAL: Verify the line maintains positive gap at endpoints
        // The line y = slope * x + intercept should satisfy:
        //   line(x) <= ceiling(x)  AND  line(x) >= floor(x)  for all x
        // This is equivalent to:
        //   ceiling(x) - line(x) >= 0  AND  line(x) - floor(x) >= 0
        
        auto verifyGapWithSlope = [&](double x) -> bool {
            double lineY = slope * x + intercept;
            double yCeilAt = ceiling.evaluateLowerHullAt(x);
            double yFloorAt = floor.evaluateUpperHullAt(x);
            // Line must be between floor and ceiling
            return (lineY <= yCeilAt + 1e-9) && (lineY >= yFloorAt - 1e-9);
        };
        
        // Check endpoints and critical point
        if (!verifyGapWithSlope(minX) || !verifyGapWithSlope(maxX) || !verifyGapWithSlope(bestX)) {
            if constexpr (STABBING_DEBUG) {
                std::cout << "Constructed line fails gap check at endpoints! Returning nullopt" << std::endl;
            }
            return std::nullopt;
        }
        
        if constexpr (STABBING_DEBUG) {
            std::cout << "Constructed line: y = " << slope << "x + " << intercept << std::endl;
        }
        
        return Line(slope, intercept);
    }
    
    bool hasStabbingLine() const {
        return findStabbingLine().has_value();
    }
    
    double getMinimumGap() const {
        if (pointCount == 0) return std::numeric_limits<double>::infinity();
        
        auto [minX, maxX] = ceiling.getXRange();
        if (minX > maxX) return std::numeric_limits<double>::infinity();
        
        auto [xCeil, _1] = ceiling.findCriticalVertexFromCeiling(floor);
        auto [xFloor, _2] = floor.findCriticalVertexFromFloor(ceiling);
        
        double gap = computeGap(xCeil);
        gap = std::min(gap, computeGap(xFloor));
        gap = std::min(gap, computeGap(minX));
        gap = std::min(gap, computeGap(maxX));
        return gap;
    }
    
    void printDebugInfo() const {
        std::cout << "=== StabbingLineStructure Debug Info ===" << std::endl;
        std::cout << "pointCount: " << pointCount << std::endl;
        std::cout << "epsilon: " << epsilon << std::endl;
        
        if (pointCount == 0) {
            std::cout << "(empty)" << std::endl;
            return;
        }
        
        auto [minX, maxX] = ceiling.getXRange();
        std::cout << "X range: [" << minX << ", " << maxX << "]" << std::endl;
        
        std::cout << "Hull samples:" << std::endl;
        for (int i = 0; i <= 10; ++i) {
            double x = minX + (maxX - minX) * i / 10.0;
            double yCeil = ceiling.evaluateLowerHullAt(x);
            double yFloor = floor.evaluateUpperHullAt(x);
            double slopeCeil = ceiling.getLowerHullSlope(x);
            double slopeFloor = floor.getUpperHullSlope(x);
            std::cout << "  x=" << x << ": ceil=" << yCeil << "(m=" << slopeCeil 
                      << "), floor=" << yFloor << "(m=" << slopeFloor 
                      << "), gap=" << (yCeil - yFloor) << std::endl;
        }
    }
    
private:
    double computeGap(double x) const {
        double yCeiling = ceiling.evaluateLowerHullAt(x);
        double yFloor = floor.evaluateUpperHullAt(x);
        return yCeiling - yFloor;
    }
};

#endif // DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

/**
 * @file StabbingLineStructure.h
 * @brief Query for epsilon-stabbing lines in O(log² n).
 * 
 * A "stabbing line" is a line that passes within vertical distance epsilon
 * of every point in the set. This structure maintains a single convex hull
 * and applies epsilon adjustments during queries.
 * 
 * Time Complexities:
 *   - insert/remove: O(log² n)
 *   - hasStabbingLine/findStabbingLine: O(log² n)
 *   - split/join: O(log² n)
 *   - build: O(n)
 */

#ifndef DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H
#define DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

#include "CHTree.h"
#include <optional>
#include <cmath>
#include <limits>

/**
 * @brief Represents a stabbing line y = slope * x + intercept.
 */
template<class Traits>
struct StabbingLine {
    using Point = typename Traits::Point_2;
    double slope;
    double intercept;
    
    StabbingLine() : slope(0), intercept(0) {}
    StabbingLine(double s, double i) : slope(s), intercept(i) {}
    
    double at(double x) const { return slope * x + intercept; }
    
    bool coversPoint(const Point& p, double epsilon) const {
        double lineY = at(p.x());
        return std::abs(p.y() - lineY) <= epsilon + 1e-9;
    }
};

/**
 * @brief Data structure for querying epsilon-stabbing lines.
 * 
 * Uses a single convex hull and applies epsilon adjustments at query time:
 * - For "ceiling" queries (upper boundary): lower hull y + epsilon
 * - For "floor" queries (lower boundary): upper hull y - epsilon
 * 
 * @tparam Traits CGAL-style kernel.
 */
template<class Traits>
class StabbingLineStructure {
public:
    using Point = typename Traits::Point_2;
    using Line = StabbingLine<Traits>;
    
private:
    CHTree<Traits> hull;  ///< Single hull storing original points
    double epsilon;
    size_t pointCount = 0;
    
public:
    explicit StabbingLineStructure(double eps = 1.0) : epsilon(eps) {}
    
    double getEpsilon() const { return epsilon; }
    size_t size() const { return pointCount; }
    bool empty() const { return pointCount == 0; }
    
    /** @brief Insert a point. O(log² n) */
    void insert(const Point& p) {
        hull.insert(p);
        ++pointCount;
    }
    
    /** @brief Remove a point. O(log² n) */
    void remove(const Point& p) {
        hull.remove(p);
        if (pointCount > 0) --pointCount;
    }
    
    /** @brief Build from sorted points. O(n) */
    void build(const std::vector<Point>& sortedPoints) {
        hull.build(sortedPoints);
        pointCount = sortedPoints.size();
    }
    
    /** @brief Split at x-coordinate. O(log² n) */
    StabbingLineStructure split(double splitX) {
        StabbingLineStructure right(epsilon);
        auto hullRight = hull.split(splitX);
        right.hull = std::move(hullRight);
        right.pointCount = right.hull.size();
        pointCount = hull.size();
        return right;
    }
    
    /** @brief Join with another structure. O(log² n) */
    void join(StabbingLineStructure& other) {
        hull.join(other.hull);
        pointCount = hull.size();
        other.pointCount = 0;
    }
    
    /**
     * @brief Find a stabbing line if one exists.
     * @return Line if exists, nullopt otherwise.
     * @note Time complexity: O(log² n)
     */
    std::optional<Line> findStabbingLine() const {
        if (pointCount == 0) {
            return Line(0.0, 0.0);
        }
        
        auto [minX, maxX] = hull.getXRange();
        
        if (minX > maxX) {
            return Line(0.0, 0.0);
        }
        
        if (std::abs(maxX - minX) < 1e-12) {
            // Single x-coordinate: ceiling = lower hull + epsilon, floor = upper hull - epsilon
            double yCeil = hull.evaluateLowerHullAt(minX) + epsilon;
            double yFloor = hull.evaluateUpperHullAt(minX) - epsilon;
            double gap = yCeil - yFloor;
            if (gap >= -1e-9) {
                return Line(0.0, (yCeil + yFloor) / 2.0);
            }
            return std::nullopt;
        }
        
        // Find critical vertices - comparing hull against itself since both 
        // ceiling and floor are derived from the same hull
        auto [xCeil, slopeCeil] = hull.findCriticalVertexFromCeiling(hull);
        auto [xFloor, slopeFloor] = hull.findCriticalVertexFromFloor(hull);
        
        double gapCeil = computeGap(xCeil);
        double gapFloor = computeGap(xFloor);
        double gapMin = computeGap(minX);
        double gapMax = computeGap(maxX);
        
        double bestGap = gapCeil;
        double bestX = xCeil;
        
        if (gapFloor < bestGap) { bestGap = gapFloor; bestX = xFloor; }
        if (gapMin < bestGap) { bestGap = gapMin; bestX = minX; }
        if (gapMax < bestGap) { bestGap = gapMax; bestX = maxX; }
        
        if (bestGap < -1e-9) {
            return std::nullopt;
        }
        
        // Ceiling = lower hull + epsilon, Floor = upper hull - epsilon
        double yCeil = hull.evaluateLowerHullAt(bestX) + epsilon;
        double yFloor = hull.evaluateUpperHullAt(bestX) - epsilon;
        double yMid = (yCeil + yFloor) / 2.0;
        
        // Get slope ranges - for ceiling use lower hull, for floor use upper hull
        auto [ceilLeft, ceilRight] = hull.getLowerHullSlopeRange(bestX);
        auto [floorLeft, floorRight] = hull.getUpperHullSlopeRange(bestX);
        
        double minValid, maxValid;
        
        if (std::isinf(ceilLeft) && ceilLeft < 0) ceilLeft = -1e15;
        if (std::isinf(ceilRight) && ceilRight > 0) ceilRight = 1e15;
        if (std::isinf(floorLeft) && floorLeft < 0) floorLeft = -1e15;
        if (std::isinf(floorRight) && floorRight > 0) floorRight = 1e15;
        
        minValid = std::max(ceilLeft, floorLeft);
        maxValid = std::min(ceilRight, floorRight);
        
        double slope;
        if (minValid <= maxValid + 1e-9) {
            slope = (minValid + maxValid) / 2.0;
        } else {
            return std::nullopt;
        }
        
        double gap = yCeil - yFloor;
        double intercept = (yCeil - gap / 2.0) - slope * bestX;
        
        auto verifyGapWithSlope = [&](double x) -> bool {
            double lineY = slope * x + intercept;
            double yCeilAt = hull.evaluateLowerHullAt(x) + epsilon;
            double yFloorAt = hull.evaluateUpperHullAt(x) - epsilon;
            return (lineY <= yCeilAt + 1e-9) && (lineY >= yFloorAt - 1e-9);
        };
        
        if (!verifyGapWithSlope(minX) || !verifyGapWithSlope(maxX) || !verifyGapWithSlope(bestX)) {
            return std::nullopt;
        }
        
        return Line(slope, intercept);
    }
    
    /** @brief Check if a stabbing line exists. O(log² n) */
    bool hasStabbingLine() const {
        return findStabbingLine().has_value();
    }
    
    /** @brief Get the minimum gap (ceiling - floor). */
    double getMinimumGap() const {
        if (pointCount == 0) return std::numeric_limits<double>::infinity();
        
        auto [minX, maxX] = hull.getXRange();
        if (minX > maxX) return std::numeric_limits<double>::infinity();
        
        auto [xCeil, _1] = hull.findCriticalVertexFromCeiling(hull);
        auto [xFloor, _2] = hull.findCriticalVertexFromFloor(hull);
        
        double gap = computeGap(xCeil);
        gap = std::min(gap, computeGap(xFloor));
        gap = std::min(gap, computeGap(minX));
        gap = std::min(gap, computeGap(maxX));
        return gap;
    }
    
private:
    /**
     * @brief Compute the gap between ceiling and floor at x.
     * 
     * Ceiling = lower hull + epsilon (upper boundary of epsilon band)
     * Floor = upper hull - epsilon (lower boundary of epsilon band)
     * Gap = ceiling - floor = (lower hull + ε) - (upper hull - ε) = lower hull - upper hull + 2ε
     */
    double computeGap(double x) const {
        double yLower = hull.evaluateLowerHullAt(x);
        double yUpper = hull.evaluateUpperHullAt(x);
        // Ceiling = yLower + epsilon, Floor = yUpper - epsilon
        // Gap = ceiling - floor = (yLower + ε) - (yUpper - ε) = yLower - yUpper + 2ε
        return (yLower + epsilon) - (yUpper - epsilon);
    }
};

#endif // DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

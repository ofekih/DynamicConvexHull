/**
 * @file StabbingLineStructure.h
 * @brief Query for epsilon-stabbing lines in O(log² n).
 * 
 * A "stabbing line" is a line that passes within vertical distance epsilon
 * of every point in the set. This structure maintains two convex hulls
 * (ceiling and floor) to efficiently query for stabbing lines.
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
 * @tparam Traits CGAL-style kernel.
 */
template<class Traits>
class StabbingLineStructure {
public:
    using Point = typename Traits::Point_2;
    using Line = StabbingLine<Traits>;
    
private:
    CHTree<Traits> ceiling;  ///< Stores (x, y+ε), query lower hull
    CHTree<Traits> floor;    ///< Stores (x, y-ε), query upper hull
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
    
    /** @brief Insert a point. O(log² n) */
    void insert(const Point& p) {
        ceiling.insert(shiftUp(p));
        floor.insert(shiftDown(p));
        ++pointCount;
    }
    
    /** @brief Remove a point. O(log² n) */
    void remove(const Point& p) {
        ceiling.remove(shiftUp(p));
        floor.remove(shiftDown(p));
        if (pointCount > 0) --pointCount;
    }
    
    /** @brief Build from sorted points. O(n) */
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
    
    /** @brief Split at x-coordinate. O(log² n) */
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
    
    /** @brief Join with another structure. O(log² n) */
    void join(StabbingLineStructure& other) {
        ceiling.join(other.ceiling);
        floor.join(other.floor);
        pointCount = ceiling.size();
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
        
        auto [minX, maxX] = ceiling.getXRange();
        
        if (minX > maxX) {
            return Line(0.0, 0.0);
        }
        
        if (std::abs(maxX - minX) < 1e-12) {
            double yCeil = ceiling.evaluateLowerHullAt(minX);
            double yFloor = floor.evaluateUpperHullAt(minX);
            double gap = yCeil - yFloor;
            if (gap >= -1e-9) {
                return Line(0.0, (yCeil + yFloor) / 2.0);
            }
            return std::nullopt;
        }
        
        auto [xCeil, slopeCeil] = ceiling.findCriticalVertexFromCeiling(floor);
        auto [xFloor, slopeFloor] = floor.findCriticalVertexFromFloor(ceiling);
        
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
        
        double yCeil = ceiling.evaluateLowerHullAt(bestX);
        double yFloor = floor.evaluateUpperHullAt(bestX);
        double yMid = (yCeil + yFloor) / 2.0;
        
        auto [ceilLeft, ceilRight] = ceiling.getLowerHullSlopeRange(bestX);
        auto [floorLeft, floorRight] = floor.getUpperHullSlopeRange(bestX);
        
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
            double yCeilAt = ceiling.evaluateLowerHullAt(x);
            double yFloorAt = floor.evaluateUpperHullAt(x);
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
    
private:
    double computeGap(double x) const {
        double yCeiling = ceiling.evaluateLowerHullAt(x);
        double yFloor = floor.evaluateUpperHullAt(x);
        return yCeiling - yFloor;
    }
};

#endif // DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

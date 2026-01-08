//
// StabbingLineStructure - Query for epsilon-stabbing lines in O(log² n)
//
// A "stabbing line" is a line that passes within vertical distance epsilon
// of every point in the set. This structure maintains two convex hulls
// (ceiling and floor) and uses binary search on the gap function to find
// the optimal separating line.
//
// Algorithm:
//   Gap(x) = LowerHull_ceiling(x) - UpperHull_floor(x)
//   A stabbing line exists iff min(Gap(x)) >= 0
//   
// Time Complexity:
//   - insert/remove: O(log² n)
//   - findStabbingLine: O(log² n)
//   - split/join: O(log² n)
//   - build: O(n)
//

#ifndef DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H
#define DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

#include "CHTree.h"
#include <optional>
#include <cmath>

template<class Traits>
struct StabbingLine {
    using Point = typename Traits::Point_2;
    double slope;
    double intercept;  // y = slope * x + intercept
    
    StabbingLine() : slope(0), intercept(0) {}
    StabbingLine(double s, double i) : slope(s), intercept(i) {}
    
    // Evaluate y-value at x
    double at(double x) const { return slope * x + intercept; }
    
    // Check if a point is within epsilon of this line
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
    
    // Helper to create shifted points
    Point shiftUp(const Point& p) const {
        return Point(p.x(), p.y() + epsilon);
    }
    
    Point shiftDown(const Point& p) const {
        return Point(p.x(), p.y() - epsilon);
    }
    
public:
    explicit StabbingLineStructure(double eps = 1.0) : epsilon(eps) {}
    
    // Get epsilon value
    double getEpsilon() const { return epsilon; }
    
    // Get number of points
    size_t size() const { return pointCount; }
    
    // Check if empty
    bool empty() const { return pointCount == 0; }
    
    // Insert a point into both trees
    void insert(const Point& p) {
        ceiling.insert(shiftUp(p));
        floor.insert(shiftDown(p));
        ++pointCount;
    }
    
    // Remove a point from both trees
    void remove(const Point& p) {
        ceiling.remove(shiftUp(p));
        floor.remove(shiftDown(p));
        if (pointCount > 0) --pointCount;
    }
    
    // Build from sorted points in O(n) time
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
    
    // Split at x-coordinate, returning the right portion
    StabbingLineStructure split(double splitX) {
        StabbingLineStructure right(epsilon);
        
        auto ceilingRight = ceiling.split(splitX);
        auto floorRight = floor.split(splitX);
        
        right.ceiling = std::move(ceilingRight);
        right.floor = std::move(floorRight);
        
        // Update point counts (approximate - actual may differ)
        right.pointCount = right.ceiling.size();
        pointCount = ceiling.size();
        
        return right;
    }
    
    // Join with another structure (other's points must have x > this's points)
    void join(StabbingLineStructure& other) {
        ceiling.join(other.ceiling);
        floor.join(other.floor);
        pointCount += other.pointCount;
        other.pointCount = 0;
    }
    
    // ========================================================================
    // Core Query: Find a stabbing line in O(log² n)
    // ========================================================================
    
    // Find a line that passes within epsilon of all points.
    // Returns std::nullopt if no such line exists.
    std::optional<Line> findStabbingLine() const {
        // Edge case: empty set - any line works
        if (pointCount == 0) {
            return Line(0.0, 0.0);  // y = 0
        }
        
        // Edge case: single point - horizontal line through it
        if (pointCount == 1) {
            auto [minX, maxX] = ceiling.getXRange();
            double y = (ceiling.evaluateLowerHullAt(minX) + 
                       floor.evaluateUpperHullAt(minX)) / 2.0;
            return Line(0.0, y);
        }
        
        // Get x-range (should be same for both trees)
        auto [minX, maxX] = ceiling.getXRange();
        if (minX > maxX) {
            return Line(0.0, 0.0);  // Empty range
        }
        
        // Binary search to find the critical x where gap is minimized
        // The gap function Gap(x) = LowerHull_ceiling(x) - UpperHull_floor(x)
        // is convex, so we can use ternary search or slope comparison
        
        double criticalX;
        double criticalGap;
        
        // Use ternary search over x-range to find minimum gap
        double lo = minX, hi = maxX;
        const int MAX_ITERATIONS = 100;
        const double TOLERANCE = 1e-10;
        
        for (int iter = 0; iter < MAX_ITERATIONS && (hi - lo) > TOLERANCE; ++iter) {
            double mid1 = lo + (hi - lo) / 3.0;
            double mid2 = hi - (hi - lo) / 3.0;
            
            double gap1 = computeGap(mid1);
            double gap2 = computeGap(mid2);
            
            if (gap1 > gap2) {
                lo = mid1;
            } else {
                hi = mid2;
            }
        }
        
        criticalX = (lo + hi) / 2.0;
        criticalGap = computeGap(criticalX);
        
        // Check if a valid stabbing line exists
        if (criticalGap < -1e-9) {
            return std::nullopt;  // No stabbing line exists
        }
        
        // Construct the stabbing line
        // Use the average slope of the two hulls at the critical point
        double slopeCeiling = ceiling.getLowerHullSlope(criticalX);
        double slopeFloor = floor.getUpperHullSlope(criticalX);
        
        // Handle infinite slopes (vertical segments)
        double slope;
        if (std::isinf(slopeCeiling) || std::isinf(slopeFloor)) {
            slope = std::isinf(slopeCeiling) ? slopeFloor : slopeCeiling;
            if (std::isinf(slope)) slope = 0.0;
        } else {
            slope = (slopeCeiling + slopeFloor) / 2.0;
        }
        
        // Compute y at critical x: midpoint of the gap
        double yCeiling = ceiling.evaluateLowerHullAt(criticalX);
        double yFloor = floor.evaluateUpperHullAt(criticalX);
        double yMid = (yCeiling + yFloor) / 2.0;
        
        // Compute intercept: y = slope * x + intercept => intercept = y - slope * x
        double intercept = yMid - slope * criticalX;
        
        return Line(slope, intercept);
    }
    
    // Convenience method: just check if a stabbing line exists
    bool hasStabbingLine() const {
        return findStabbingLine().has_value();
    }
    
    // Get the minimum gap between the hulls (for debugging/analysis)
    double getMinimumGap() const {
        if (pointCount == 0) return std::numeric_limits<double>::infinity();
        if (pointCount == 1) return 2.0 * epsilon;
        
        auto [minX, maxX] = ceiling.getXRange();
        if (minX > maxX) return std::numeric_limits<double>::infinity();
        
        // Ternary search for minimum
        double lo = minX, hi = maxX;
        for (int iter = 0; iter < 100 && (hi - lo) > 1e-10; ++iter) {
            double mid1 = lo + (hi - lo) / 3.0;
            double mid2 = hi - (hi - lo) / 3.0;
            
            if (computeGap(mid1) > computeGap(mid2)) {
                lo = mid1;
            } else {
                hi = mid2;
            }
        }
        
        return computeGap((lo + hi) / 2.0);
    }
    
private:
    // Compute gap at x: LowerHull_ceiling(x) - UpperHull_floor(x)
    // Positive gap means valid stabbing line exists at this x
    double computeGap(double x) const {
        double yCeiling = ceiling.evaluateLowerHullAt(x);
        double yFloor = floor.evaluateUpperHullAt(x);
        return yCeiling - yFloor;
    }
};

#endif // DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

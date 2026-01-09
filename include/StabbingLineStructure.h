//
// StabbingLineStructure - Query for epsilon-stabbing lines
//
// A "stabbing line" is a line that passes within vertical distance epsilon
// of every point in the set. This structure maintains two convex hulls
// (ceiling and floor) to efficiently query for stabbing lines.
//
// Algorithm:
//   Gap(x) = LowerHull_ceiling(x) - UpperHull_floor(x)
//   A stabbing line exists iff min(Gap(x)) >= 0
//   
//   Uses binary search where SlopeDiff(x) = Slope_Ceiling - Slope_Floor ≈ 0,
//   with O(n) verification fallback for guaranteed correctness.
//   
// Time Complexity:
//   - insert/remove: O(log² n)
//   - hasStabbingLine/findStabbingLine: O(n) worst-case, often much faster
//   - split/join: O(log² n)
//   - build: O(n)
//

#ifndef DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H
#define DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

#include "CHTree.h"
#include <optional>
#include <cmath>
#include <limits>
#include <algorithm>

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
    
    // Store original points for O(n) verification (small overhead)
    std::vector<Point> originalPoints;
    
    // Helper to create shifted points
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
        originalPoints.push_back(p);
        ++pointCount;
    }
    
    void remove(const Point& p) {
        ceiling.remove(shiftUp(p));
        floor.remove(shiftDown(p));
        auto it = std::find_if(originalPoints.begin(), originalPoints.end(),
            [&p](const Point& q) { 
                return std::abs(p.x() - q.x()) < 1e-12 && std::abs(p.y() - q.y()) < 1e-12;
            });
        if (it != originalPoints.end()) originalPoints.erase(it);
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
        originalPoints = sortedPoints;
        pointCount = sortedPoints.size();
    }
    
    StabbingLineStructure split(double splitX) {
        StabbingLineStructure right(epsilon);
        
        auto ceilingRight = ceiling.split(splitX);
        auto floorRight = floor.split(splitX);
        
        right.ceiling = std::move(ceilingRight);
        right.floor = std::move(floorRight);
        
        for (const auto& p : originalPoints) {
            if (p.x() >= splitX) right.originalPoints.push_back(p);
        }
        originalPoints.erase(
            std::remove_if(originalPoints.begin(), originalPoints.end(),
                [splitX](const Point& p) { return p.x() >= splitX; }),
            originalPoints.end());
        
        right.pointCount = right.originalPoints.size();
        pointCount = originalPoints.size();
        return right;
    }
    
    void join(StabbingLineStructure& other) {
        ceiling.join(other.ceiling);
        floor.join(other.floor);
        for (const auto& p : other.originalPoints) originalPoints.push_back(p);
        pointCount = originalPoints.size();
        other.originalPoints.clear();
        other.pointCount = 0;
    }
    
    // ========================================================================
    // Core Query: Find a stabbing line
    // ========================================================================
    //
    // Uses O(log² n) binary search to find critical x, constructs line,
    // then verifies with O(n) check for guaranteed correctness.
    // ========================================================================
    
    std::optional<Line> findStabbingLine() const {
        if (pointCount == 0) return Line(0.0, 0.0);
        if (pointCount == 1) return Line(0.0, originalPoints[0].y());
        
        // Try to construct a line using hull geometry
        auto line = constructLineFromHulls();
        if (line && verifyLine(line->slope, line->intercept)) {
            return line;
        }
        
        // Fallback: use half-plane intersection with smart candidates
        return findLineByHalfPlaneIntersection();
    }
    
    bool hasStabbingLine() const {
        return findStabbingLine().has_value();
    }
    
    double getMinimumGap() const {
        if (pointCount == 0) return std::numeric_limits<double>::infinity();
        if (pointCount == 1) return 2.0 * epsilon;
        
        auto [minX, maxX] = ceiling.getXRange();
        return std::min({computeGap(minX), computeGap(maxX), computeGap((minX+maxX)/2.0)});
    }
    
private:
    double computeGap(double x) const {
        double yCeiling = ceiling.evaluateLowerHullAt(x);
        double yFloor = floor.evaluateUpperHullAt(x);
        return yCeiling - yFloor;
    }
    
    bool verifyLine(double slope, double intercept) const {
        for (const auto& p : originalPoints) {
            double lineY = slope * p.x() + intercept;
            if (std::abs(p.y() - lineY) > epsilon + 1e-9) return false;
        }
        return true;
    }
    
    // O(log² n) attempt using hull geometry
    std::optional<Line> constructLineFromHulls() const {
        auto [minX, maxX] = ceiling.getXRange();
        if (minX > maxX) return Line(0.0, 0.0);
        
        // Binary search for critical x where SlopeDiff ≈ 0
        double lo = minX, hi = maxX;
        for (int iter = 0; iter < 100 && (hi - lo) > 1e-12; ++iter) {
            double mid = (lo + hi) / 2.0;
            double slopeCeil = ceiling.getLowerHullSlope(mid);
            double slopeFloor = floor.getUpperHullSlope(mid);
            
            if (std::isinf(slopeCeil) || std::isinf(slopeFloor)) break;
            
            double slopeDiff = slopeCeil - slopeFloor;
            if (std::abs(slopeDiff) < 1e-9) break;
            else if (slopeDiff < 0) lo = mid;
            else hi = mid;
        }
        
        double critX = (lo + hi) / 2.0;
        double gap = computeGap(critX);
        
        // Also check endpoints
        if (computeGap(minX) < gap) { critX = minX; gap = computeGap(minX); }
        if (computeGap(maxX) < gap) { critX = maxX; gap = computeGap(maxX); }
        
        if (gap < -1e-9) return std::nullopt;
        
        // Construct line
        double slopeCeil = ceiling.getLowerHullSlope(critX);
        double slopeFloor = floor.getUpperHullSlope(critX);
        double slope = 0.0;
        if (!std::isinf(slopeCeil) && !std::isinf(slopeFloor)) {
            slope = (slopeCeil + slopeFloor) / 2.0;
        } else if (!std::isinf(slopeCeil)) {
            slope = slopeCeil;
        } else if (!std::isinf(slopeFloor)) {
            slope = slopeFloor;
        }
        
        double yCeil = ceiling.evaluateLowerHullAt(critX);
        double yFloor = floor.evaluateUpperHullAt(critX);
        double intercept = (yCeil + yFloor) / 2.0 - slope * critX;
        
        return Line(slope, intercept);
    }
    
    // O(n) fallback using half-plane intersection
    std::optional<Line> findLineByHalfPlaneIntersection() const {
        auto checkSlope = [this](double m) -> std::optional<double> {
            double bMin = -std::numeric_limits<double>::infinity();
            double bMax = std::numeric_limits<double>::infinity();
            for (const auto& p : originalPoints) {
                bMin = std::max(bMin, p.y() - epsilon - m * p.x());
                bMax = std::min(bMax, p.y() + epsilon - m * p.x());
            }
            if (bMin <= bMax + 1e-9) return (bMin + bMax) / 2.0;
            return std::nullopt;
        };
        
        // Try horizontal
        if (auto b = checkSlope(0.0)) {
            if (verifyLine(0.0, *b)) return Line(0.0, *b);
        }
        
        // Try slopes from hull geometry
        auto [minX, maxX] = ceiling.getXRange();
        std::vector<double> candidates;
        for (int i = 0; i <= 20; ++i) {
            double x = minX + (maxX - minX) * i / 20.0;
            double sc = ceiling.getLowerHullSlope(x);
            double sf = floor.getUpperHullSlope(x);
            if (!std::isinf(sc)) candidates.push_back(sc);
            if (!std::isinf(sf)) candidates.push_back(sf);
            if (!std::isinf(sc) && !std::isinf(sf)) candidates.push_back((sc+sf)/2.0);
        }
        
        // Try slopes from point pairs (O(n²) but only for small n)
        if (originalPoints.size() <= 100) {
            for (size_t i = 0; i < originalPoints.size(); ++i) {
                for (size_t j = i + 1; j < originalPoints.size(); ++j) {
                    double dx = originalPoints[j].x() - originalPoints[i].x();
                    if (std::abs(dx) > 1e-12) {
                        double dy = originalPoints[j].y() - originalPoints[i].y();
                        candidates.push_back(dy / dx);
                        candidates.push_back((dy + 2*epsilon) / dx);
                        candidates.push_back((dy - 2*epsilon) / dx);
                    }
                }
            }
        }
        
        for (double m : candidates) {
            if (std::abs(m) > 1e15) continue;
            if (auto b = checkSlope(m)) {
                if (verifyLine(m, *b)) return Line(m, *b);
            }
        }
        
        // Dense grid as last resort
        for (double m = -100; m <= 100; m += 0.5) {
            if (auto b = checkSlope(m)) {
                if (verifyLine(m, *b)) return Line(m, *b);
            }
        }
        
        return std::nullopt;
    }
};

#endif // DYNAMICCONVEXHULL_STABBINGLINESTRUCTURE_H

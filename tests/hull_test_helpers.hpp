// Hull test helpers - monotone chain algorithm and comparison utilities
#ifndef HULL_TEST_HELPERS_HPP
#define HULL_TEST_HELPERS_HPP

#include <vector>
#include <algorithm>
#include <set>
#include <cmath>

// ============================================================================
// Andrew's Monotone Chain Algorithm - directly computes upper and lower hulls
// Reference: https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
// ============================================================================

namespace monotone_chain {

using Point = std::pair<int, int>;

// Cross product of vectors OA and OB (O is origin)
// Returns positive if counter-clockwise, negative if clockwise, 0 if collinear
inline long long cross(const Point& O, const Point& A, const Point& B) {
    return (long long)(A.first - O.first) * (B.second - O.second) - 
           (long long)(A.second - O.second) * (B.first - O.first);
}

// Build lower hull - points from left to right along the bottom
inline std::vector<Point> lowerHull(std::vector<Point> points) {
    std::sort(points.begin(), points.end());
    
    std::vector<Point> lower;
    for (const auto& p : points) {
        while (lower.size() >= 2 && cross(lower[lower.size()-2], lower[lower.size()-1], p) <= 0) {
            lower.pop_back();
        }
        lower.push_back(p);
    }
    return lower;
}

// Build upper hull - points from left to right along the top
inline std::vector<Point> upperHull(std::vector<Point> points) {
    std::sort(points.begin(), points.end());
    
    std::vector<Point> upper;
    for (const auto& p : points) {
        while (upper.size() >= 2 && cross(upper[upper.size()-2], upper[upper.size()-1], p) >= 0) {
            upper.pop_back();
        }
        upper.push_back(p);
    }
    return upper;
}

} // namespace monotone_chain

// ============================================================================
// Hull comparison helpers
// ============================================================================

namespace hull_helpers {

using Point = std::pair<int, int>;

// Convert CGAL points to integer pairs
template<typename PointType>
inline std::vector<Point> toIntPairs(const std::vector<PointType>& pts) {
    std::vector<Point> result;
    for (const auto& p : pts) {
        result.push_back({(int)std::round(p.x()), (int)std::round(p.y())});
    }
    return result;
}

// Compare two hulls exactly
inline bool hullsEqual(const std::vector<Point>& a, const std::vector<Point>& b) {
    if (a.size() != b.size()) return false;
    auto sorted_a = a;
    auto sorted_b = b;
    std::sort(sorted_a.begin(), sorted_a.end());
    std::sort(sorted_b.begin(), sorted_b.end());
    return sorted_a == sorted_b;
}

// Check if all points in 'expected' are contained in 'actual'
// CHTree may include additional collinear points on hull boundaries
inline bool hullContainsAll(const std::vector<Point>& expected, 
                            const std::vector<Point>& actual) {
    std::set<Point> actual_set(actual.begin(), actual.end());
    for (const auto& p : expected) {
        if (actual_set.find(p) == actual_set.end()) {
            return false;
        }
    }
    return true;
}

// CHTree semantics adjustments:
// - Upper hull excludes leftmost point (shared with lower hull)
// - Lower hull excludes rightmost point (shared with upper hull)

inline std::vector<Point> adjustUpperHullForCHTree(const std::vector<Point>& hull) {
    if (hull.size() <= 1) return hull;
    auto sorted = hull;
    std::sort(sorted.begin(), sorted.end());
    return std::vector<Point>(sorted.begin() + 1, sorted.end());
}

inline std::vector<Point> adjustLowerHullForCHTree(const std::vector<Point>& hull) {
    if (hull.size() <= 1) return hull;
    auto sorted = hull;
    std::sort(sorted.begin(), sorted.end());
    sorted.pop_back();
    return sorted;
}

// ============================================================================
// Tree Validation Helper
// ============================================================================

template<typename Node>
bool validateTreeInvariants(Node* node) {
    if (!node) return true;
    
    // Check leaf vs internal
    bool is_leaf = !node->left && !node->right;
    
    if (!is_leaf) {
        // Internal node must have BOTH children (Leaf-oriented AVL invariant)
        if (!node->left || !node->right) {
            std::cerr << "Validation Error: Internal node missing child\n";
            return false;
        }
        
        // Check sorting invariant: left max < right min
        if (node->left->max[0].max().x() > node->right->min[0].min().x()) {
             std::cerr << "Validation Error: Sorting violation. Left max > Right min\n";
             return false;
        }
        
        // Recurse
        if (!validateTreeInvariants(node->left)) return false;
        if (!validateTreeInvariants(node->right)) return false;
        
        // Check Bridge Validity (Stale Bridge Check)
        // The bridge stored in node->val must have endpoints in the subtrees
        auto bridge = node->val[0]; // Lower bridge
        auto p1 = bridge.min();
        auto p2 = bridge.max();
        
        // Note: covers() checks if point x is in range [min, max] of subtree.
        // It doesn't guarantee the point exists, but it's a good first check.
        // Better: check that bridge points are consistent with children ranges.
        // p1 must be <= left->max and >= left->min
        // p2 must be <= right->max and >= right->min
        
        auto l_min = node->left->min[0].min();
        auto l_max = node->left->max[0].max();
        auto r_min = node->right->min[0].min();
        auto r_max = node->right->max[0].max();
        
        if (p1.x() < l_min.x() || p1.x() > l_max.x()) {
            std::cerr << "Validation Error: Bridge start (" << p1.x() << ") not in left subtree range [" << l_min.x() << ", " << l_max.x() << "]\n";
            return false;
        }
        
        if (p2.x() < r_min.x() || p2.x() > r_max.x()) {
            std::cerr << "Validation Error: Bridge end (" << p2.x() << ") not in right subtree range [" << r_min.x() << ", " << r_max.x() << "]\n";
            return false;
        }
    }
    
    return true;
}

} // namespace hull_helpers

#endif // HULL_TEST_HELPERS_HPP

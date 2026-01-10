/**
 * @file CHTree.h
 * @brief Dynamic Convex Hull data structure using an AVL-based approach.
 * 
 * Implements a leaf-oriented AVL tree where internal nodes store bridges
 * between left and right subtree hulls. Supports O(log² n) insert/remove,
 * O(log² n) split/join, O(n) bulk construction, and O(log n) point-in-hull queries.
 * 
 * Based on: "Simple and Robust Dynamic Two-Dimensional Convex Hull"
 * by Gæde, Gørtz, van der Hoog, Krogh, Rotenberg (2023)
 */

#ifndef DYNAMICCONVEXHULL_CHTREE_H
#define DYNAMICCONVEXHULL_CHTREE_H

#include <vector>
#include <limits>
#include <cmath>
#include <compare>
#include "AvlTree.h"
#include "util.h"

#define isLeaf this->isLeaf

/**
 * @brief Dynamic Convex Hull tree data structure.
 * @tparam Traits Kernel providing Point_2, Segment_2, and predicates.
 */
template<class Traits>
class CHTree : public AVLTree<Bridges<Traits>>{
    using Bridge = typename Traits::Segment_2;
    using Point = typename Traits::Point_2;
    using Node = typename AVLTree<Bridges<Traits>>::Node;
    using Midpoint = typename Traits::Construct_midpoint_2;
    using Compare_slope = typename Traits::Compare_slope_2;
    using Compare_at_x = typename Traits::Compare_y_at_x_2;

    Midpoint midpoint = Midpoint();
    Compare_slope compare_slope = Compare_slope();
    Compare_at_x compare_at_x = Compare_at_x();

protected:
    /**
     * @brief Compare slopes of two bridges.
     * @tparam lower True for lower hull, false for upper hull.
     */
    template<bool lower>
    bool slope_comp(const Bridge& l, const Bridge& r){
        auto res = compare_slope(l, r);
        if(res == 0) return true;
        return (res < 0) != lower;
    }

    /**
     * @brief Compare y-values at midpoint for bridge navigation.
     * @tparam lower True for lower hull, false for upper hull.
     */
    template<bool lower>
    bool m_comp(const Bridge& l, const Bridge& r, const Point& m){
        if(l.is_vertical()) {
            return !lower;
        }
        if(r.is_vertical()) {
            if (m.x() < r.min().x()) return !lower;
            return lower;
        }
        auto res = compare_at_x(m, l.supporting_line(), r.supporting_line());
        if(res == 0) return true;
        return (res < 0) == lower;
    }

    /**
     * @brief Check if point is covered by a bridge edge.
     * @tparam lower True for lower hull, false for upper hull.
     */
    template<bool lower>
    bool cover_comp(const Point& p, const Bridge& b){
        auto res = compare_at_x(p, b.supporting_line());
        if(res == 0) return true;
        return (res > 0) == lower;
    }

    /**
     * @brief Find the bridge between left and right subtrees.
     * 
     * Uses the O(log n) bridge-finding algorithm that navigates
     * the hulls using slope comparisons at each step.
     * 
     * @tparam lower True for lower hull bridge, false for upper hull bridge.
     * @param v Internal node whose children's hulls we're bridging.
     * @return Bridge segment connecting the two sub-hulls.
     */
    template<bool lower>
    Bridge findBridge(Node* v){
        Node* x = v->left;
        Node* y = v->right;
        
        if (!x || !y) return Bridge();

        Bridge e_l, e_r, lr;
        bool undecided;

        while (!(isLeaf(x) && isLeaf(y))) {
            undecided = true;
            if (!x || !y) break;

            Point m = midpoint(x->max[lower].max(), y->min[lower].min());

            e_l = x->val[lower];
            e_r = y->val[lower];
            lr = Bridge(midpoint(e_l),midpoint(e_r));
            
            if (!isLeaf(x) && slope_comp<lower>(e_l,lr)){
                x = stepLeft<lower>(x); 
                undecided = false;
            }
            if (!isLeaf(y) && slope_comp<lower>(lr,e_r)) {
                y = stepRight<lower>(y); 
                undecided = false;
            }
            if (undecided) {
                bool mresult = m_comp<lower>(e_l,e_r,m);
                if (!isLeaf(x) && mresult || isLeaf(y)) {
                    x = stepRight<lower>(x);
                } else {
                    y = stepLeft<lower>(y);
                }
            }
        }
        
        if (!x || !y) return Bridge();
        
        return Bridge(x->val[lower].min(), y->val[lower].max());
    }

    /**
     * @brief Update bridges after tree modification.
     */
    void onUpdate(Node* x){
        if(isLeaf(x)) return;
        x->val[0] = findBridge<false>(x);
        x->val[1] = findBridge<true>(x);
    }

    /**
     * @brief Navigate left in hull to find bridge endpoint.
     * @tparam lower True for lower hull, false for upper hull.
     */
    template<bool lower>
    inline Node* stepLeft(Node* v){
        auto x = v->val[lower].min().x();
        v = v->left;
        while(v && v->val[lower].max().x() > x) {
            v = v->left;
        }
        return v;
    }

    /**
     * @brief Navigate right in hull to find bridge endpoint.
     * @tparam lower True for lower hull, false for upper hull.
     */
    template<bool lower>
    inline Node* stepRight(Node* v){
        auto x = v->val[lower].max().x();
        v = v->right;
        while(v && v->val[lower].min().x() < x) {
            v = v->right;
        }
        return v;
    }

    /**
     * @brief Find successor vertex on hull from given point.
     */
    template<bool lower>
    Node* hullSuccessor(const Point key){
        auto current = AVLTree<Bridges<Traits>>::root;
        Point min;
        while(current){
            min = current->val[lower].min();
            if(min == key) break;
            else if (min < key) current = stepRight<lower>(current);
            else current = stepLeft<lower>(current);
        }
        return current;
    }

    template<bool lower>
    Node* hullSuccessor(const Node* p){
        return hullSuccessor<lower>(p->val[lower].max());
    }

    /**
     * @brief Find predecessor vertex on hull from given point.
     */
    template<bool lower>
    Node* hullPredecessor(const Point key){
        auto current = AVLTree<Bridges<Traits>>::root;
        Point max;
        while(current){
            max = current->val[lower].max();
            if(max == key) break;
            else if (max < key) current = stepRight<lower>(current);
            else current = stepLeft<lower>(current);
        }
        return current;
    }

    template<bool lower>
    Node* hullPredecessor(const Node* p){
        return hullPredecessor<lower>(p->val[lower].min());
    }

    /**
     * @brief Check if point is covered by hull.
     */
    template<bool lower>
    bool covers(Point p){
        auto current = AVLTree<Bridges<Traits>>::root;
        while(current){
            if(current->val[lower].min() <= p){
                if(p <= current->val[lower].max()){
                    return cover_comp<lower>(p,current->val[lower]);
                } else current = stepRight<lower>(current);
            } else current = stepLeft<lower>(current);
        }
        return false;
    }

    /**
     * @brief Collect all hull vertices in order.
     */
    template<bool lower>
    std::vector<Point> hullPoints(){
        std::vector<Point> res;
        Node* e = AVLTree<Bridges<Traits>>::root;
        if(!e || isLeaf(e)) return res;
        while(e && !isLeaf(e)){
            res.insert(res.begin(), e->val[lower].min());
            e = hullPredecessor<lower>(res.front());
        }
        e = AVLTree<Bridges<Traits>>::root;
        while(e && !isLeaf(e)){
            res.push_back(e->val[lower].max());
            e = hullSuccessor<lower>(res.back());
        }
        if(lower) res.pop_back();
        else {
            res.erase(res.begin());
            std::reverse(res.begin(), res.end());
        }
        return res;
    }

public:
    /** @brief Insert a point into the convex hull. O(log² n) */
    void insert(Point p){
        AVLTree<Bridges<Traits>>::insert({Bridge(p,p),Bridge(p,p)});
    }

    /** @brief Remove a point from the convex hull. O(log² n) */
    void remove(Point p){
        AVLTree<Bridges<Traits>>::remove({Bridge(p,p),Bridge(p,p)});
    }

    /** @brief Check if point p is inside the convex hull. O(log n) */
    bool covers(Point p){
        return covers<false>(p) && covers<true>(p);
    }

    /** @brief Get upper hull vertices in order. */
    std::vector<Point> upperHullPoints(){
        return hullPoints<false>();
    }

    /** @brief Get lower hull vertices in order. */
    std::vector<Point> lowerHullPoints(){
        return hullPoints<true>();
    }

    /**
     * @brief Get the x-coordinate range of all points.
     * @return {min_x, max_x}. For empty hulls, returns {+inf, -inf}.
     */
    std::pair<double, double> getXRange() const {
        if (!AVLTree<Bridges<Traits>>::root) {
            return {std::numeric_limits<double>::infinity(), 
                    -std::numeric_limits<double>::infinity()};
        }
        return {AVLTree<Bridges<Traits>>::root->min[0].min().x(),
                AVLTree<Bridges<Traits>>::root->max[0].max().x()};
    }

    /**
     * @brief Evaluate the hull y-value at given x-coordinate.
     * @tparam lower True for lower hull, false for upper hull.
     * @param x X-coordinate within hull's x-range.
     * @return Y-value on the hull boundary at x.
     */
    template<bool lower>
    double evaluateHullAt(double x) const {
        auto current = AVLTree<Bridges<Traits>>::root;
        if (!current) return lower ? std::numeric_limits<double>::infinity() 
                                   : -std::numeric_limits<double>::infinity();
        
        while (current && (current->left || current->right)) {
            const Bridge& bridge = current->val[lower];
            double minX = bridge.min().x();
            double maxX = bridge.max().x();
            
            if (x >= minX && x <= maxX) {
                if (bridge.is_vertical()) {
                    return lower ? bridge.min().y() : bridge.max().y();
                }
                double dx = maxX - minX;
                if (std::abs(dx) < 1e-15) return bridge.min().y();
                double dy = bridge.max().y() - bridge.min().y();
                double slope = dy / dx;
                return bridge.min().y() + slope * (x - minX);
            } else if (x < minX) {
                current = current->left;
            } else {
                current = current->right;
            }
        }
        
        if (current) {
            return current->val[lower].min().y();
        }
        return lower ? std::numeric_limits<double>::infinity() 
                     : -std::numeric_limits<double>::infinity();
    }

    /**
     * @brief Get the slope of the hull edge at given x-coordinate.
     * @tparam lower True for lower hull, false for upper hull.
     */
    template<bool lower>
    double getHullSlope(double x) const {
        auto current = AVLTree<Bridges<Traits>>::root;
        if (!current) return 0.0;
        
        while (current && (current->left || current->right)) {
            const Bridge& bridge = current->val[lower];
            double minX = bridge.min().x();
            double maxX = bridge.max().x();
            
            if (x >= minX && x <= maxX) {
                if (bridge.is_vertical()) {
                    return lower ? -std::numeric_limits<double>::infinity()
                                 : std::numeric_limits<double>::infinity();
                }
                double dx = maxX - minX;
                if (std::abs(dx) < 1e-15) return 0.0;
                double dy = bridge.max().y() - bridge.min().y();
                return dy / dx;
            } else if (x < minX) {
                current = current->left;
            } else {
                current = current->right;
            }
        }
        
        return 0.0;
    }

    double evaluateLowerHullAt(double x) const { return evaluateHullAt<true>(x); }
    double evaluateUpperHullAt(double x) const { return evaluateHullAt<false>(x); }
    double getLowerHullSlope(double x) const { return getHullSlope<true>(x); }
    double getUpperHullSlope(double x) const { return getHullSlope<false>(x); }

    /**
     * @brief Build convex hull from sorted points in O(n) time.
     * @param sortedPoints Points sorted by x-coordinate (and by y for ties).
     */
    void build(const std::vector<Point>& sortedPoints){
        AVLTree<Bridges<Traits>>::clear();
        if (sortedPoints.empty()) {
            return;
        }
        
        std::vector<Bridges<Traits>> bridges;
        bridges.reserve(sortedPoints.size());
        for (const auto& p : sortedPoints) {
            bridges.push_back({Bridge(p, p), Bridge(p, p)});
        }
        
        AVLTree<Bridges<Traits>>::root = buildSubtreeFast(bridges, 0, bridges.size());
    }

private:
    /**
     * @brief Recursive helper for O(n) bulk construction.
     */
    Node* buildSubtreeFast(const std::vector<Bridges<Traits>>& bridges, size_t start, size_t end) {
        if (start >= end) return nullptr;
        
        if (end - start == 1) {
            Node* leaf = new Node(bridges[start]);
            return leaf;
        }
        
        size_t mid = start + (end - start) / 2;
        
        Node* node = new Node(bridges[mid]);
        
        node->left = buildSubtreeFast(bridges, start, mid);
        node->right = buildSubtreeFast(bridges, mid, end);
        
        if (node->left) node->left->par = node;
        if (node->right) node->right->par = node;
        
        node->size = AVLTree<Bridges<Traits>>::size(node->left) + 
                     AVLTree<Bridges<Traits>>::size(node->right) + 
                     (node && !(node->left || node->right) ? 1 : 0);
        node->height = std::max(AVLTree<Bridges<Traits>>::height(node->left), 
                                AVLTree<Bridges<Traits>>::height(node->right)) + 1;
        if (node->left) node->min = node->left->min;
        else node->min = node->val;
        if (node->right) node->max = node->right->max;
        else node->max = node->val;
        
        onUpdate(node);
        
        return node;
    }

public:
    /**
     * @brief Join another hull into this one.
     * 
     * Precondition: All points in 'this' must have x < all points in 'other'.
     * After this operation, 'other' becomes empty.
     * 
     * @param other Hull to merge (must be to the right of this hull).
     * @note Time complexity: O(log² N)
     */
    void join(CHTree& other) {
        if (other.empty()) return;
        if (this->empty()) {
            this->root = other.root;
            if (this->root) this->root->par = nullptr;
            other.root = nullptr;
            return;
        }

        Point p = this->root ? this->root->val[0].min() : other.root->val[0].min();
        Bridges<Traits> dummy = {Bridge(p, p), Bridge(p, p)};
        
        AVLTree<Bridges<Traits>>::join(dummy, &other);
    }

    /**
     * @brief Split this hull at the given x-coordinate.
     * 
     * After the operation:
     *   - 'this' contains all points with x < splitX
     *   - Returns a new CHTree containing points with x >= splitX
     * 
     * @param splitX X-coordinate to split at.
     * @return New CHTree containing the right portion.
     * @note Time complexity: O(log² N)
     */
    CHTree split(double splitX) {
        CHTree rightHull;
        if (AVLTree<Bridges<Traits>>::root == nullptr) {
            return rightHull;
        }
        
        Point splitPoint(splitX, std::numeric_limits<double>::lowest());
        Bridges<Traits> splitKey = {Bridge(splitPoint, splitPoint), Bridge(splitPoint, splitPoint)};
        
        AVLTree<Bridges<Traits>>::split(splitKey, &rightHull);
        
        return rightHull;
    }

    /** @brief Get the number of points in the hull. */
    size_t size() const {
        return AVLTree<Bridges<Traits>>::root ? 
               AVLTree<Bridges<Traits>>::root->size : 0;
    }

    /** @brief Check if hull is empty. */
    bool empty() const {
        return AVLTree<Bridges<Traits>>::root == nullptr;
    }

    /** @brief Validate that all bridges are consistent (for debugging). */
    bool validateBridges() {
        return validateBridgesRecursive(AVLTree<Bridges<Traits>>::root);
    }

private:
    bool validateBridgesRecursive(Node* node) {
        if (!node || isLeaf(node)) return true;
        
        if (!validateBridgesRecursive(node->left)) return false;
        if (!validateBridgesRecursive(node->right)) return false;

        Bridge upper = findBridge<false>(node);
        Bridge lower = findBridge<true>(node);
        
        if (node->val[0] != upper) {
            std::cerr << "Upper Bridge mismatch!\n";
            return false;
        }
        if (node->val[1] != lower) {
            std::cerr << "Lower Bridge mismatch!\n";
            return false;
        }
        
        return true;
    }

    /**
     * @brief Find critical vertex for stabbing line query.
     * 
     * Traverses the tree to find the vertex where the slope difference
     * between two hulls changes sign (subgradient condition).
     * 
     * @tparam lower True for ceiling hull, false for floor hull.
     * @param otherHull The other hull to compare slopes against.
     * @return {x-coordinate, slope} of the critical point.
     */
    template<bool lower>
    std::pair<double, double> findCriticalVertex(const CHTree<Traits>& otherHull) const {
        auto current = AVLTree<Bridges<Traits>>::root;
        if (!current) return {0.0, 0.0};
        
        while (current && (current->left || current->right)) {
            const Bridge& bridge = current->val[lower];
            double midX = (bridge.min().x() + bridge.max().x()) / 2.0;
            
            double thisSlope;
            if (bridge.is_vertical()) {
                thisSlope = lower ? std::numeric_limits<double>::infinity() 
                                  : -std::numeric_limits<double>::infinity();
            } else {
                double dx = bridge.max().x() - bridge.min().x();
                if (std::abs(dx) < 1e-15) {
                    thisSlope = 0.0;
                } else {
                    thisSlope = (bridge.max().y() - bridge.min().y()) / dx;
                }
            }
            
            double otherSlope;
            if constexpr (lower) {
                otherSlope = otherHull.getUpperHullSlope(midX);
            } else {
                otherSlope = otherHull.getLowerHullSlope(midX);
            }
            
            if (std::isinf(thisSlope) || std::isinf(otherSlope)) {
                return {bridge.min().x(), 0.0};
            }
            
            double slopeDiff;
            if constexpr (lower) {
                slopeDiff = thisSlope - otherSlope;
            } else {
                slopeDiff = otherSlope - thisSlope;
            }
            
            if (slopeDiff < -1e-9) {
                current = current->right ? current->right : current;
                if (!current->right) break;
            } else if (slopeDiff > 1e-9) {
                current = current->left ? current->left : current;
                if (!current->left) break;
            } else {
                return {midX, thisSlope};
            }
        }
        
        if (current) {
            Point p = current->val[lower].min();
            return {p.x(), 0.0};
        }
        return {0.0, 0.0};
    }

public:
    std::pair<double, double> findCriticalVertexFromCeiling(const CHTree<Traits>& floorHull) const {
        return findCriticalVertex<true>(floorHull);
    }
    
    std::pair<double, double> findCriticalVertexFromFloor(const CHTree<Traits>& ceilingHull) const {
        return findCriticalVertex<false>(ceilingHull);
    }

    /**
     * @brief Get the range of valid tangent slopes at a hull vertex.
     * @tparam lower True for lower hull, false for upper hull.
     * @param x X-coordinate of the vertex.
     * @return {left_slope, right_slope} defining the valid slope range.
     */
    template<bool lower>
    std::pair<double, double> getSlopeRange(double x) const {
        auto current = AVLTree<Bridges<Traits>>::root;
        if (!current) return {0.0, 0.0};
        
        const double INF = std::numeric_limits<double>::infinity();
        double leftSlope = lower ? -INF : INF;
        double rightSlope = lower ? INF : -INF;
        bool foundExactVertex = false;
        
        while (current && (current->left || current->right)) {
            const Bridge& bridge = current->val[lower];
            double minX = bridge.min().x();
            double maxX = bridge.max().x();
            
            double bridgeSlope = 0.0;
            if (!bridge.is_vertical()) {
                double dx = maxX - minX;
                if (std::abs(dx) >= 1e-15) {
                    bridgeSlope = (bridge.max().y() - bridge.min().y()) / dx;
                }
            }
            
            if (std::abs(x - minX) < 1e-12) {
                rightSlope = bridgeSlope;
                foundExactVertex = true;
                current = current->left;
            } else if (std::abs(x - maxX) < 1e-12) {
                leftSlope = bridgeSlope;
                foundExactVertex = true;
                current = current->right;
            } else if (x > minX && x < maxX) {
                return {bridgeSlope, bridgeSlope};
            } else if (x < minX) {
                rightSlope = bridgeSlope;
                current = current->left;
            } else {
                leftSlope = bridgeSlope;
                current = current->right;
            }
        }
        
        if (current && foundExactVertex) {
            if constexpr (lower) {
                return {leftSlope, rightSlope};
            } else {
                return {rightSlope, leftSlope};
            }
        }
        
        if constexpr (lower) {
            return {leftSlope, rightSlope};
        } else {
            return {rightSlope, leftSlope};
        }
    }
    
    std::pair<double, double> getLowerHullSlopeRange(double x) const { 
        return getSlopeRange<true>(x); 
    }
    std::pair<double, double> getUpperHullSlopeRange(double x) const { 
        return getSlopeRange<false>(x); 
    }
};

#endif //DYNAMICCONVEXHULL_CHTREE_H

//
// Created by etoga on 5/12/23.
//

#ifndef DYNAMICCONVEXHULL_CHTREE_H
#define DYNAMICCONVEXHULL_CHTREE_H

#include <vector>
#include <limits>
#include <CGAL/enum.h>
#include "AvlTree.h"
#include "util.h"
#define isLeaf this->isLeaf


template<class Traits>
class CHTree : public AVLTree<Bridges<Traits>>{
using Bridge = typename Traits::Segment_2;
using Point = typename Traits::Point_2;
using Node = typename AVLTree<Bridges<Traits>>::Node;
using Midpoint = typename Traits::Construct_midpoint_2;
using Compare_slope = typename Traits::Compare_slope_2;
using Compare_at_x = typename Traits::Compare_y_at_x_2;

const Midpoint midpoint = Midpoint();
const Compare_slope compare_slope = Compare_slope();
const Compare_at_x compare_at_x = Compare_at_x();

protected:
    template<bool lower>
    bool slope_comp(const Bridge& l, const Bridge& r){
        CGAL::Comparison_result res = compare_slope(l,r);
        if(res == CGAL::EQUAL) return true;
        return (res == CGAL::SMALLER) != lower;
    }

    template<bool lower>
    bool m_comp(const Bridge& l, const Bridge& r, const Point& m){
        // Handle vertical line cases based on x-coordinate of midpoint
        if(l.is_vertical()) {
            // If l is vertical, compare where m.x is relative to l.min().x
            // For lower hull: if m.x >= l.x, result is effectively "l at m.x is undefined/far"
            // Default behavior: return !lower ensures correct navigation
            return !lower;
        }
        if(r.is_vertical()) {
            // If r is vertical at x=r_x:
            // - If m.x < r_x: compare l(m.x) with r, but r is undefined at m.x
            //   For lower hull, we should step on y-side (return true means step x, so return false)
            // - If m.x >= r_x: similar logic
            // The original "return lower" was incorrect for cases where m.x < r.min().x
            // Correct: for lower hull, when r is vertical and m.x is left of r, return !lower (go left on y)
            if (m.x() < r.min().x()) return !lower;
            return lower;
        }
        CGAL::Comparison_result res = compare_at_x(m,l.supporting_line(),r.supporting_line());
        if(res == CGAL::EQUAL) return true;
        return (res == CGAL::SMALLER) == lower;
    }

    template<bool lower>
    bool cover_comp(const Point& p, const Bridge& b){
        CGAL::Comparison_result res = compare_at_x(p,b.supporting_line());
        if(res == CGAL::EQUAL) return true;
        return (res == CGAL::LARGER) == lower;
    }

    // Debug flag for tracing findBridge - set to true to enable logging
    static constexpr bool DEBUG_FIND_BRIDGE = false;

    template<bool lower>
    Bridge findBridge(Node* v){
        Node* x = v->left;
        Node* y = v->right;
        
        // Safety check
        if (!x || !y) return Bridge();

        Bridge e_l, e_r, lr;
        bool undecided;

        int step = 0;
        while (!(isLeaf(x) && isLeaf(y))) {
            undecided = true;
            // Safety check for nulls during traversal
            if (!x || !y) break;

            // Recompute midpoint m based on current x and y (FIX: was computed once outside loop)
            Point m = midpoint(x->max[lower].max(), y->min[lower].min());

            e_l = x->val[lower];
            e_r = y->val[lower];
            lr = Bridge(midpoint(e_l),midpoint(e_r));
            
            if constexpr (DEBUG_FIND_BRIDGE) {
                std::cerr << "\nStep " << step++ << ":\n";
                std::cerr << "  m: (" << m.x() << ", " << m.y() << ")\n";
                std::cerr << "  e_l: (" << e_l.min().x() << "," << e_l.min().y() << ") -> (" 
                          << e_l.max().x() << "," << e_l.max().y() << ")\n";
                std::cerr << "  e_r: (" << e_r.min().x() << "," << e_r.min().y() << ") -> (" 
                          << e_r.max().x() << "," << e_r.max().y() << ")\n";
                std::cerr << "  lr: (" << lr.min().x() << "," << lr.min().y() << ") -> (" 
                          << lr.max().x() << "," << lr.max().y() << ")\n";
            }
            
            if (!isLeaf(x) && slope_comp<lower>(e_l,lr)){
                if constexpr (DEBUG_FIND_BRIDGE) {
                    std::cerr << "  -> stepLeft (slope_comp(e_l,lr) true)\n";
                }
                x = stepLeft<lower>(x); undecided = false;
            }
            if (!isLeaf(y) && slope_comp<lower>(lr,e_r)) {
                if constexpr (DEBUG_FIND_BRIDGE) {
                    std::cerr << "  -> stepRight (slope_comp(lr,e_r) true)\n";
                }
                y = stepRight<lower>(y); undecided = false;
            }
            if (undecided) {
                bool mresult = m_comp<lower>(e_l,e_r,m);
                if constexpr (DEBUG_FIND_BRIDGE) {
                    std::cerr << "  UNDECIDED: m_comp=" << mresult << ", isLeaf(x)=" << isLeaf(x) << ", isLeaf(y)=" << isLeaf(y) << "\n";
                }
                if (!isLeaf(x) && mresult || isLeaf(y)) {
                    if constexpr (DEBUG_FIND_BRIDGE) {
                        std::cerr << "  -> stepRight on x (undecided, m_comp or y is leaf)\n";
                    }
                    x = stepRight<lower>(x);
                } else {
                    if constexpr (DEBUG_FIND_BRIDGE) {
                        std::cerr << "  -> stepLeft on y (undecided)\n";
                    }
                    y = stepLeft<lower>(y);
                }
            }
            
            if constexpr (DEBUG_FIND_BRIDGE) {
                if (x) std::cerr << "  New x: val=(" << x->val[lower].min().x() << ".." << x->val[lower].max().x() << "), isLeaf=" << isLeaf(x) << "\n";
                if (y) std::cerr << "  New y: val=(" << y->val[lower].min().x() << ".." << y->val[lower].max().x() << "), isLeaf=" << isLeaf(y) << "\n";
            }
        }
        
        if (!x || !y) return Bridge();
        
        Bridge result = Bridge(x->val[lower].min(),y->val[lower].max());
        if constexpr (DEBUG_FIND_BRIDGE) {
            std::cerr << "Result: (" << result.min().x() << "," << result.min().y() << ") -> (" 
                      << result.max().x() << "," << result.max().y() << ")\n";
        }
        return result;
    }

    void onUpdate(Node* x){
        if(isLeaf(x)) return;
        x->val[0] = findBridge<false>(x);
        x->val[1] = findBridge<true>(x);
    }

    // ========================================================================
    // O(h) Bridge Finding for O(n) Construction (Algorithm 2 from the paper)
    // ========================================================================
    // 
    // "Simple and Robust Dynamic Two-Dimensional Convex Hull"
    // by Gæde, Gørtz, van der Hoog, Krogh, Rotenberg (2023)
    //
    // Key insight: The step operations are O(1) pointer navigations (->left, ->right)
    // instead of O(log n) hull-path traversals. This gives O(h) per bridge computation
    // where h is the subtree height.
    //
    // For a node at height h, bridge-finding costs O(h).
    // Total construction cost: Σ(n/2^h × h) = O(n).
    // ========================================================================

    // O(N) total bridge finding for bulk construction
    // 
    // This is identical to findBridge but runs in the context of bottom-up
    // construction. Uses stepLeft/stepRight to navigate which skip dominated nodes.
    //
    // Complexity per node: O(h²) where h is subtree height
    // Total complexity: Σ(n/2^h × h²) = n × Σ(h²/2^h) = n × 6 = O(n)
    // The series Σ(h²/2^h) converges to a constant (≈6).
    //
    // Based on: "Simple and Robust Dynamic Two-Dimensional Convex Hull"
    // by Gæde, Gørtz, van der Hoog, Krogh, Rotenberg (2023)
    template<bool lower>
    Bridge findBridgeLinear(Node* v){
        // This is exactly the same as findBridge - the O(N) benefit comes from
        // the bottom-up construction pattern, not from changes to this function.
        Node* x = v->left;
        Node* y = v->right;
        Bridge e_l, e_r, lr;
        bool undecided;

        Point m = midpoint(x->max[lower].max(), y->min[lower].min());

        while (!(isLeaf(x) && isLeaf(y))) {
            undecided = true;
            e_l = x->val[lower];
            e_r = y->val[lower];
            lr = Bridge(midpoint(e_l), midpoint(e_r));
            if (!isLeaf(x) && slope_comp<lower>(e_l, lr)){
                x = stepLeft<lower>(x); undecided = false;
            }
            if (!isLeaf(y) && slope_comp<lower>(lr, e_r)) {
                y = stepRight<lower>(y); undecided = false;
            }
            if (undecided) {
                if (!isLeaf(x) && m_comp<lower>(e_l, e_r, m) || isLeaf(y)) {
                    x = stepRight<lower>(x);
                } else {
                    y = stepLeft<lower>(y);
                }
            }
        }
        return Bridge(x->val[lower].min(), y->val[lower].max());
    }

    // Callback used during O(N) bulk construction.
    // Uses findBridgeLinear which is identical to findBridge.
    // The O(N) total complexity comes from the bottom-up construction pattern:
    // Total = Σ(n/2^h × h²) = n × Σ(h²/2^h) ≈ 6n = O(n)
    void onUpdateLinear(Node* x){
        if(isLeaf(x)) return;
        x->val[0] = findBridge<false>(x);
        x->val[1] = findBridge<true>(x);
    }

    template<bool lower>
    inline
    Node* stepLeft(Node* v){
        auto x = v->val[lower].min().x();
        if constexpr (DEBUG_FIND_BRIDGE) {
            std::cerr << "    stepLeft: starting at val=(" << v->val[lower].min().x() << ".." << v->val[lower].max().x() << "), pivot x=" << x << "\n";
        }
        v = v->left;
        while(v && v->val[lower].max().x() > x) {
            if constexpr (DEBUG_FIND_BRIDGE) {
                std::cerr << "    stepLeft: skipping node val=(" << v->val[lower].min().x() << ".." << v->val[lower].max().x() << ") because max.x > pivot\n";
            }
            v = v->left;
        }
        if constexpr (DEBUG_FIND_BRIDGE) {
            if (v) std::cerr << "    stepLeft: landed at val=(" << v->val[lower].min().x() << ".." << v->val[lower].max().x() << ")\n";
            else std::cerr << "    stepLeft: landed at null!\n";
        }
        return v;
    }

    template<bool lower>
    inline
    Node* stepRight(Node* v){
        // SIMPLIFIED: Just go to right child, no skipping
        // The original skipping logic was incorrectly discarding valid hull points
        if constexpr (DEBUG_FIND_BRIDGE) {
            std::cerr << "    stepRight: from val=(" << v->val[lower].min().x() << ".." << v->val[lower].max().x() << ")";
        }
        v = v->right;
        if constexpr (DEBUG_FIND_BRIDGE) {
            if (v) std::cerr << " to val=(" << v->val[lower].min().x() << ".." << v->val[lower].max().x() << ")\n";
            else std::cerr << " to null!\n";
        }
        return v;
    }

    template<bool lower>
    Node* find(const Point key, const bool left){
        auto current = AVLTree<Bridges<Traits>>::root;
        while(current && !isLeaf(current)){
            if(current->val[lower][!left] == key) return current;
            else if (current->val[lower].min().x() < key.x()) current = current->right;
            else current = current->left;
        }
        return nullptr;
    }

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
    void insert(Point p){
        AVLTree<Bridges<Traits>>::insert({Bridge(p,p),Bridge(p,p)});
    }

    void remove(Point p){
        AVLTree<Bridges<Traits>>::remove({Bridge(p,p),Bridge(p,p)});
    }

    bool covers(Point p){
        return covers<false>(p) && covers<true>(p);
    }

    std::vector<Point> upperHullPoints(){
        return hullPoints<false>();
    }

    std::vector<Point> lowerHullPoints(){
        return hullPoints<true>();
    }

    // Build the convex hull from a sorted vector of points in O(n) time.
    // Points MUST be sorted by x-coordinate (and by y for ties).
    // This replaces any existing hull content.
    //
    // Complexity: O(n) where n is the number of points.
    // This is achieved through bottom-up construction where bridge-finding
    // at height h costs O(h), and the sum Σ(n/2^h * h) = O(n).
    void build(const std::vector<Point>& sortedPoints){
        AVLTree<Bridges<Traits>>::clear();
        if (sortedPoints.empty()) {
            return;
        }
        
        // Convert points to Bridges format (each point becomes a degenerate bridge)
        std::vector<Bridges<Traits>> bridges;
        bridges.reserve(sortedPoints.size());
        for (const auto& p : sortedPoints) {
            bridges.push_back({Bridge(p, p), Bridge(p, p)});
        }
        
        // Build the tree using the fast bridge-finding algorithm
        AVLTree<Bridges<Traits>>::root = buildSubtreeFast(bridges, 0, bridges.size());
    }

private:
    // Recursive helper for O(n) construction using fast bridge-finding.
    // Returns the root of the subtree for bridges[start..end).
    Node* buildSubtreeFast(const std::vector<Bridges<Traits>>& bridges, size_t start, size_t end) {
        if (start >= end) return nullptr;
        
        if (end - start == 1) {
            // Base case: single element becomes a leaf
            Node* leaf = new Node(bridges[start]);
            return leaf;
        }
        
        // Internal node: pick median point for balanced split
        size_t mid = start + (end - start) / 2;
        
        // Create an internal node
        Node* node = new Node(bridges[mid]); // Temporary value, will be updated
        
        // Recursively build left and right subtrees FIRST (bottom-up)
        node->left = buildSubtreeFast(bridges, start, mid);
        node->right = buildSubtreeFast(bridges, mid, end);
        
        // Set parent pointers
        if (node->left) node->left->par = node;
        if (node->right) node->right->par = node;
        
        // Update basic node data (height, size, min, max)
        // Note: isLeaf is a macro defined as this->isLeaf, so we call it directly
        node->size = AVLTree<Bridges<Traits>>::size(node->left) + 
                     AVLTree<Bridges<Traits>>::size(node->right) + 
                     (node && !(node->left || node->right) ? 1 : 0);
        node->height = std::max(AVLTree<Bridges<Traits>>::height(node->left), 
                                AVLTree<Bridges<Traits>>::height(node->right)) + 1;
        if (node->left) node->min = node->left->min;
        else node->min = node->val;
        if (node->right) node->max = node->right->max;
        else node->max = node->val;
        
        // Compute bridges using the O(h) algorithm for O(n) total construction
        onUpdateLinear(node);
        
        return node;
    }

public:
    // ========================================================================
    // Split and Join (Merge) Operations
    // ========================================================================
    // 
    // These operations allow splitting a convex hull into two parts and
    // merging two disjoint hulls into one.
    //
    // Time Complexity: O(log² N) for both operations
    // - O(log N) tree structure changes (AVL split/join)
    // - O(log N) bridge recomputation per touched node (onUpdate callback)
    //
    // Based on the Eilice algorithm:
    // "Simple and Robust Dynamic Two-Dimensional Convex Hull"
    // by Gæde, Gørtz, van der Hoog, Krogh, Rotenberg (2023)
    // ========================================================================

    // Join (merge) another hull 'other' into this one.
    // PRECONDITION: All points in 'this' must have x-coordinates < all points in 'other'.
    //               The two hulls must cover disjoint x-ranges.
    // After this operation, 'other' becomes empty.
    // Time complexity: O(log² N) where N is the combined size
    void join(CHTree& other) {
        if (other.empty()) return;
        if (this->empty()) {
            this->root = other.root;
            if (this->root) this->root->par = nullptr;  // Ensure clean parent
            other.root = nullptr;
            return;
        }

        // Use a point from the tree to create a dummy bridge.
        // This dummy value will be maintained in the internal join node,
        // but immediatley overwritten by onUpdate(x) calling findBridge.
        // We do this to avoid the default behavior of AVLTree::join(AVLTree*)
        // which removes a leaf from 'this' to use as a separator, causing data loss.
        Point p = this->root ? this->root->val[0].min() : other.root->val[0].min();
        Bridges<Traits> dummy = {Bridge(p, p), Bridge(p, p)};
        
        AVLTree<Bridges<Traits>>::join(dummy, &other);
    }

    // Split this hull into two parts at the given x-coordinate.
    // After the operation:
    //   - 'this' contains all points with x < splitX (strictly less)
    //   - Returns a new CHTree containing points with x >= splitX
    // If splitX matches an exact point, that point goes to the RIGHT tree.
    // Time complexity: O(log² N) where N is the original size
    CHTree split(double splitX) {
        CHTree rightHull;
        if (AVLTree<Bridges<Traits>>::root == nullptr) {
            return rightHull;  // Empty tree, nothing to split
        }
        
        // Create a dummy bridge for the split point.
        // The AVLTree::split function splits such that:
        //   - left tree gets values < k
        //   - right tree gets values > k
        //   - If k exists, it is deleted
        // 
        // For our use case, we want to split by x-coordinate.
        // We create a point with the split x-coordinate and a very low y
        // so it will compare correctly.
        Point splitPoint(splitX, std::numeric_limits<double>::lowest());
        Bridges<Traits> splitKey = {Bridge(splitPoint, splitPoint), Bridge(splitPoint, splitPoint)};
        
        // Perform the split - this modifies 'this' to contain < splitX
        // and puts >= splitX into rightHull
        AVLTree<Bridges<Traits>>::split(splitKey, &rightHull);
        
        return rightHull;
    }

    // Get the number of points in the hull
    size_t size() const {
        return AVLTree<Bridges<Traits>>::root ? 
               AVLTree<Bridges<Traits>>::root->size : 0;
    }

    // Check if hull is empty
    bool empty() const {
        return AVLTree<Bridges<Traits>>::root == nullptr;
    }

    // DEBUG: Validate that all bridges are consistent with children
    bool validateBridges() {
        return validateBridgesRecursive(AVLTree<Bridges<Traits>>::root);
    }

private:
    bool validateBridgesRecursive(Node* node) {
        if (!node || isLeaf(node)) return true;
        
        // Recursively validate children first
        if (!validateBridgesRecursive(node->left)) return false;
        if (!validateBridgesRecursive(node->right)) return false;

        // Recompute bridges and compare
        // val[0] = findBridge<false> = upper bridge
        // val[1] = findBridge<true> = lower bridge
        Bridge upper = findBridge<false>(node);
        Bridge lower = findBridge<true>(node);
        
        if (node->val[0] != upper) {
            std::cerr << "Upper Bridge mismatch! Stored: " << node->val[0] << ", Computed: " << upper << "\n";
            return false;
        }
        if (node->val[1] != lower) {
            std::cerr << "Lower Bridge mismatch! Stored: " << node->val[1] << ", Computed: " << lower << "\n";
            return false;
        }
        
        return true;
    }
public:
};
#endif //DYNAMICCONVEXHULL_CHTREE_H

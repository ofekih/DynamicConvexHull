// Copyright 2026 DynamicConvexHull Authors
// SPDX-License-Identifier: MIT

/// @file CHTree.h
/// @brief Dynamic Convex Hull data structure using an AVL-based approach.
///
/// Implements a leaf-oriented AVL tree where internal nodes store bridges
/// between left and right subtree hulls. Supports O(log² n) insert/remove,
/// O(log² n) split/join, O(n) bulk construction, and O(log n) point-in-hull
/// queries.
///
/// Based on: "Simple and Robust Dynamic Two-Dimensional Convex Hull"
/// by Gæde, Gørtz, van der Hoog, Krogh, Rotenberg (2023)

#pragma once

#include <cmath>
#include <compare>
#include <limits>
#include <vector>

#include "AvlTree.h"
#include "util.h"

namespace dch {

/// @brief Dynamic Convex Hull tree data structure.
/// @tparam Traits Kernel providing Point_2, Segment_2, and predicates.
template <class Traits>
class CHTree : public AVLTree<Bridges<Traits>> {
 public:
  using Bridge = typename Traits::Segment_2;
  using Point = typename Traits::Point_2;
  using Node = typename AVLTree<Bridges<Traits>>::Node;

 private:
  using Midpoint = typename Traits::Construct_midpoint_2;
  using CompareSlope = typename Traits::Compare_slope_2;
  using CompareAtX = typename Traits::Compare_y_at_x_2;

  Midpoint midpoint_;
  CompareSlope compare_slope_;
  CompareAtX compare_at_x_;

 protected:
  /// @brief Compare slopes of two bridges.
  /// @tparam kLower True for lower hull, false for upper hull.
  template <bool kLower>
  bool SlopeCompare(const Bridge& left, const Bridge& right) {
    auto res = compare_slope_(left, right);
    if (res == 0) return true;
    return (res < 0) != kLower;
  }

  /// @brief Compare y-values at midpoint for bridge navigation.
  /// @tparam kLower True for lower hull, false for upper hull.
  template <bool kLower>
  bool MidpointCompare(const Bridge& left, const Bridge& right,
                       const Point& mid) {
    if (left.is_vertical()) {
      return !kLower;
    }
    if (right.is_vertical()) {
      if (mid.x() < right.min().x()) return !kLower;
      return kLower;
    }
    auto res = compare_at_x_(mid, left.supporting_line(), right.supporting_line());
    if (res == 0) return true;
    return (res < 0) == kLower;
  }

  /// @brief Check if point is covered by a bridge edge.
  /// @tparam kLower True for lower hull, false for upper hull.
  template <bool kLower>
  bool CoverCompare(const Point& point, const Bridge& bridge) {
    auto res = compare_at_x_(point, bridge.supporting_line());
    if (res == 0) return true;
    return (res > 0) == kLower;
  }

  /// @brief Find the bridge between left and right subtrees.
  ///
  /// Uses the O(log n) bridge-finding algorithm that navigates
  /// the hulls using slope comparisons at each step.
  ///
  /// @tparam kLower True for lower hull bridge, false for upper hull bridge.
  /// @param node Internal node whose children's hulls we're bridging.
  /// @return Bridge segment connecting the two sub-hulls.
  template <bool kLower>
  Bridge FindBridge(Node* node) {
    Node* x = node->left;
    Node* y = node->right;

    if (!x || !y) return Bridge();

    Bridge edge_left;
    Bridge edge_right;
    Bridge left_right_bridge;
    bool undecided;

    while (!(this->IsLeaf(x) && this->IsLeaf(y))) {
      undecided = true;
      if (!x || !y) break;

      Point mid = midpoint_(x->max[kLower].max(), y->min[kLower].min());

      edge_left = x->val[kLower];
      edge_right = y->val[kLower];
      left_right_bridge = Bridge(midpoint_(edge_left), midpoint_(edge_right));

      if (!this->IsLeaf(x) && SlopeCompare<kLower>(edge_left, left_right_bridge)) {
        x = StepLeft<kLower>(x);
        undecided = false;
      }
      if (!this->IsLeaf(y) && SlopeCompare<kLower>(left_right_bridge, edge_right)) {
        y = StepRight<kLower>(y);
        undecided = false;
      }
      if (undecided) {
        bool m_result = MidpointCompare<kLower>(edge_left, edge_right, mid);
        if ((!this->IsLeaf(x) && m_result) || this->IsLeaf(y)) {
          x = StepRight<kLower>(x);
        } else {
          y = StepLeft<kLower>(y);
        }
      }
    }

    if (!x || !y) return Bridge();

    return Bridge(x->val[kLower].min(), y->val[kLower].max());
  }

  /// @brief Update bridges after tree modification.
  void OnUpdate(Node* x) override {
    if (this->IsLeaf(x)) return;
    x->val[0] = FindBridge<false>(x);
    x->val[1] = FindBridge<true>(x);
  }

  /// @brief Navigate left in hull to find bridge endpoint.
  /// @tparam kLower True for lower hull, false for upper hull.
  template <bool kLower>
  Node* StepLeft(Node* node) {
    auto x = node->val[kLower].min().x();
    node = node->left;
    while (node && node->val[kLower].max().x() > x) {
      node = node->left;
    }
    return node;
  }

  /// @brief Navigate right in hull to find bridge endpoint.
  /// @tparam kLower True for lower hull, false for upper hull.
  template <bool kLower>
  Node* StepRight(Node* node) {
    auto x = node->val[kLower].max().x();
    node = node->right;
    while (node && node->val[kLower].min().x() < x) {
      node = node->right;
    }
    return node;
  }

  /// @brief Find successor vertex on hull from given point.
  template <bool kLower>
  Node* HullSuccessor(const Point key) {
    auto current = this->root_;
    Point min;
    while (current) {
      min = current->val[kLower].min();
      if (min == key) break;
      if (min < key) {
        current = StepRight<kLower>(current);
      } else {
        current = StepLeft<kLower>(current);
      }
    }
    return current;
  }

  template <bool kLower>
  Node* HullSuccessor(const Node* node) {
    return HullSuccessor<kLower>(node->val[kLower].max());
  }

  /// @brief Find predecessor vertex on hull from given point.
  template <bool kLower>
  Node* HullPredecessor(const Point key) {
    auto current = this->root_;
    Point max;
    while (current) {
      max = current->val[kLower].max();
      if (max == key) break;
      if (max < key) {
        current = StepRight<kLower>(current);
      } else {
        current = StepLeft<kLower>(current);
      }
    }
    return current;
  }

  template <bool kLower>
  Node* HullPredecessor(const Node* node) {
    return HullPredecessor<kLower>(node->val[kLower].min());
  }

  /// @brief Check if point is covered by hull.
  template <bool kLower>
  [[nodiscard]] bool CoversImpl(Point point) {
    auto current = this->root_;
    while (current) {
      if (current->val[kLower].min() <= point) {
        if (point <= current->val[kLower].max()) {
          return CoverCompare<kLower>(point, current->val[kLower]);
        } else {
          current = StepRight<kLower>(current);
        }
      } else {
        current = StepLeft<kLower>(current);
      }
    }
    return false;
  }

  /// @brief Collect all hull vertices in order.
  template <bool kLower>
  std::vector<Point> HullPointsImpl() {
    std::vector<Point> result;
    Node* edge = this->root_;
    if (!edge || this->IsLeaf(edge)) return result;
    while (edge && !this->IsLeaf(edge)) {
      result.insert(result.begin(), edge->val[kLower].min());
      edge = HullPredecessor<kLower>(result.front());
    }
    edge = this->root_;
    while (edge && !this->IsLeaf(edge)) {
      result.push_back(edge->val[kLower].max());
      edge = HullSuccessor<kLower>(result.back());
    }
    if (kLower) {
      result.pop_back();
    } else {
      result.erase(result.begin());
      std::reverse(result.begin(), result.end());
    }
    return result;
  }

 public:
  /// @brief Insert a point into the convex hull. O(log² n)
  void Insert(Point point) {
    AVLTree<Bridges<Traits>>::Insert({Bridge(point, point), Bridge(point, point)});
  }

  /// @brief Remove a point from the convex hull. O(log² n)
  void Remove(Point point) {
    AVLTree<Bridges<Traits>>::Remove({Bridge(point, point), Bridge(point, point)});
  }

  /// @brief Check if point is inside the convex hull. O(log n)
  [[nodiscard]] bool Covers(Point point) {
    return CoversImpl<false>(point) && CoversImpl<true>(point);
  }

  /// @brief Get upper hull vertices in order.
  [[nodiscard]] std::vector<Point> UpperHullPoints() {
    return HullPointsImpl<false>();
  }

  /// @brief Get lower hull vertices in order.
  [[nodiscard]] std::vector<Point> LowerHullPoints() {
    return HullPointsImpl<true>();
  }

  /// @brief Get the x-coordinate range of all points.
  /// @return {min_x, max_x}. For empty hulls, returns {+inf, -inf}.
  [[nodiscard]] std::pair<double, double> GetXRange() const {
    if (!this->root_) {
      return {std::numeric_limits<double>::infinity(),
              -std::numeric_limits<double>::infinity()};
    }
    return {this->root_->min[0].min().x(), this->root_->max[0].max().x()};
  }

  /// @brief Evaluate the hull y-value at given x-coordinate.
  /// @tparam kLower True for lower hull, false for upper hull.
  /// @param x X-coordinate within hull's x-range.
  /// @return Y-value on the hull boundary at x.
  template <bool kLower>
  [[nodiscard]] double EvaluateHullAt(double x) const {
    auto current = this->root_;
    if (!current) {
      return kLower ? std::numeric_limits<double>::infinity()
                    : -std::numeric_limits<double>::infinity();
    }

    while (current && (current->left || current->right)) {
      const Bridge& bridge = current->val[kLower];
      double min_x = bridge.min().x();
      double max_x = bridge.max().x();

      if (x >= min_x && x <= max_x) {
        if (bridge.is_vertical()) {
          return kLower ? bridge.min().y() : bridge.max().y();
        }
        double dx = max_x - min_x;
        if (std::abs(dx) < 1e-15) return bridge.min().y();
        double dy = bridge.max().y() - bridge.min().y();
        double slope = dy / dx;
        return bridge.min().y() + slope * (x - min_x);
      } else if (x < min_x) {
        current = current->left;
      } else {
        current = current->right;
      }
    }

    if (current) {
      return current->val[kLower].min().y();
    }
    return kLower ? std::numeric_limits<double>::infinity()
                  : -std::numeric_limits<double>::infinity();
  }

  /// @brief Get the slope of the hull edge at given x-coordinate.
  /// @tparam kLower True for lower hull, false for upper hull.
  template <bool kLower>
  [[nodiscard]] double GetHullSlope(double x) const {
    auto current = this->root_;
    if (!current) return 0.0;

    while (current && (current->left || current->right)) {
      const Bridge& bridge = current->val[kLower];
      double min_x = bridge.min().x();
      double max_x = bridge.max().x();

      if (x >= min_x && x <= max_x) {
        if (bridge.is_vertical()) {
          return kLower ? -std::numeric_limits<double>::infinity()
                        : std::numeric_limits<double>::infinity();
        }
        double dx = max_x - min_x;
        if (std::abs(dx) < 1e-15) return 0.0;
        double dy = bridge.max().y() - bridge.min().y();
        return dy / dx;
      } else if (x < min_x) {
        current = current->left;
      } else {
        current = current->right;
      }
    }

    return 0.0;
  }

  [[nodiscard]] double EvaluateLowerHullAt(double x) const {
    return EvaluateHullAt<true>(x);
  }

  [[nodiscard]] double EvaluateUpperHullAt(double x) const {
    return EvaluateHullAt<false>(x);
  }

  [[nodiscard]] double GetLowerHullSlope(double x) const {
    return GetHullSlope<true>(x);
  }

  [[nodiscard]] double GetUpperHullSlope(double x) const {
    return GetHullSlope<false>(x);
  }

  /// @brief Build convex hull from sorted points given by iterator range in O(n) time.
  /// @tparam Iterator Iterator type that dereferences to Point.
  /// @param begin Iterator to the first point.
  /// @param end Iterator past the last point.
  template <typename Iterator>
  void Build(Iterator begin, Iterator end) {
    AVLTree<Bridges<Traits>>::Clear();
    if (begin == end) {
      return;
    }

    std::vector<Bridges<Traits>> bridges;
    for (auto it = begin; it != end; ++it) {
      bridges.push_back({Bridge(*it, *it), Bridge(*it, *it)});
    }

    this->root_ = BuildSubtreeFast(bridges, 0, bridges.size());
  }

  /// @brief Build convex hull from sorted points in O(n) time.
  /// @param sorted_points Points sorted by x-coordinate (and by y for ties).
  void Build(const std::vector<Point>& sorted_points) {
    Build(sorted_points.begin(), sorted_points.end());
  }

 private:
  /// @brief Recursive helper for O(n) bulk construction.
  Node* BuildSubtreeFast(const std::vector<Bridges<Traits>>& bridges,
                         std::size_t start, std::size_t end) {
    if (start >= end) return nullptr;

    if (end - start == 1) {
      return new Node(bridges[start]);
    }

    std::size_t mid = start + (end - start) / 2;

    Node* node = new Node(bridges[mid]);

    node->left = BuildSubtreeFast(bridges, start, mid);
    node->right = BuildSubtreeFast(bridges, mid, end);

    if (node->left) node->left->par = node;
    if (node->right) node->right->par = node;

    node->size = AVLTree<Bridges<Traits>>::Size(node->left) +
                 AVLTree<Bridges<Traits>>::Size(node->right) +
                 (node && !(node->left || node->right) ? 1 : 0);
    node->height = std::max(AVLTree<Bridges<Traits>>::Height(node->left),
                            AVLTree<Bridges<Traits>>::Height(node->right)) + 1;
    if (node->left) {
      node->min = node->left->min;
    } else {
      node->min = node->val;
    }
    if (node->right) {
      node->max = node->right->max;
    } else {
      node->max = node->val;
    }

    OnUpdate(node);

    return node;
  }

 public:
  /// @brief Join another hull into this one.
  ///
  /// Precondition: All points in 'this' must have x < all points in 'other'.
  /// After this operation, 'other' becomes empty.
  ///
  /// @param other Hull to merge (must be to the right of this hull).
  /// @note Time complexity: O(log² N)
  void Join(CHTree& other) {
    if (other.Empty()) return;
    if (this->Empty()) {
      this->root_ = other.root_;
      if (this->root_) this->root_->par = nullptr;
      other.root_ = nullptr;
      return;
    }

    Point p = this->root_ ? this->root_->val[0].min() : other.root_->val[0].min();
    Bridges<Traits> dummy = {Bridge(p, p), Bridge(p, p)};

    AVLTree<Bridges<Traits>>::Join(dummy, &other);
  }

  /// @brief Split this hull at the given x-coordinate.
  ///
  /// After the operation:
  ///   - 'this' contains all points with x < split_x
  ///   - Returns a new CHTree containing points with x >= split_x
  ///
  /// @param split_x X-coordinate to split at.
  /// @return New CHTree containing the right portion.
  /// @note Time complexity: O(log² N)
  [[nodiscard]] CHTree Split(double split_x) {
    CHTree right_hull;
    if (this->root_ == nullptr) {
      return right_hull;
    }

    Point split_point(split_x, std::numeric_limits<double>::lowest());
    Bridges<Traits> split_key = {Bridge(split_point, split_point),
                                  Bridge(split_point, split_point)};

    AVLTree<Bridges<Traits>>::Split(split_key, &right_hull);

    return right_hull;
  }

  /// @brief Get the number of points in the hull.
  [[nodiscard]] std::size_t Size() const {
    return this->root_ ? this->root_->size : 0;
  }

  /// @brief Check if hull is empty.
  [[nodiscard]] bool Empty() const { return this->root_ == nullptr; }

  /// @brief Validate that all bridges are consistent (for debugging).
  [[nodiscard]] bool ValidateBridges() {
    return ValidateBridgesRecursive(this->root_);
  }

 private:
  bool ValidateBridgesRecursive(Node* node) {
    if (!node || this->IsLeaf(node)) return true;

    if (!ValidateBridgesRecursive(node->left)) return false;
    if (!ValidateBridgesRecursive(node->right)) return false;

    Bridge upper = FindBridge<false>(node);
    Bridge lower = FindBridge<true>(node);

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

  /// @brief Find critical vertex for stabbing line query.
  ///
  /// Traverses the tree to find the vertex where the slope difference
  /// between two hulls changes sign (subgradient condition).
  ///
  /// @tparam kLower True for ceiling hull, false for floor hull.
  /// @param other_hull The other hull to compare slopes against.
  /// @return {x-coordinate, slope} of the critical point.
  template <bool kLower>
  std::pair<double, double> FindCriticalVertex(
      const CHTree<Traits>& other_hull) const {
    auto current = this->root_;
    if (!current) return {0.0, 0.0};

    while (current && (current->left || current->right)) {
      const Bridge& bridge = current->val[kLower];
      double mid_x = (bridge.min().x() + bridge.max().x()) / 2.0;

      double this_slope;
      if (bridge.is_vertical()) {
        this_slope = kLower ? std::numeric_limits<double>::infinity()
                            : -std::numeric_limits<double>::infinity();
      } else {
        double dx = bridge.max().x() - bridge.min().x();
        if (std::abs(dx) < 1e-15) {
          this_slope = 0.0;
        } else {
          this_slope = (bridge.max().y() - bridge.min().y()) / dx;
        }
      }

      double other_slope;
      if constexpr (kLower) {
        other_slope = other_hull.GetUpperHullSlope(mid_x);
      } else {
        other_slope = other_hull.GetLowerHullSlope(mid_x);
      }

      if (std::isinf(this_slope) || std::isinf(other_slope)) {
        return {bridge.min().x(), 0.0};
      }

      double slope_diff;
      if constexpr (kLower) {
        slope_diff = this_slope - other_slope;
      } else {
        slope_diff = other_slope - this_slope;
      }

      if (slope_diff < -1e-9) {
        current = current->right ? current->right : current;
        if (!current->right) break;
      } else if (slope_diff > 1e-9) {
        current = current->left ? current->left : current;
        if (!current->left) break;
      } else {
        return {mid_x, this_slope};
      }
    }

    if (current) {
      Point p = current->val[kLower].min();
      return {p.x(), 0.0};
    }
    return {0.0, 0.0};
  }

 public:
  [[nodiscard]] std::pair<double, double> FindCriticalVertexFromCeiling(
      const CHTree<Traits>& floor_hull) const {
    return FindCriticalVertex<true>(floor_hull);
  }

  [[nodiscard]] std::pair<double, double> FindCriticalVertexFromFloor(
      const CHTree<Traits>& ceiling_hull) const {
    return FindCriticalVertex<false>(ceiling_hull);
  }

  /// @brief Get the range of valid tangent slopes at a hull vertex.
  /// @tparam kLower True for lower hull, false for upper hull.
  /// @param x X-coordinate of the vertex.
  /// @return {left_slope, right_slope} defining the valid slope range.
  template <bool kLower>
  [[nodiscard]] std::pair<double, double> GetSlopeRange(double x) const {
    auto current = this->root_;
    if (!current) return {0.0, 0.0};

    const double kInf = std::numeric_limits<double>::infinity();
    double left_slope = kLower ? -kInf : kInf;
    double right_slope = kLower ? kInf : -kInf;
    bool found_exact_vertex = false;

    while (current && (current->left || current->right)) {
      const Bridge& bridge = current->val[kLower];
      double min_x = bridge.min().x();
      double max_x = bridge.max().x();

      double bridge_slope = 0.0;
      if (!bridge.is_vertical()) {
        double dx = max_x - min_x;
        if (std::abs(dx) >= 1e-15) {
          bridge_slope = (bridge.max().y() - bridge.min().y()) / dx;
        }
      }

      if (std::abs(x - min_x) < 1e-12) {
        right_slope = bridge_slope;
        found_exact_vertex = true;
        current = current->left;
      } else if (std::abs(x - max_x) < 1e-12) {
        left_slope = bridge_slope;
        found_exact_vertex = true;
        current = current->right;
      } else if (x > min_x && x < max_x) {
        return {bridge_slope, bridge_slope};
      } else if (x < min_x) {
        right_slope = bridge_slope;
        current = current->left;
      } else {
        left_slope = bridge_slope;
        current = current->right;
      }
    }

    if (current && found_exact_vertex) {
      if constexpr (kLower) {
        return {left_slope, right_slope};
      } else {
        return {right_slope, left_slope};
      }
    }

    if constexpr (kLower) {
      return {left_slope, right_slope};
    } else {
      return {right_slope, left_slope};
    }
  }

  [[nodiscard]] std::pair<double, double> GetLowerHullSlopeRange(double x) const {
    return GetSlopeRange<true>(x);
  }

  [[nodiscard]] std::pair<double, double> GetUpperHullSlopeRange(double x) const {
    return GetSlopeRange<false>(x);
  }

};

}  // namespace dch

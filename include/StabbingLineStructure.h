// Copyright 2026 DynamicConvexHull Authors
// SPDX-License-Identifier: MIT

/// @file StabbingLineStructure.h
/// @brief Query for epsilon-stabbing lines in O(log² n).
///
/// A "stabbing line" is a line that passes within vertical distance epsilon
/// of every point in the set. This structure maintains a single convex hull
/// and applies epsilon adjustments during queries.
///
/// Time Complexities:
///   - insert/remove: O(log² n)
///   - HasStabbingLine/FindStabbingLine: O(log² n)
///   - split/join: O(log² n)
///   - build: O(n)

#pragma once

#include <cmath>
#include <limits>
#include <optional>
#include <vector>

#include "CHTree.h"

namespace dch {

/// @brief Represents a stabbing line y = slope * x + intercept.
template <class Traits>
struct StabbingLine {
  using Point = typename Traits::Point_2;
  double slope;
  double intercept;

  StabbingLine() : slope(0), intercept(0) {}
  StabbingLine(double s, double i) : slope(s), intercept(i) {}

  [[nodiscard]] double At(double x) const { return slope * x + intercept; }

  [[nodiscard]] bool CoversPoint(const Point& point, double epsilon) const {
    double line_y = At(point.x());
    return std::abs(point.y() - line_y) <= epsilon + 1e-9;
  }
};

/// @brief Data structure for querying epsilon-stabbing lines.
///
/// Uses a single convex hull and applies epsilon adjustments at query time:
/// - For "ceiling" queries (upper boundary): lower hull y + epsilon
/// - For "floor" queries (lower boundary): upper hull y - epsilon
///
/// @tparam Traits CGAL-style kernel.
template <class Traits>
class StabbingLineStructure {
 public:
  using Point = typename Traits::Point_2;
  using Line = StabbingLine<Traits>;

  explicit StabbingLineStructure(double eps = 1.0) : epsilon_(eps) {}

  [[nodiscard]] double GetEpsilon() const { return epsilon_; }
  [[nodiscard]] std::size_t Size() const { return point_count_; }
  [[nodiscard]] bool Empty() const { return point_count_ == 0; }

  /// @brief Insert a point. O(log² n)
  void Insert(const Point& point) {
    hull_.Insert(point);
    ++point_count_;
  }

  /// @brief Remove a point. O(log² n)
  void Remove(const Point& point) {
    hull_.Remove(point);
    if (point_count_ > 0) --point_count_;
  }

  /// @brief Build from sorted points given by iterator range. O(n)
  /// @tparam Iterator Iterator type that dereferences to Point.
  /// @param begin Iterator to the first point.
  /// @param end Iterator past the last point.
  template <typename Iterator>
  void Build(Iterator begin, Iterator end) {
    hull_.Build(begin, end);
    point_count_ = hull_.Size();
  }

  /// @brief Build from sorted points. O(n)
  void Build(const std::vector<Point>& sorted_points) {
    Build(sorted_points.begin(), sorted_points.end());
  }

  /// @brief Split at x-coordinate. O(log² n)
  [[nodiscard]] StabbingLineStructure Split(double split_x) {
    StabbingLineStructure right(epsilon_);
    auto hull_right = hull_.Split(split_x);
    right.hull_ = std::move(hull_right);
    right.point_count_ = right.hull_.Size();
    point_count_ = hull_.Size();
    return right;
  }

  /// @brief Join with another structure. O(log² n)
  void Join(StabbingLineStructure& other) {
    hull_.Join(other.hull_);
    point_count_ = hull_.Size();
    other.point_count_ = 0;
  }

  /// @brief Find a stabbing line if one exists.
  /// @return Line if exists, nullopt otherwise.
  /// @note Time complexity: O(log² n)
  [[nodiscard]] std::optional<Line> FindStabbingLine() const {
    if (point_count_ == 0) {
      return Line(0.0, 0.0);
    }

    auto [min_x, max_x] = hull_.GetXRange();

    if (min_x > max_x) {
      return Line(0.0, 0.0);
    }

    if (std::abs(max_x - min_x) < 1e-12) {
      // Single x-coordinate: ceiling = lower hull + epsilon,
      //                      floor = upper hull - epsilon
      double y_ceil = hull_.EvaluateLowerHullAt(min_x) + epsilon_;
      double y_floor = hull_.EvaluateUpperHullAt(min_x) - epsilon_;
      double gap = y_ceil - y_floor;
      if (gap >= -1e-9) {
        return Line(0.0, (y_ceil + y_floor) / 2.0);
      }
      return std::nullopt;
    }

    // Find critical vertices - comparing hull against itself since both
    // ceiling and floor are derived from the same hull.
    auto [x_ceil, slope_ceil] = hull_.FindCriticalVertexFromCeiling(hull_);
    auto [x_floor, slope_floor] = hull_.FindCriticalVertexFromFloor(hull_);

    double gap_ceil = ComputeGap(x_ceil);
    double gap_floor = ComputeGap(x_floor);
    double gap_min = ComputeGap(min_x);
    double gap_max = ComputeGap(max_x);

    double best_gap = gap_ceil;
    double best_x = x_ceil;

    if (gap_floor < best_gap) {
      best_gap = gap_floor;
      best_x = x_floor;
    }
    if (gap_min < best_gap) {
      best_gap = gap_min;
      best_x = min_x;
    }
    if (gap_max < best_gap) {
      best_gap = gap_max;
      best_x = max_x;
    }

    if (best_gap < -1e-9) {
      return std::nullopt;
    }

    // Ceiling = lower hull + epsilon, Floor = upper hull - epsilon
    double y_ceil = hull_.EvaluateLowerHullAt(best_x) + epsilon_;
    double y_floor = hull_.EvaluateUpperHullAt(best_x) - epsilon_;
    double y_mid = (y_ceil + y_floor) / 2.0;

    // Get slope ranges - for ceiling use lower hull, for floor use upper hull
    auto [ceil_left, ceil_right] = hull_.GetLowerHullSlopeRange(best_x);
    auto [floor_left, floor_right] = hull_.GetUpperHullSlopeRange(best_x);

    if (std::isinf(ceil_left) && ceil_left < 0) ceil_left = -1e15;
    if (std::isinf(ceil_right) && ceil_right > 0) ceil_right = 1e15;
    if (std::isinf(floor_left) && floor_left < 0) floor_left = -1e15;
    if (std::isinf(floor_right) && floor_right > 0) floor_right = 1e15;

    double min_valid = std::max(ceil_left, floor_left);
    double max_valid = std::min(ceil_right, floor_right);

    double slope;
    if (min_valid <= max_valid + 1e-9) {
      slope = (min_valid + max_valid) / 2.0;
    } else {
      return std::nullopt;
    }

    double gap = y_ceil - y_floor;
    double intercept = (y_ceil - gap / 2.0) - slope * best_x;

    auto verify_gap_with_slope = [&](double x) -> bool {
      double line_y = slope * x + intercept;
      double y_ceil_at = hull_.EvaluateLowerHullAt(x) + epsilon_;
      double y_floor_at = hull_.EvaluateUpperHullAt(x) - epsilon_;
      return (line_y <= y_ceil_at + 1e-9) && (line_y >= y_floor_at - 1e-9);
    };

    if (!verify_gap_with_slope(min_x) || !verify_gap_with_slope(max_x) ||
        !verify_gap_with_slope(best_x)) {
      return std::nullopt;
    }

    return Line(slope, intercept);
  }

  /// @brief Check if a stabbing line exists. O(log² n)
  [[nodiscard]] bool HasStabbingLine() const {
    return FindStabbingLine().has_value();
  }

  /// @brief Get the minimum gap (ceiling - floor).
  [[nodiscard]] double GetMinimumGap() const {
    if (point_count_ == 0) return std::numeric_limits<double>::infinity();

    auto [min_x, max_x] = hull_.GetXRange();
    if (min_x > max_x) return std::numeric_limits<double>::infinity();

    auto [x_ceil, unused1] = hull_.FindCriticalVertexFromCeiling(hull_);
    auto [x_floor, unused2] = hull_.FindCriticalVertexFromFloor(hull_);

    double gap = ComputeGap(x_ceil);
    gap = std::min(gap, ComputeGap(x_floor));
    gap = std::min(gap, ComputeGap(min_x));
    gap = std::min(gap, ComputeGap(max_x));
    return gap;
  }

  [[nodiscard]] std::pair<double, double> GetXRange() const {
    return hull_.GetXRange();
  }

 private:
  CHTree<Traits> hull_;  ///< Single hull storing original points
  double epsilon_;
  std::size_t point_count_ = 0;

  /// @brief Compute the gap between ceiling and floor at x.
  ///
  /// Ceiling = lower hull + epsilon (upper boundary of epsilon band)
  /// Floor = upper hull - epsilon (lower boundary of epsilon band)
  /// Gap = ceiling - floor = (lower + ε) - (upper - ε) = lower - upper + 2ε
  [[nodiscard]] double ComputeGap(double x) const {
    double y_lower = hull_.EvaluateLowerHullAt(x);
    double y_upper = hull_.EvaluateUpperHullAt(x);
    return (y_lower + epsilon_) - (y_upper - epsilon_);
  }
};

}  // namespace dch

// Copyright 2026 DynamicConvexHull Authors
// SPDX-License-Identifier: MIT

/// @file inexact.h
/// @brief Inexact (floating-point) kernel for convex hull operations.
///
/// Provides a lightweight alternative to CGAL exact kernels for cases where
/// floating-point precision is acceptable. Uses C++20 std::partial_ordering.

#pragma once

#include <compare>
#include <cstddef>

namespace dch {

/// @brief Inexact geometric kernel using floating-point arithmetic.
/// @tparam T Numeric type (typically double).
template <typename T>
class InexactKernel {
 public:
  struct Point_2 {
    T x_;
    T y_;

    Point_2() = default;
    Point_2(T x, T y) : x_(x), y_(y) {}

    [[nodiscard]] T x() const { return x_; }
    [[nodiscard]] T y() const { return y_; }

    bool operator==(const Point_2& other) const {
      return x_ == other.x_ && y_ == other.y_;
    }

    bool operator<(const Point_2& other) const {
      if (x_ < other.x_) return true;
      return (x_ == other.x_ && y_ < other.y_);
    }

    bool operator<=(const Point_2& other) const {
      if (x_ < other.x_) return true;
      return (x_ == other.x_ && y_ <= other.y_);
    }
  };

  struct Line_2 {
    T slope;
    T intercept;

    [[nodiscard]] T Eval(T x) const { return intercept + x * slope; }
  };

  struct Segment_2 {
    Point_2 source;
    Point_2 target;

    [[nodiscard]] const Point_2& min() const { return source; }
    [[nodiscard]] const Point_2& max() const { return target; }

    Segment_2() = default;
    Segment_2(Point_2 a, Point_2 b) : source(a), target(b) {}
    Segment_2(const Segment_2& s) = default;

    [[nodiscard]] bool is_vertical() const { return source.x_ == target.x_; }

    [[nodiscard]] T slope() const {
      return (target.y_ - source.y_) / (target.x_ - source.x_);
    }

    [[nodiscard]] Line_2 supporting_line() const {
      T s = slope();
      return {s, source.y_ - s * source.x_};
    }

    bool operator==(const Segment_2& other) const {
      return source == other.source && target == other.target;
    }

    bool operator!=(const Segment_2& other) const { return !(*this == other); }

    Point_2& operator[](std::size_t idx) { return idx % 2 ? target : source; }

    const Point_2& operator[](std::size_t idx) const {
      return idx % 2 ? target : source;
    }
  };

  class Construct_midpoint_2 {
   public:
    [[nodiscard]] Point_2 operator()(const Point_2& a, const Point_2& b) const {
      return {(b.x_ - a.x_) / 2.0 + a.x_, (b.y_ - a.y_) / 2.0 + a.y_};
    }

    [[nodiscard]] Point_2 operator()(const Segment_2& s) const {
      return this->operator()(s.source, s.target);
    }
  };

  class Compare_slope_2 {
   public:
    [[nodiscard]] std::partial_ordering operator()(const Segment_2& left,
                                                   const Segment_2& right) const {
      T left_slope = left.slope();
      T right_slope = right.slope();
      return left_slope <=> right_slope;
    }
  };

  class Compare_y_at_x_2 {
   public:
    [[nodiscard]] std::partial_ordering operator()(const Point_2& point,
                                                   const Line_2& left,
                                                   const Line_2& right) const {
      T left_val = left.Eval(point.x_);
      T right_val = right.Eval(point.x_);
      return left_val <=> right_val;
    }

    [[nodiscard]] std::partial_ordering operator()(const Point_2& point,
                                                   const Line_2& line) const {
      T val = line.Eval(point.x_);
      return val <=> point.y_;
    }
  };

  class Compare_xy_2 {
   public:
    [[nodiscard]] std::partial_ordering operator()(const Point_2& left,
                                                   const Point_2& right) const {
      if (auto cmp = left.x_ <=> right.x_; cmp != 0) return cmp;
      return left.y_ <=> right.y_;
    }
  };
};

// Type alias for backward compatibility and convenience.
template <typename T>
using Inexact_kernel = InexactKernel<T>;

/// @brief Integer coordinate kernel for convex hull operations.
///
/// Uses integer types for x/y coordinates to avoid precision loss with large
/// values (e.g., uint64_t keys > 2^53). Slope and intercept remain as floats.
///
/// @tparam IntType Integer type for coordinates (e.g., int64_t, uint64_t).
/// @tparam FloatType Floating-point type for slope/intercept (typically double).
template <typename IntType, typename FloatType = double>
class IntegerKernel {
 public:
  struct Point_2 {
    IntType x_;
    IntType y_;

    Point_2() = default;
    Point_2(IntType x, IntType y) : x_(x), y_(y) {}
    
    // Allow construction from compatible types
    template <typename T1, typename T2>
    Point_2(T1 x, T2 y) : x_(static_cast<IntType>(x)), y_(static_cast<IntType>(y)) {}

    [[nodiscard]] IntType x() const { return x_; }
    [[nodiscard]] IntType y() const { return y_; }
    
    // Provide double versions for compatibility with existing code
    [[nodiscard]] FloatType x_float() const { return static_cast<FloatType>(x_); }
    [[nodiscard]] FloatType y_float() const { return static_cast<FloatType>(y_); }

    bool operator==(const Point_2& other) const {
      return x_ == other.x_ && y_ == other.y_;
    }

    bool operator<(const Point_2& other) const {
      if (x_ < other.x_) return true;
      return (x_ == other.x_ && y_ < other.y_);
    }

    bool operator<=(const Point_2& other) const {
      if (x_ < other.x_) return true;
      return (x_ == other.x_ && y_ <= other.y_);
    }
  };

  struct Line_2 {
    FloatType slope;
    FloatType intercept;

    [[nodiscard]] FloatType Eval(FloatType x) const { return intercept + x * slope; }
    [[nodiscard]] FloatType Eval(IntType x) const { 
      return intercept + static_cast<FloatType>(x) * slope; 
    }
  };

  struct Segment_2 {
    Point_2 source;
    Point_2 target;

    [[nodiscard]] const Point_2& min() const { return source; }
    [[nodiscard]] const Point_2& max() const { return target; }

    Segment_2() = default;
    Segment_2(Point_2 a, Point_2 b) : source(a), target(b) {}
    Segment_2(const Segment_2& s) = default;

    [[nodiscard]] bool is_vertical() const { return source.x_ == target.x_; }

    [[nodiscard]] FloatType slope() const {
      return static_cast<FloatType>(target.y_ - source.y_) / 
             static_cast<FloatType>(target.x_ - source.x_);
    }

    [[nodiscard]] Line_2 supporting_line() const {
      FloatType s = slope();
      return {s, static_cast<FloatType>(source.y_) - s * static_cast<FloatType>(source.x_)};
    }

    bool operator==(const Segment_2& other) const {
      return source == other.source && target == other.target;
    }

    bool operator!=(const Segment_2& other) const { return !(*this == other); }

    Point_2& operator[](std::size_t idx) { return idx % 2 ? target : source; }

    const Point_2& operator[](std::size_t idx) const {
      return idx % 2 ? target : source;
    }
  };

  class Construct_midpoint_2 {
   public:
    [[nodiscard]] Point_2 operator()(const Point_2& a, const Point_2& b) const {
      // Use integer arithmetic for midpoint (floor division)
      return {a.x_ + (b.x_ - a.x_) / 2, a.y_ + (b.y_ - a.y_) / 2};
    }

    [[nodiscard]] Point_2 operator()(const Segment_2& s) const {
      return this->operator()(s.source, s.target);
    }
  };

  class Compare_slope_2 {
   public:
    [[nodiscard]] std::partial_ordering operator()(const Segment_2& left,
                                                   const Segment_2& right) const {
      FloatType left_slope = left.slope();
      FloatType right_slope = right.slope();
      return left_slope <=> right_slope;
    }
  };

  class Compare_y_at_x_2 {
   public:
    [[nodiscard]] std::partial_ordering operator()(const Point_2& point,
                                                   const Line_2& left,
                                                   const Line_2& right) const {
      FloatType left_val = left.Eval(point.x_);
      FloatType right_val = right.Eval(point.x_);
      return left_val <=> right_val;
    }

    [[nodiscard]] std::partial_ordering operator()(const Point_2& point,
                                                   const Line_2& line) const {
      FloatType val = line.Eval(point.x_);
      return val <=> static_cast<FloatType>(point.y_);
    }
  };

  class Compare_xy_2 {
   public:
    [[nodiscard]] std::partial_ordering operator()(const Point_2& left,
                                                   const Point_2& right) const {
      if (auto cmp = left.x_ <=> right.x_; cmp != 0) return cmp;
      return left.y_ <=> right.y_;
    }
  };
};

}  // namespace dch

// Copyright 2026 DynamicConvexHull Authors
// SPDX-License-Identifier: MIT

/// @file util.h
/// @brief Utility types for the Dynamic Convex Hull structure.

#pragma once

#include <array>
#include <compare>

namespace dch {

/// @brief Pair of bridges (upper and lower) stored at each internal node.
/// @tparam Traits Kernel providing Point_2, Segment_2, and comparators.
template <class Traits>
struct Bridges {
  using Compare = typename Traits::Compare_xy_2;
  using Bridge = typename Traits::Segment_2;

  static constexpr Compare compare = Compare();

  std::array<Bridge, 2> data;

  Bridges(Bridge upper, Bridge lower) : data({upper, lower}) {}

  Bridge& operator[](std::size_t idx) { return data[idx % 2]; }
  const Bridge& operator[](std::size_t idx) const { return data[idx % 2]; }

  bool operator==(const Bridges& other) const {
    return compare(data[0][0], other[0][0]) == 0;
  }

  bool operator!=(const Bridges& other) const { return !(*this == other); }

  bool operator<(const Bridges& other) const {
    return compare(data[0][0], other[0][0]) < 0;
  }

  bool operator<=(const Bridges& other) const {
    return compare(data[0][0], other[0][0]) <= 0;
  }
};

}  // namespace dch

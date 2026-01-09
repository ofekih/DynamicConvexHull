/**
 * @file util.h
 * @brief Utility types for the Dynamic Convex Hull structure.
 */

#ifndef DYNAMICCONVEXHULL_UTIL_H
#define DYNAMICCONVEXHULL_UTIL_H

#include <CGAL/enum.h>
#include <array>

/**
 * @brief Pair of bridges (upper and lower) stored at each internal node.
 * @tparam Traits CGAL-style kernel.
 */
template<class Traits>
struct Bridges{
    using Compare = typename Traits::Compare_xy_2;
    using Bridge = typename Traits::Segment_2;

    static constexpr Compare compare = Compare();

    std::array<Bridge, 2> data;

    Bridges(Bridge x, Bridge y) : data({x, y}) {}

    Bridge& operator[](size_t idx) { return data[idx % 2]; }
    const Bridge& operator[](size_t idx) const { return data[idx % 2]; }

    bool operator==(const Bridges& b) const {
        return (compare(data[0][0], b[0][0]) == CGAL::EQUAL);
    }

    bool operator!=(const Bridges& b) const {
        return !(*this == b);
    }

    bool operator<(const Bridges& b) const {
        return (compare(data[0][0], b[0][0]) == CGAL::SMALLER);
    }

    bool operator<=(const Bridges& b) const {
        return (compare(data[0][0], b[0][0]) != CGAL::LARGER);
    }
};

#endif //DYNAMICCONVEXHULL_UTIL_H

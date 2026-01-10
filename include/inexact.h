/**
 * @file inexact.h
 * @brief Inexact (floating-point) kernel for convex hull operations.
 * 
 * Provides a lightweight alternative to CGAL exact kernels for cases where
 * floating-point precision is acceptable. Uses C++20 std::partial_ordering.
 */

#ifndef DYNAMICCONVEXHULL_INEXACT_H
#define DYNAMICCONVEXHULL_INEXACT_H

#include <compare>

/**
 * @brief Inexact geometric kernel using floating-point arithmetic.
 * @tparam T Numeric type (typically double).
 */
template<typename T>
class Inexact_kernel {
public:
    struct Point_2 {
        T a, b;

        Point_2() = default;
        Point_2(T x, T y) : a(x), b(y) {}

        T x() const { return a; }
        T y() const { return b; }

        bool operator==(const Point_2& o) const {
            return a == o.a && b == o.b;
        }
        
        bool operator<(const Point_2& o) const {
            if (a < o.a) return true;
            return (a == o.a && b < o.b);
        }
        
        bool operator<=(const Point_2& o) const {
            if (a < o.a) return true;
            return (a == o.a && b <= o.b);
        }
    };

    struct Line_2 {
        T slope, intercept;

        T eval(T x) const {
            return intercept + x * slope;
        }
    };

    struct Segment_2 {
        Point_2 a, b;

        const Point_2& min() const { return a; }
        const Point_2& max() const { return b; }

        Segment_2() = default;
        Segment_2(Point_2 a, Point_2 b) : a(a), b(b) {}
        Segment_2(const Segment_2& s) : a(s.a), b(s.b) {}

        [[nodiscard]] bool is_vertical() const {
            return a.a == b.a;
        }

        T slope() const {
            return (b.b - a.b) / (b.a - a.a);
        }

        Line_2 supporting_line() const {
            T s = slope();
            return {s, a.b - s * a.a};
        }

        bool operator==(const Segment_2& o) const {
            return a == o.a && b == o.b;
        }
        
        bool operator!=(const Segment_2& o) const {
            return !(*this == o);
        }
        
        Point_2& operator[](size_t idx) { return idx % 2 ? b : a; }
        const Point_2& operator[](size_t idx) const { return idx % 2 ? b : a; }
    };

    class Construct_midpoint_2 {
    public:
        Point_2 operator()(const Point_2& a, const Point_2& b) const {
            return {(b.a - a.a) / 2. + a.a, (b.b - a.b) / 2. + a.b};
        }
        
        Point_2 operator()(const Segment_2& s) const {
            return this->operator()(s.a, s.b);
        }
    };

    class Compare_slope_2 {
    public:
        std::partial_ordering operator()(const Segment_2& l, const Segment_2& r) const {
            T lslope = l.slope();
            T rslope = r.slope();
            return lslope <=> rslope;
        }
    };

    class Compare_y_at_x_2 {
    public:
        std::partial_ordering operator()(const Point_2& x, const Line_2& l, const Line_2& r) const {
            T lval = l.eval(x.a);
            T rval = r.eval(x.a);
            return lval <=> rval;
        }
        
        std::partial_ordering operator()(const Point_2& p, const Line_2& l) const {
            T val = l.eval(p.a);
            return val <=> p.b;
        }
    };

    class Compare_xy_2 {
    public:
        std::partial_ordering operator()(const Point_2& l, const Point_2& r) const {
            if (auto cmp = l.a <=> r.a; cmp != 0) return cmp;
            return l.b <=> r.b;
        }
    };
};

#endif //DYNAMICCONVEXHULL_INEXACT_H

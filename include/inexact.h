#ifndef DYNAMICCONVEXHULL_INEXACT_H
#define DYNAMICCONVEXHULL_INEXACT_H

#include <CGAL/enum.h>

template<typename T>
class Inexact_kernel {
public:
    struct Point_2{
        T a,b;

        Point_2()= default;
        Point_2(T x, T y) : a(x),b(y) {};

        T x(){return a;}
        T y(){return b;}

        bool operator==(const Point_2& o) const {
            return a == o.a && b == o.b;
        }
        bool operator<(const Point_2& o) const {
            if(a < o.a) return true;
            else return (a == o.a && b < o.b);
        }
        bool operator<=(const Point_2& o) const{
            if(a < o.a) return true;
            else return (a == o.a && b <= o.b);
        }
    };
    struct Line_2{
        T slope,intercept;

        T eval(T x) const {
            return intercept + x*slope;
        }

    };
    struct Segment_2{
        Point_2 a,b;

        Point_2 min(){return a;}
        Point_2 max(){return b;}

        Segment_2()= default;
        Segment_2(Point_2 a, Point_2 b) : a(a),b(b) {}
        Segment_2(Segment_2& s): a(s.a),b(s.b) {};

        [[nodiscard]] bool is_vertical() const {
            return a.a == b.a;
        }


        T slope() const {
            return (b.b-a.b)/(b.a-a.a);
        }

        Line_2 supporting_line() {
            T s = slope();
            return {s, a.b - s*a.a};
        }
        Line_2 supporting_line() const {
            T s = slope();
            return {s,a.b-s*a.a};
        }

        bool operator==(const Segment_2& o){
            return a == o.a && b == o.b;
        }
        Point_2& operator[](size_t idx){return idx%2 ? b:a;}
        const Point_2& operator[](size_t idx) const {return idx%2 ? b:a;}
    };


    class Construct_midpoint_2{
    public:
        Point_2 operator()(Point_2& a, Point_2& b){
            return {(b.a-a.a)/2.+a.a,(b.b-a.b)/2.+a.b};
        }
        Point_2 operator()(const Point_2& a, const Point_2& b) const {
            return {(b.a-a.a)/2.+a.a,(b.b-a.b)/2.+a.b};
        }
        Point_2 operator()(Segment_2& s){
            return this->operator()(s.a,s.b);
        }
        Point_2 operator()(const Segment_2& s) const {
            return this->operator()(s.a,s.b);
        }
    };

    class Compare_slope_2{
    public:
        CGAL::Comparison_result operator()(const Segment_2& l, const Segment_2& r) const {
            T lslope = l.slope();
            T rslope = r.slope();
            if(lslope < rslope) return CGAL::SMALLER;
            else if(lslope > rslope) return CGAL::LARGER;
            else return CGAL::EQUAL;
        };

    };
    class Compare_y_at_x_2{
    public:
        CGAL::Comparison_result operator()(const Point_2& x, const Line_2& l, const Line_2& r) const {
            T lval = l.eval(x.a);
            T rval = r.eval(x.a);
            if(lval < rval) return CGAL::SMALLER;
            else if(lval > rval) return CGAL::LARGER;
            else return CGAL::EQUAL;
        }
        CGAL::Comparison_result operator()(const Point_2& p, const Line_2& l) const {
            T val = l.eval(p.a);
            if(val < p.b) return CGAL::LARGER;
            if(val > p.b) return CGAL::SMALLER;
            else return CGAL::EQUAL;
        }
    };

    class Compare_xy_2{
    public:
        CGAL::Comparison_result operator()( Point_2& l, Point_2& r) {
            if(l.a < r.a) return CGAL::SMALLER;
            else if(l.a > r.a) return CGAL::LARGER;
            else {
                if(l.b < r.b) return CGAL::SMALLER;
                else if( l.b > r.b) return CGAL::LARGER;
                else {
                    return CGAL::EQUAL;
                }
            }
        }
        CGAL::Comparison_result operator()(const Point_2& l, const Point_2& r) const {
            if(l.a < r.a) return CGAL::SMALLER;
            else if(l.a > r.a) return CGAL::LARGER;
            else {
                if(l.b < r.b) return CGAL::SMALLER;
                else if( l.b > r.b) return CGAL::LARGER;
                else {
                    return CGAL::EQUAL;
                }
            }
        }
    };
};

#endif //DYNAMICCONVEXHULL_INEXACT_H

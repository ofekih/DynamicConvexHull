/**
 * @file main.cpp
 * @brief Demo and testing program for Dynamic Convex Hull
 */

#include <vector>
#include <random>
#include <chrono>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "include/CHTree.h"
#include "inexact.h"

using namespace dch;

using K = Inexact_kernel<double>;
using Point_2 = K::Point_2;

const double PI = std::acos(-1);

// Simple monotone chain reference implementation
namespace reference {

inline long long cross(const std::pair<double,double>& O, 
                       const std::pair<double,double>& A, 
                       const std::pair<double,double>& B) {
    double dx1 = A.first - O.first;
    double dy1 = A.second - O.second;
    double dx2 = B.first - O.first;
    double dy2 = B.second - O.second;
    return dx1 * dy2 - dy1 * dx2;
}

std::vector<std::pair<double,double>> upperHull(std::vector<std::pair<double,double>> points) {
    std::sort(points.begin(), points.end());
    std::vector<std::pair<double,double>> upper;
    for (const auto& p : points) {
        while (upper.size() >= 2 && cross(upper[upper.size()-2], upper[upper.size()-1], p) >= 0) {
            upper.pop_back();
        }
        upper.push_back(p);
    }
    return upper;
}

std::vector<std::pair<double,double>> lowerHull(std::vector<std::pair<double,double>> points) {
    std::sort(points.begin(), points.end());
    std::vector<std::pair<double,double>> lower;
    for (const auto& p : points) {
        while (lower.size() >= 2 && cross(lower[lower.size()-2], lower[lower.size()-1], p) <= 0) {
            lower.pop_back();
        }
        lower.push_back(p);
    }
    return lower;
}

} // namespace reference

/**
 * @brief Generate test data with various distributions
 */
std::vector<std::pair<double,double>> generate_data(const size_t n, const int mode, const int seed) {
    std::vector<std::pair<double,double>> data(n);
    std::mt19937 engine(seed);

    if(mode == 0){
        double r = 2147483647;
        std::normal_distribution<double> d {0.,r/(2+std::log(n))};
        std::generate(data.begin(),data.end(),[&](){return std::make_pair(d(engine),d(engine));});
    } else if(mode == 1){
        double r = 1000;
        std::uniform_real_distribution<double> d {-r,r};
        std::generate(data.begin(),data.end(),[&](){return std::make_pair(d(engine),d(engine));});
    } else if(mode == 2){
        double r = 1000;
        std::uniform_real_distribution<double> d {-r,r};
        std::generate(data.begin(),data.end(),[&](){
            auto x = d(engine);
            auto y = d(engine);
            while(x*x + y*y > r*r){
                x = d(engine);
                y = d(engine);
            }
            return std::make_pair(x,y);
        });
    } else if(mode == 3){
        double r = 1000;
        std::uniform_real_distribution<double> d {0.,2*PI};
        std::generate(data.begin(),data.end(),[&](){
            auto a = d(engine);
            return std::make_pair(r*std::cos(a),r*std::sin(a));
        });
    }
    return data;
}

/**
 * @brief Runtime test measuring insert/remove performance
 */
void runtimeTest(int size, int mode, int seed) {
    std::cout << "MODE = " << mode << std::endl;
    std::vector<std::pair<double,double>> data = generate_data(size,mode,seed);
    
    CHTree<K> CH;
    std::vector<Point_2> out;
    auto start = std::chrono::steady_clock::now();
    for(const auto& p : data) CH.Insert(Point_2(p.first,p.second));
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "CHTree insert time: " << elapsed.count() << "ms\n";

    start = std::chrono::steady_clock::now();
    for(const auto& p : data) CH.Remove(Point_2(p.first,p.second));
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "CHTree remove time: " << elapsed.count() << "ms\n";
}

/**
 * @brief Verify CHTree against reference implementation
 */
bool verify(CHTree<K>& CH, std::vector<std::pair<double,double>>& data, int size){
    std::vector<std::pair<double,double>> temp(data.begin(), data.begin() + size);
    
    auto refUpper = reference::upperHull(temp);
    auto res = CH.UpperHullPoints();
    
    // Skip leftmost point in CHTree upper hull for comparison
    if (res.size() > 1 && refUpper.size() > 1) {
        for(size_t i = 1; i < res.size() && i < refUpper.size(); i++){
            if(std::abs(res[i].x() - refUpper[i].first) > 1e-9 || 
               std::abs(res[i].y() - refUpper[i].second) > 1e-9){
                return false;
            }
        }
    }
    
    auto refLower = reference::lowerHull(temp);
    res = CH.LowerHullPoints();
    
    // Skip rightmost point in CHTree lower hull for comparison
    if (res.size() > 1 && refLower.size() > 1) {
        for(size_t i = 0; i < res.size() - 1 && i < refLower.size() - 1; i++){
            if(std::abs(res[i].x() - refLower[i].first) > 1e-9 || 
               std::abs(res[i].y() - refLower[i].second) > 1e-9){
                return false;
            }
        }
    }
    
    return true;
}

bool verificationTest(int verify_step, bool shuffle) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::vector<std::pair<double,double>> data;
    auto CH = CHTree<K>();
    double x;
    double y;
    for(int i = 0; i <= verify_step; i++){
        x = g()%1000;
        y = g()%1000;
        CH.Insert(Point_2(x,y));
        data.emplace_back(x,y);
    }
    if (!verify(CH, data, data.size())) {
        std::cout << "Initial verification failed\n";
        return false;
    }
    if(shuffle) {
        std::shuffle(data.begin(), data.end(), g);
    }
    for(int i = data.size()-1; i > 0 && i > (int)data.size()-verify_step-1; i--){
        CH.Remove(Point_2(data[i].first, data[i].second));
        if (!verify(CH, data, i)) {
            std::cout << "Verification failed after remove " << i << "\n";
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    int size = 10000;
    int mode = 0;
    int seed = 0;
    int verify_step = 0;
    bool shuffle = false;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            size = std::stoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-mode") == 0 && i + 1 < argc) {
            mode = std::stoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-seed") == 0 && i + 1 < argc) {
            seed = std::stoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-verify") == 0 && i + 1 < argc) {
            verify_step = std::stoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-shuffle") == 0) {
            shuffle = true;
        }
    }

    if (verify_step > 0) {
        std::cout << "Running verification against reference convex hull: ";
        if (verificationTest(verify_step, shuffle)) {
            std::cout << "PASSED\n";
        } else {
            std::cout << "FAILED\n";
            return 1;
        }
    } else {
        runtimeTest(size, mode, seed);
    }

    return 0;
}

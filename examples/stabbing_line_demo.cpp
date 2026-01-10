/**
 * @file stabbing_line_demo.cpp
 * @brief Interactive demo for StabbingLineStructure
 * 
 * Compile: g++ -std=c++17 -O2 -I../include -o stabbing_line_demo stabbing_line_demo.cpp -lCGAL -lgmp
 * Run: ./stabbing_line_demo
 */

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "StabbingLineStructure.h"
#include "inexact.h"

using K = Inexact_kernel<double>;
using Point_2 = K::Point_2;
typedef StabbingLineStructure<K> SLS;

void printResult(const SLS& sls, const std::vector<Point_2>& points, double epsilon) {
    std::cout << "  Points: ";
    for (const auto& p : points) {
        std::cout << "(" << p.x() << "," << p.y() << ") ";
    }
    std::cout << "\n  Epsilon: " << epsilon << "\n";
    
    if (auto line = sls.findStabbingLine()) {
        std::cout << "  ✓ Stabbing line EXISTS: y = " << line->slope << "x + " << line->intercept << "\n";
        
        std::cout << "  Verification:\n";
        for (const auto& p : points) {
            double lineY = line->at(p.x());
            double error = std::abs(p.y() - lineY);
            std::cout << "    Point (" << p.x() << "," << p.y() << "): line_y=" << lineY 
                      << ", error=" << error << (error <= epsilon ? " ✓" : " ✗") << "\n";
        }
    } else {
        std::cout << "  ✗ No stabbing line exists\n";
    }
    std::cout << "\n";
}

int main() {
    std::cout << "=== Epsilon-Stabbing Line Demo ===\n\n";
    
    // Example 1: Simple case - points near y = x
    std::cout << "Example 1: Points near y = x\n";
    {
        std::vector<Point_2> points = {Point_2(0, 0.3), Point_2(5, 4.8), Point_2(10, 10.2)};
        SLS sls(1.0);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 1.0);
    }
    
    // Example 2: Horizontal points
    std::cout << "Example 2: Horizontal points\n";
    {
        std::vector<Point_2> points = {Point_2(0, 5.1), Point_2(3, 5.3), Point_2(7, 4.9), Point_2(10, 5.0)};
        SLS sls(0.5);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 0.5);
    }
    
    // Example 3: Outlier point - NO stabbing line
    std::cout << "Example 3: Outlier point (NO stabbing line)\n";
    {
        std::vector<Point_2> points = {Point_2(0, 0), Point_2(5, 5), Point_2(10, 10), Point_2(5, 50)};
        SLS sls(1.0);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 1.0);
    }
    
    // Example 4: Zig-zag pattern - NO stabbing line
    std::cout << "Example 4: Zig-zag pattern (NO stabbing line)\n";
    {
        std::vector<Point_2> points = {Point_2(0, 0), Point_2(1, 20), Point_2(2, 0), Point_2(3, 20)};
        SLS sls(1.0);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 1.0);
    }
    
    // Example 5: Same x-coordinate, far apart - NO stabbing line
    std::cout << "Example 5: Same x, far apart vertically (NO stabbing line)\n";
    {
        std::vector<Point_2> points = {Point_2(5, 0), Point_2(5, 10)};
        SLS sls(1.0);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 1.0);
    }
    
    // Example 6: Same x-coordinate, close together
    std::cout << "Example 6: Same x, close together (HAS stabbing line)\n";
    {
        std::vector<Point_2> points = {Point_2(5, 0), Point_2(5, 1.5)};
        SLS sls(1.0);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 1.0);
    }
    
    // Example 7: Triangle pattern
    std::cout << "Example 7: Triangle pattern (NO stabbing line with small epsilon)\n";
    {
        std::vector<Point_2> points = {Point_2(0, 0), Point_2(5, 10), Point_2(10, 0)};
        SLS sls(0.1);
        for (const auto& p : points) sls.insert(p);
        printResult(sls, points, 0.1);
    }
    
    // Interactive mode
    std::cout << "=== Try Your Own Points ===\n";
    std::cout << "Enter epsilon value: ";
    double epsilon;
    std::cin >> epsilon;
    
    std::cout << "Enter points as (x y), one per line. Enter 'done' when finished:\n";
    
    SLS sls(epsilon);
    std::vector<Point_2> userPoints;
    std::string input;
    while (true) {
        std::cout << "> ";
        std::cin >> input;
        if (input == "done" || input == "q" || input == "quit") break;
        
        try {
            double x = std::stod(input);
            double y;
            std::cin >> y;
            sls.insert(Point_2(x, y));
            userPoints.emplace_back(x, y);
            std::cout << "  Added (" << x << ", " << y << ")\n";
        } catch (...) {
            break;
        }
    }
    
    if (!userPoints.empty()) {
        std::cout << "\nYour result:\n";
        printResult(sls, userPoints, epsilon);
    }
    
    return 0;
}

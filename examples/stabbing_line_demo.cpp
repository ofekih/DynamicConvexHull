// Simple interactive example for StabbingLineStructure
// Compile: g++ -std=c++17 -O2 -I../include -o stabbing_line_demo stabbing_line_demo.cpp -lCGAL -lgmp
// Run: ./stabbing_line_demo

#include <iostream>
#include <vector>
#include <iomanip>
#include "StabbingLineStructure.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef StabbingLineStructure<K> SLS;

void printResult(const SLS& sls, double epsilon) {
    std::cout << "  Points: ";
    for (const auto& p : sls.getOriginalPoints()) {
        std::cout << "(" << p.x() << "," << p.y() << ") ";
    }
    std::cout << "\n  Epsilon: " << epsilon << "\n";
    
    if (auto line = sls.findStabbingLine()) {
        std::cout << "  ✓ Stabbing line EXISTS: y = " << line->slope << "x + " << line->intercept << "\n";
        
        // Verify coverage
        std::cout << "  Verification:\n";
        for (const auto& p : sls.getOriginalPoints()) {
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
    
    // -------------------------------------------------------------------------
    // Example 1: Simple case - points near y = x (should have stabbing line)
    // -------------------------------------------------------------------------
    std::cout << "Example 1: Points near y = x\n";
    {
        SLS sls(1.0);
        sls.insert(Point_2(0, 0.3));
        sls.insert(Point_2(5, 4.8));
        sls.insert(Point_2(10, 10.2));
        printResult(sls, 1.0);
    }
    
    // -------------------------------------------------------------------------
    // Example 2: Horizontal points (should have horizontal stabbing line)
    // -------------------------------------------------------------------------
    std::cout << "Example 2: Horizontal points\n";
    {
        SLS sls(0.5);
        sls.insert(Point_2(0, 5.1));
        sls.insert(Point_2(3, 5.3));
        sls.insert(Point_2(7, 4.9));
        sls.insert(Point_2(10, 5.0));
        printResult(sls, 0.5);
    }
    
    // -------------------------------------------------------------------------
    // Example 3: Outlier point - NO stabbing line
    // -------------------------------------------------------------------------
    std::cout << "Example 3: Outlier point (NO stabbing line)\n";
    {
        SLS sls(1.0);
        sls.insert(Point_2(0, 0));
        sls.insert(Point_2(5, 5));
        sls.insert(Point_2(10, 10));
        sls.insert(Point_2(5, 50));  // FAR outlier
        printResult(sls, 1.0);
    }
    
    // -------------------------------------------------------------------------
    // Example 4: Zig-zag pattern - NO stabbing line
    // -------------------------------------------------------------------------
    std::cout << "Example 4: Zig-zag pattern (NO stabbing line)\n";
    {
        SLS sls(1.0);
        sls.insert(Point_2(0, 0));
        sls.insert(Point_2(1, 20));
        sls.insert(Point_2(2, 0));
        sls.insert(Point_2(3, 20));
        printResult(sls, 1.0);
    }
    
    // -------------------------------------------------------------------------
    // Example 5: Same x-coordinate, far apart - NO stabbing line
    // -------------------------------------------------------------------------
    std::cout << "Example 5: Same x, far apart vertically (NO stabbing line)\n";
    {
        SLS sls(1.0);
        sls.insert(Point_2(5, 0));
        sls.insert(Point_2(5, 10));  // 10 units apart at same x, epsilon=1
        printResult(sls, 1.0);
    }
    
    // -------------------------------------------------------------------------
    // Example 6: Same x-coordinate, close together - DOES have stabbing line
    // -------------------------------------------------------------------------
    std::cout << "Example 6: Same x, close together (HAS stabbing line)\n";
    {
        SLS sls(1.0);
        sls.insert(Point_2(5, 0));
        sls.insert(Point_2(5, 1.5));  // 1.5 apart, less than 2*epsilon
        printResult(sls, 1.0);
    }
    
    // -------------------------------------------------------------------------
    // Example 7: Triangle pattern - NO stabbing line with small epsilon
    // -------------------------------------------------------------------------
    std::cout << "Example 7: Triangle pattern (NO stabbing line with small epsilon)\n";
    {
        SLS sls(0.1);
        sls.insert(Point_2(0, 0));
        sls.insert(Point_2(5, 10));
        sls.insert(Point_2(10, 0));
        printResult(sls, 0.1);
    }
    
    // -------------------------------------------------------------------------
    // Interactive: Try your own points!
    // -------------------------------------------------------------------------
    std::cout << "=== Try Your Own Points ===\n";
    std::cout << "Enter epsilon value: ";
    double epsilon;
    std::cin >> epsilon;
    
    std::cout << "Enter points as (x y), one per line. Enter 'done' when finished:\n";
    
    SLS sls(epsilon);
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
            std::cout << "  Added (" << x << ", " << y << ")\n";
        } catch (...) {
            break;
        }
    }
    
    if (sls.size() > 0) {
        std::cout << "\nYour result:\n";
        printResult(sls, epsilon);
    }
    
    return 0;
}

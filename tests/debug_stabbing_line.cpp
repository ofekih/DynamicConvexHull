// Debug test to find minimal failing example for stabbing line algorithm
// Compile with: make debug_stabbing_line
// Run with: ./debug_stabbing_line

#define STABBING_DEBUG true

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include "StabbingLineStructure.h"
#include "inexact.h"

using namespace dch;

using K = Inexact_kernel<double>;
using Point_2 = K::Point_2;
typedef StabbingLineStructure<K> SLS;
typedef StabbingLine<K> Line;

// Independent O(n) verification
bool verifyLineCoversAllPoints(const std::vector<Point_2>& points, 
                                 double slope, double intercept, double epsilon) {
    for (const auto& p : points) {
        double lineY = slope * p.x() + intercept;
        double error = std::abs(p.y() - lineY);
        if (error > epsilon + 1e-9) {
            std::cout << "VERIFY FAIL: Point (" << p.x() << "," << p.y() 
                      << ") has error " << error << " > epsilon " << epsilon << std::endl;
            return false;
        }
    }
    return true;
}

// Reference O(n) half-plane intersection
struct HalfPlaneIntersection {
    std::vector<Point_2> points;
    double epsilon;
    
    HalfPlaneIntersection(const std::vector<Point_2>& pts, double eps) 
        : points(pts), epsilon(eps) {}
    
    std::optional<std::pair<double, double>> checkSlope(double m) const {
        double bMin = -1e308, bMax = 1e308;
        for (const auto& p : points) {
            bMin = std::max(bMin, p.y() - epsilon - m * p.x());
            bMax = std::min(bMax, p.y() + epsilon - m * p.x());
        }
        if (bMin <= bMax + 1e-9) {
            return std::make_pair(m, (bMin + bMax) / 2.0);
        }
        return std::nullopt;
    }
    
    std::optional<std::pair<double, double>> findLine() const {
        if (points.empty()) return std::make_pair(0.0, 0.0);
        
        // Try horizontal
        if (auto r = checkSlope(0.0)) return r;
        
        // Try slopes from point pairs
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                double dx = points[j].x() - points[i].x();
                if (std::abs(dx) < 1e-12) continue;
                
                double dy1 = (points[j].y() - epsilon) - (points[i].y() + epsilon);
                double dy2 = (points[j].y() + epsilon) - (points[i].y() - epsilon);
                
                if (auto r = checkSlope(dy1 / dx)) return r;
                if (auto r = checkSlope(dy2 / dx)) return r;
            }
        }
        
        return std::nullopt;
    }
    
    bool hasStabbingLine() const {
        return findLine().has_value();
    }
};

void testCase(const std::vector<Point_2>& points, double epsilon, const std::string& name) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST: " << name << std::endl;
    std::cout << "Points (" << points.size() << "): ";
    for (const auto& p : points) std::cout << "(" << p.x() << "," << p.y() << ") ";
    std::cout << "\nEpsilon: " << epsilon << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Build structure
    SLS sls(epsilon);
    std::vector<Point_2> sorted = points;
    std::sort(sorted.begin(), sorted.end(), 
        [](const Point_2& a, const Point_2& b) {
            if (a.x() != b.x()) return a.x() < b.x();
            return a.y() < b.y();
        });
    sls.Build(sorted);
    
    // Debug info was removed during cleanup
    
    // Run algorithm
    std::cout << "\n--- Running findStabbingLine() ---" << std::endl;
    auto algResult = sls.FindStabbingLine();
    
    // Run reference
    HalfPlaneIntersection ref(points, epsilon);
    auto refResult = ref.findLine();
    
    std::cout << "\n--- Results ---" << std::endl;
    std::cout << "Algorithm: " << (algResult ? "YES" : "NO");
    if (algResult) {
        std::cout << " (y = " << algResult->slope << "x + " << algResult->intercept << ")";
    }
    std::cout << std::endl;
    
    std::cout << "Reference: " << (refResult ? "YES" : "NO");
    if (refResult) {
        std::cout << " (y = " << refResult->first << "x + " << refResult->second << ")";
    }
    std::cout << std::endl;
    
    // Verify
    bool algValid = true;
    if (algResult) {
        algValid = verifyLineCoversAllPoints(points, algResult->slope, algResult->intercept, epsilon);
        std::cout << "Algorithm line covers all points: " << (algValid ? "YES" : "NO") << std::endl;
    }
    
    // Check for discrepancy
    bool hasDiscrepancy = false;
    if (algResult.has_value() != refResult.has_value()) {
        hasDiscrepancy = true;
        std::cout << "*** DISCREPANCY: Algorithm and reference disagree on existence! ***" << std::endl;
    } else if (algResult && !algValid) {
        hasDiscrepancy = true;
        std::cout << "*** DISCREPANCY: Algorithm returned invalid line! ***" << std::endl;
    }
    
    if (hasDiscrepancy) {
        std::cout << "\n*** FAILURE ***" << std::endl;
    } else {
        std::cout << "\nPASS" << std::endl;
    }
}

void testTrial29() {
    // Reproduce trial 29 from ReturnedLineAlwaysValid test
    std::cout << "\n=== Reproducing Trial 29 ===\n" << std::endl;
    
    std::mt19937 rng(42);  // Same seed as test
    std::uniform_real_distribution<double> noise(-0.8, 0.8);
    
    // Fast-forward to trial 29
    for (int trial = 0; trial < 29; ++trial) {
        (void)rng();  // slope
        (void)rng();  // intercept
        int n = 5 + (rng() % 20);
        for (int i = 0; i < n; ++i) {
            (void)noise(rng);
        }
    }
    
    // Now generate trial 29
    double slope = ((int)(rng() % 100) - 50) / 10.0;
    double intercept = ((int)(rng() % 100) - 50);
    int n = 5 + (rng() % 20);
    double epsilon = 1.0;
    
    std::vector<Point_2> points;
    for (int i = 0; i < n; ++i) {
        double x = i * 2.0;
        double y = slope * x + intercept + noise(rng);
        points.emplace_back(x, y);
    }
    
    testCase(points, epsilon, "Trial 29 (seed 42)");
}

void findMinimalFailingCase() {
    std::cout << "\n=== Finding minimal failing case ===\n" << std::endl;
    
    // Try increasingly large random tests with different seeds
    for (int seed = 0; seed < 1000; ++seed) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> noise(-0.5, 0.5);
        std::uniform_real_distribution<double> wild(-50, 50);
        std::uniform_real_distribution<double> epsDist(0.5, 1.5);
        std::uniform_int_distribution<int> nDist(2, 5);  // Start with very small n
        std::uniform_int_distribution<int> slopeDist(-50, 49);
        std::uniform_int_distribution<int> interceptDist(-50, 49);
        
        std::vector<Point_2> points;
        double epsilon = epsDist(rng);
        int n = nDist(rng);
        
        double slope = slopeDist(rng) / 10.0;
        double intercept = interceptDist(rng);
        bool addNoise = rng() % 2 == 0;
        
        for (int i = 0; i < n; ++i) {
            double x = i * 2.0;
            double y = addNoise ? (slope * x + intercept + noise(rng)) : wild(rng);
            points.emplace_back(x, y);
        }
        
        // Sort for build
        std::vector<Point_2> sorted = points;
        std::sort(sorted.begin(), sorted.end(), 
            [](const Point_2& a, const Point_2& b) {
                if (a.x() != b.x()) return a.x() < b.x();
                return a.y() < b.y();
            });
        
        // Quick check without debug
        SLS sls(epsilon);
        sls.Build(sorted);
        auto algResult = sls.FindStabbingLine();
        
        HalfPlaneIntersection ref(points, epsilon);
        auto refResult = ref.findLine();
        
        bool algValid = true;
        if (algResult) {
            algValid = verifyLineCoversAllPoints(points, algResult->slope, algResult->intercept, epsilon);
        }
        
        bool hasFailure = (algResult.has_value() != refResult.has_value()) || (algResult && !algValid);
        
        if (hasFailure) {
            std::cout << "*** Found failing case at seed " << seed << " ***" << std::endl;
            testCase(points, epsilon, "Seed " + std::to_string(seed));
            return;
        }
    }
    
    std::cout << "No failures found in first 1000 seeds with small n" << std::endl;
}

int main() {
    std::cout << "Stabbing Line Debug Test" << std::endl;
    std::cout << "========================" << std::endl;
    
    // First, test some known simple cases
    testCase({{0, 0}, {5, 5}, {10, 10}}, 1.0, "Collinear points");
    testCase({{0, 0}, {5, 5}}, 1.0, "Two points");
    testCase({{0, 0}, {5, 5}, {5, 50}}, 1.0, "Outlier (should fail)");
    
    // Test the specific failing trial
    testTrial29();
    
    // Find minimal failing case
    findMinimalFailingCase();
    
    return 0;
}


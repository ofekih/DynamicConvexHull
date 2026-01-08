add_executable(test_convex_hull test_convex_hull.cpp)
target_link_libraries(test_convex_hull PRIVATE dch_lib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_convex_hull)

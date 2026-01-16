cmake_minimum_required(VERSION 3.22)
project(DynamicConvexHull)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3")

include(FetchContent)

# Google Test
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Google Benchmark
FetchContent_Declare(
  benchmark
  URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.zip
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(benchmark)

# Create an interface library for the headers (no external dependencies)
add_library(dch_lib INTERFACE)
target_include_directories(dch_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(DynamicConvexHull main.cpp)
target_link_libraries(DynamicConvexHull dch_lib)

enable_testing()

add_subdirectory(tests)
add_subdirectory(benchmarks)
add_subdirectory(examples)
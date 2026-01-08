.PHONY: all build test benchmark clean

# Default target
all: build

# Configure and build the project
build:
	cmake -S . -B build
	cmake --build build

# Run unit tests
test: build
	ctest --test-dir build --output-on-failure

# Run benchmarks
benchmark: build
	./build/benchmarks/benchmark_convex_hull

# Clean build directory
clean:
	rm -rf build

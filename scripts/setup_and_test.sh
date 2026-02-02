#!/bin/bash
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_ROOT/build"

# Limit parallel jobs to avoid memory overflow in WSL
# Set to number of cores minus 1, or 2 for low-memory systems
PARALLEL_JOBS="${PARALLEL_JOBS:-2}"

echo -e "${GREEN}=== MKF Build and Test Script ===${NC}"
echo "Project root: $PROJECT_ROOT"
echo "Build directory: $BUILD_DIR"

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install dependencies
install_deps() {
    echo -e "\n${YELLOW}Installing build dependencies...${NC}"

    if command_exists apt; then
        sudo apt update
        sudo apt install -y cmake ninja-build libopenblas-dev gnuplot
    else
        echo -e "${RED}apt not found. Please install cmake, ninja-build, libopenblas-dev, gnuplot manually.${NC}"
        exit 1
    fi

    if command_exists npm; then
        echo -e "${YELLOW}Installing quicktype...${NC}"
        npm install -g quicktype || sudo npm install -g quicktype
    else
        echo -e "${RED}npm not found. Please install Node.js and npm first.${NC}"
        exit 1
    fi
}

# Function to verify dependencies
check_deps() {
    echo -e "\n${YELLOW}Checking dependencies...${NC}"
    local missing=0

    for cmd in cmake ninja g++ node npm quicktype; do
        if command_exists "$cmd"; then
            echo -e "  ${GREEN}✓${NC} $cmd found"
        else
            echo -e "  ${RED}✗${NC} $cmd missing"
            missing=1
        fi
    done

    return $missing
}

# Function to configure and build
build_project() {
    echo -e "\n${YELLOW}Configuring CMake...${NC}"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    cmake .. -G "Ninja" \
        -DCMAKE_BUILD_TYPE=Release \
        -DINCLUDE_MKF_TESTS=ON \
        -DEMBED_MAS_DATA=ON

    echo -e "\n${YELLOW}Building project (parallel jobs: $PARALLEL_JOBS)...${NC}"
    ninja -j"$PARALLEL_JOBS"

    echo -e "${GREEN}Build complete!${NC}"
}

# Function to run tests
run_tests() {
    local test_filter="${1:-}"

    cd "$BUILD_DIR"

    if [[ ! -f "MKF_tests" ]]; then
        echo -e "${RED}MKF_tests executable not found. Run build first.${NC}"
        exit 1
    fi

    if [[ -n "$test_filter" ]]; then
        echo -e "\n${YELLOW}Running tests with filter: $test_filter${NC}"
        ./MKF_tests "$test_filter" --reporter console
    else
        echo -e "\n${YELLOW}Running smoke tests...${NC}"
        ./MKF_tests "[smoke-test]" --reporter console
    fi
}

# Function to run all tests
run_all_tests() {
    cd "$BUILD_DIR"

    if [[ ! -f "MKF_tests" ]]; then
        echo -e "${RED}MKF_tests executable not found. Run build first.${NC}"
        exit 1
    fi

    echo -e "\n${YELLOW}Running all tests...${NC}"
    ./MKF_tests --reporter console
}

# Function to list available tests
list_tests() {
    cd "$BUILD_DIR"

    if [[ ! -f "MKF_tests" ]]; then
        echo -e "${RED}MKF_tests executable not found. Run build first.${NC}"
        exit 1
    fi

    echo -e "\n${YELLOW}Available tests:${NC}"
    ./MKF_tests --list-tests
}

# Print usage
usage() {
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  install     Install build dependencies (requires sudo)"
    echo "  check       Check if all dependencies are installed"
    echo "  build       Configure and build the project"
    echo "  test        Run smoke tests (default)"
    echo "  test-all    Run all tests"
    echo "  test-filter <filter>  Run tests matching filter (e.g., '[adviser]')"
    echo "  list        List all available tests"
    echo "  full        Full setup: install deps, build, and run smoke tests"
    echo "  rebuild     Clean and rebuild"
    echo ""
    echo "Environment variables:"
    echo "  PARALLEL_JOBS  Number of parallel compile jobs (default: 2, for WSL memory limits)"
    echo ""
    echo "Examples:"
    echo "  $0 full                    # Complete setup and test"
    echo "  $0 build && $0 test        # Build then test"
    echo "  $0 test-filter '[adviser]' # Run adviser tests only"
    echo "  $0 test-filter 'Settings'  # Run tests with 'Settings' in name"
    echo "  PARALLEL_JOBS=4 $0 build   # Build with 4 parallel jobs"
}

# Main logic
case "${1:-test}" in
    install)
        install_deps
        ;;
    check)
        check_deps
        ;;
    build)
        check_deps || { echo -e "${RED}Missing dependencies. Run '$0 install' first.${NC}"; exit 1; }
        build_project
        ;;
    test)
        run_tests "[smoke-test]"
        ;;
    test-all)
        run_all_tests
        ;;
    test-filter)
        if [[ -z "${2:-}" ]]; then
            echo -e "${RED}Please provide a test filter.${NC}"
            usage
            exit 1
        fi
        run_tests "$2"
        ;;
    list)
        list_tests
        ;;
    full)
        install_deps
        check_deps
        build_project
        run_tests "[smoke-test]"
        echo -e "\n${GREEN}=== Setup complete! ===${NC}"
        ;;
    rebuild)
        echo -e "${YELLOW}Cleaning build directory...${NC}"
        rm -rf "$BUILD_DIR"
        build_project
        ;;
    help|--help|-h)
        usage
        ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        usage
        exit 1
        ;;
esac

#!/bin/bash
# Build script for GCTB Python interface

set -e  # Exit on error

echo "=========================================="
echo "Building GCTB Python Interface"
echo "=========================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if running on macOS or Linux
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo -e "${GREEN}Detected macOS${NC}"
    SYSTEM="macOS"
else
    echo -e "${GREEN}Detected Linux${NC}"
    SYSTEM="Linux"
fi

# Check for required tools
echo ""
echo "Checking prerequisites..."

command -v cmake >/dev/null 2>&1 || { 
    echo -e "${RED}Error: cmake is required but not installed.${NC}" >&2
    exit 1
}
echo -e "${GREEN}✓ CMake found${NC}"

command -v python3 >/dev/null 2>&1 || { 
    echo -e "${RED}Error: python3 is required but not installed.${NC}" >&2
    exit 1
}
echo -e "${GREEN}✓ Python3 found${NC}"

# Check for pybind11
python3 -c "import pybind11" 2>/dev/null || {
    echo -e "${YELLOW}Warning: pybind11 not found. Installing...${NC}"
    pip3 install pybind11
}
echo -e "${GREEN}✓ pybind11 found${NC}"

# Set environment variables for Eigen and Boost if not set
if [ -z "$EIGEN3_INCLUDE_DIR" ]; then
    if [[ "$SYSTEM" == "macOS" ]]; then
        export EIGEN3_INCLUDE_DIR="/usr/local/include/eigen3"
    else
        export EIGEN3_INCLUDE_DIR="/usr/include/eigen3"
    fi
fi

if [ -z "$BOOST_LIB" ]; then
    if [[ "$SYSTEM" == "macOS" ]]; then
        export BOOST_LIB="/usr/local/include"
    else
        export BOOST_LIB="/usr/include"
    fi
fi

echo ""
echo "Environment:"
echo "  EIGEN3_INCLUDE_DIR: $EIGEN3_INCLUDE_DIR"
echo "  BOOST_LIB: $BOOST_LIB"

# Clean old build
echo ""
echo "Cleaning old build files..."
rm -rf build/
rm -rf python/gctb/_core*.so
rm -rf python/gctb.egg-info/

# Install in development mode
echo ""
echo -e "${GREEN}Building and installing GCTB...${NC}"
pip3 install -e . -v

# Run tests
echo ""
echo -e "${GREEN}Running basic tests...${NC}"
python3 -c "import gctb; print(f'GCTB version: {gctb.__version__}')" || {
    echo -e "${RED}Failed to import gctb${NC}"
    exit 1
}

echo ""
echo -e "${GREEN}=========================================="
echo "Build completed successfully!"
echo "==========================================${NC}"
echo ""
echo "Try it out:"
echo "  python3 -c 'import gctb; data = gctb.Data(); print(data)'"
echo ""
echo "Run tests:"
echo "  cd python && pytest tests/ -v"


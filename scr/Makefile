# Makefile for gctb
# Supported platforms: UNIX/Linux, macOS

# Set the output filename
OUTPUT = gctb

# Set the path to dependent libraries
EIGEN_PATH = $(EIGEN3_INCLUDE_DIR)
BOOST_PATH = $(BOOST_LIB)

# Set the C++ compiler based on the system
ifeq ($(shell uname -s),Linux)
    CXX = g++
    CXXFLAGS = -I $(EIGEN_PATH) -I $(BOOST_PATH) -DUNIX -DNDEBUG -msse2 -m64 -fopenmp -O3 -std=c++0x -Wall -w
    LIB = -ldl
else ifeq ($(shell uname -s),Darwin)
    CXX = /usr/local/opt/llvm/bin/clang++
    CXXFLAGS = -I /usr/local/include -I $(EIGEN_PATH) -I $(BOOST_PATH) -DUNIX -Dfopen64=fopen -stdlib=libc++ -m64 -fopenmp -O3 -std=c++0x -Wall -w
    LIB = -lz -lm -L/usr/local/opt/libomp/lib
endif

# Set the source files and headers
HDR = gctb.hpp data.hpp gadgets.hpp hsq.hpp mcmc.hpp model.hpp options.hpp stat.hpp vgmaf.hpp predict.hpp xci.hpp stratify.hpp
SRC = gctb.cpp data.cpp gadgets.cpp hsq.cpp main.cpp mcmc.cpp model.cpp options.cpp stat.cpp vgmaf.cpp predict.cpp xci.cpp stratify.cpp eigen.cpp
OBJ = $(SRC:.cpp=.o)

# Define the rules
.PHONY: all clean

all: $(OUTPUT)

$(OUTPUT): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB)

%.o: %.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(OUTPUT) *~

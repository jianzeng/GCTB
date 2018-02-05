
# -----------------------------------------------------------------
#   Makefile for gctb 
#   
#   Supported platforms
#       Unix / Linux                UNIX
#       MacOS                       MAC
# ---------------------------------------------------------------------

# Set this variable to either LINUX, MAC or WIN
SYS = LINUX
OUTPUT = gctb

# Put C++ compiler here
CXX_UNIX = /Users/l.lloydjones/homebrew/Cellar/llvm/5.0.1/bin/clang++   #g++
CXX_MAC = clang-omp   #g++

# Any other compiler flags here ( -Wall, -g, etc)
CXXFLAGS = -w -O3 -m64 -DEIGEN_NO_DEBUG -msse2 -std=c++0x -I.

# Some system specific flags
ifeq ($(SYS),LINUX)
 CXXFLAGS += -fopenmp -I/Users/l.lloydjones/homebrew/Cellar/eigen/3.3.4/include/eigen3/ -I/Users/l.lloydjones/homebrew/Cellar/boost/1.62.0/include -I /Users/l.lloydjones/homebrew/Cellar/open-mpi/3.0.0_2/include  -DUNIX -Dfopen64=fopen
 LIB += -L /Users/l.lloydjones/homebrew/Cellar/open-mpi/3.0.0_2/lib -L /Users/l.lloydjones/homebrew/Cellar/llvm/5.0.1/lib -lz -lmpi -lm
 CXX = $(CXX_UNIX)
endif

ifeq ($(SYS),MAC)
 CXXFLAGS += -fopenmp -I /usr/local/include -I /usr/local/lib -I /Users/j.zeng/XCode/Eigen3 -I /Users/j.zeng/XCode -I /usr/local/include/libiomp -DUNIX -Dfopen64=fopen 
 LIB += -lz -lmpi -lm
 CXX = $(CXX_MAC)
endif


HDR = gctb.hpp \
      data.hpp \
      gadgets.hpp \
      hsq.hpp \
      mcmc.hpp \
      model.hpp \
      options.hpp \
      predict.hpp \
      stat.hpp \
      mympi.hpp \
      xci.hpp \
      vgmaf.hpp

SRC = gctb.cpp \
      data.cpp \
      gadgets.cpp \
      hsq.cpp \
      main.cpp \
      mcmc.cpp \
      model.cpp \
      options.cpp \
      predict.cpp \
      stat.cpp \
      mympi.cpp \
      xci.cpp \
      vgmaf.cpp

	   
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	rm -f *.o *~

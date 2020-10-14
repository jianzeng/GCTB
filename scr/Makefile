
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

# Set path to the dependent libraries
EIGEN = $(EIGEN3_INCLUDE_DIR)
BOOST = $(BOOST_LIB)

# Put C++ compiler here
CXX_UNIX = g++
CXX_MAC = g++

# Any other compiler flags here ( -Wall, -g, etc)
CXXFLAGS = -fopenmp -w -O3 -m64 -DEIGEN_NO_DEBUG -msse2 -std=c++0x -I.

# Some system specific flags
ifeq ($(SYS),LINUX)
 CXXFLAGS += -I $(EIGEN) -I $(BOOST) -DUNIX -m64 -ldl
 #LIB += -static-libgcc -static-libstdc++
 LIB += -static
# LIB += -static -Wl,--start-group -Wl,--end-group -lpthread -lm
 CXX = $(CXX_UNIX)
endif

ifeq ($(SYS),MAC)
 CXXFLAGS += -I /usr/local/include -I $(EIGEN) -I $(BOOST) -DUNIX -Dfopen64=fopen 
 LIB += -lz -lm
 CXX = $(CXX_MAC)
endif


HDR = gctb.hpp \
      data.hpp \
      gadgets.hpp \
      hsq.hpp \
      mcmc.hpp \
      model.hpp \
      options.hpp \
      stat.hpp \
      vgmaf.hpp \
      predict.hpp \
      xci.hpp \
      stratify.hpp

SRC = gctb.cpp \
      data.cpp \
      gadgets.cpp \
      hsq.cpp \
      main.cpp \
      mcmc.cpp \
      model.cpp \
      options.cpp \
      stat.cpp \
      vgmaf.cpp \
      predict.cpp \
      xci.cpp \
      stratify.cpp

	   
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

# Set this variable to either LINUX, MAC or WIN
OUTPUT = greml

# Set path to library dependencies
EIGEN = ./eigen-3.4.0/

# Put C++ compiler here
CXX = g++

# Any other compiler flags here ( -Wall, -g, etc) -march=native
# Other back-end for eigen  -llapack -llapacke -lblas
# -std=c++11 used before and worked instead of -std=gnu++11

#CXXFLAGS   = -fopenmp -Wall -std=c++11 -g -Ofast -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE -I $(EIGEN) -llapack -llapacke -lblas -march=native
CXXFLAGS   = -fopenmp -Wall -std=c++11 -g -O3 -I $(EIGEN) -march=native

LIB += -static-libgcc

HDR = REML_f.hpp \
	REML_d.hpp

SRC = greml.cpp \
	REML_f.cpp \
	REML_d.cpp

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
	rm -f *.o *~ greml

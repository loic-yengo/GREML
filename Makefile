#ldscdom: 	ldscdom.cpp LDSC.h
#	g++ -I"./eigen/" -Wall -std=c++11 -g -O3 ldscdom.cpp -o ldscdom -lm	
#clean:
#	rm ldscdom
# Set this variable to either LINUX, MAC or WIN
SYS = LINUX
OUTPUT = greml

# Set path to library dependencies
EIGEN = ./eigen-3.4.0/

# Put C++ compiler here
CXX = g++

# Any other compiler flags here ( -Wall, -g, etc)
CXXFLAGS = -Wall -std=c++11 -g -O3 -DNDEBUG -I. -march=native

# Some system specific flags
ifeq ($(SYS),LINUX)
 CXXFLAGS += -I $(EIGEN)
endif

HDR = REML.hpp
SRC = greml.cpp \
			REML.cpp

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

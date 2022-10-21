# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: Oct 15, 2022
#
# This makefile creates exception objects in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++11 -c    

OPTIM  = -O3

DEBUG_FLAG  = -g
# The directory where the include files needed to create the library objects are
INC_DIR = ../../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib

# targets this makefile builds
EXP_OBJS = baseexception.o

DEBUG_OBJS = debug_baseexception.o

baseexception.o: $(INC_DIR)/baseexception.h 
	$(CXX) $(CXXFLAGS) baseexception.cxx -I$(INC_DIR)
	mv baseexception.o $(LIB_DIR)

debug_baseexception.o: $(INC_DIR)/baseexception.h 
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) baseexception.cxx -I$(INC_DIR)
	mv baseexception.o $(LIB_DIR)/debug_baseexception.o

all: $(EXP_OBJS)

debug: $(DEBUG_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(EXP_OBJS); rm $(DEBUG_OBJS)



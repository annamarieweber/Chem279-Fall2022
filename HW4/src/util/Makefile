# MSSE Chem 279- Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: November 5 2022
#
# This makefile creates objects from utility functions in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++11 -c 

OPTIM       = -O3

DEBUG_FLAG  = -g

# The directory where the include files needed to create the library objects are
INC_DIR = ../../include
# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib

# The next line contains the list of utility objects created by this Makefile.
UTILS_OBJS = filereader.o combinatorics.o

DEBUG_OBJS = debug_filereader.o debug_combinatorics.o

filereader.o: $(INC_DIR)/filereader.h
	$(CXX) $(CXXFLAGS) $(OPTIM) filereader.cxx -I$(INC_DIR)
	mv filereader.o $(LIB_DIR)

combinatorics.o: $(INC_DIR)/combinatorics.h
	$(CXX) $(CXXFLAGS) $(OPTIM) combinatorics.cxx -I$(INC_DIR)
	mv combinatorics.o $(LIB_DIR)

debug_filereader.o: $(INC_DIR)/filereader.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) filereader.cxx -I$(INC_DIR)
	mv filereader.o $(LIB_DIR)/debug_filereader.o

debug_combinatorics.o: $(INC_DIR)/combinatorics.h
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) combinatorics.cxx -I$(INC_DIR)
	mv combinatorics.o $(LIB_DIR)/debug_combinatorics.o

debug: $(DEBUG_OBJS)

all: $(UTILS_OBJS)

cleanall:
	cd $(LIB_DIR); rm $(UTILS_OBJS); rm $(DEBUG_OBJS)


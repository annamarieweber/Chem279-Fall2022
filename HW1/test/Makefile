# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: Sept 18 2022
#
# This makefile creates test objects in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object must use c++17 for sciplot
CXXFLAGS    = -std=c++17 -O3 -c    

# flags for building wtih the pyhton library to ue matplotc++

# The directory where the include files needed to create the library objects are
INC_DIR = ../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../lib

# targets this makefile builds
TEST_OBJS = evaluation.o

evaluation.o: $(INC_DIR)/evaluation.h 
	$(CXX) $(CXXFLAGS) evaluation.cxx -I$(INC_DIR) 
	mv evaluation.o $(LIB_DIR)

all: $(TEST_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(TEST_OBJS)



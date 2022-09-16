# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: Sept 10, 2022
#
# This makefile creates molecules objects in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++11 -O3 -c    

# The directory where the include files needed to create the library objects are
INC_DIR = ../../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib

# targets this makefile builds
CLUSTER_OBJS = cluster.o

cluster.o: $(INC_DIR)/cluster.h 
	$(CXX) $(CXXFLAGS) cluster.cxx -I$(INC_DIR)
	mv cluster.o $(LIB_DIR)

all: $(CLUSTER_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(CLUSTER_OBJS)



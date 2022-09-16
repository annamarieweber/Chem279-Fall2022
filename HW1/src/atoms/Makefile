# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: Sept 10, 2022
#
# This makefile creates atoms objects in C++
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
ATOM_OBJS = atomcoord.o atominfos.o clustermap.o

atomcoord.o: $(INC_DIR)/atomcoord.h 
	$(CXX) $(CXXFLAGS) atomcoord.cxx -I$(INC_DIR)
	mv atomcoord.o $(LIB_DIR)

atominfos.o: $(INC_DIR)/atominfos.h $(INC_DIR)/atomcoord.h
	$(CXX) $(CXXFLAGS) atominfos.cxx -I$(INC_DIR)
	mv atominfos.o $(LIB_DIR)

clustermap.o: $(INC_DIR)/atominfos.h $(INC_DIR)/clustermap.h $(INC_DIR)/constants.h $(INC_DIR)/baseexception.h
	$(CXX) $(CXXFLAGS) clustermap.cxx -I$(INC_DIR)
	mv clustermap.o $(LIB_DIR)

all: $(ATOM_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(ATOM_OBJS)


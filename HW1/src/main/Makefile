# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator: Anna Weber
# Date Created: Sept 10, 2022
#
# This makefile creates the main executable and runs it

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++11 -O3 -o

# armadillo flags
ARMAFLAGS = -larmadillo

# The directory where the include files needed to create the library objects are
INC_DIR = ../../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib
BIN_DIR     = ../../bin
OUTPUT_DIR  = ../../output

# the lib objects to include
INC_OBJS = $(INC_DIR) $(LIB_DIR)/filereader.o $(LIB_DIR)/cluster.o $(LIB_DIR)/clustermap.o $(LIB_DIR)/atominfos.o $(LIB_DIR)/atomcoord.o $(LIB_DIR)/baseexception.o $(LIB_DIR)/evaluation.o

# The next line contains the list of utility objects created by this Makefile.
EXECUTABLES = main_program

main_program: $(INC_DIR)/filereader.h $(INC_DIR)/cluster.h $(INC_DIR)/clustermap.h $(INC_DIR)/atominfos.h $(INC_DIR)/atomcoord.h $(INC_DIR)/baseexception.h $(INC_DIR)/evaluation.h
	$(CXX) $(CXXFLAGS) main_program main.cxx -I$(INC_OBJS) $(ARMAFLAGS)
	mv main_program $(BIN_DIR)

all: $(EXECUTABLES)

run:
	make all
	echo Running Program ....
	$(BIN_DIR)/main_program ../../sample-input.txt | $(OUTPUT_DIR)/main_program.out


cleanall:
	cd $(BIN_DIR); rm $(EXECUTABLES)


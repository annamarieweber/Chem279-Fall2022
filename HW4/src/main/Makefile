# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator: Anna Weber
# Date Created: November 5, 2022
#
# This makefile creates the main executable and runs it

# GNU C++ Compiler
CXX         = g++      

# compiler flags, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++11 -o

OPTIM       = -O3

DEBUG_FLAG   = -g

# armadillo flags
ARMAFLAGS = -larmadillo

# The directory where the include files needed to create the library objects are
INC_DIR = ../../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../../lib
BIN_DIR     = ../../bin
OUTPUT_DIR  = ../../output

# the lib objects to include
INC_OBJS = $(INC_DIR) $(LIB_DIR)/filereader.o $(LIB_DIR)/combinatorics.o $(LIB_DIR)/baseexception.o $(LIB_DIR)/shell.o $(LIB_DIR)/shelloverlapintegral.o $(LIB_DIR)/cluster.o 

DEBUG_OBJS = $(INC_DIR) $(LIB_DIR)/debug_filereader.o $(LIB_DIR)/debug_combinatorics.o $(LIB_DIR)/debug_baseexception.o $(LIB_DIR)/debug_shell.o $(LIB_DIR)/debug_shelloverlapintegral.o $(LIB_DIR)/debug_cluster.o 

# The next line contains the list of utility objects created by this Makefile.
EXECUTABLES = main_program

main_program: $(INC_DIR)/filereader.h $(INC_DIR)/combinatorics.h $(INC_DIR)/baseexception.h $(INC_DIR)/shell.h $(INC_DIR)/shelloverlapintegral.h $(INC_DIR)/cluster.h 
	$(CXX) $(CXXFLAGS) main_program $(OPTIM) main.cxx -I$(INC_OBJS) $(ARMAFLAGS)
	mv main_program $(BIN_DIR)

all: $(EXECUTABLES)

debug:$(INC_DIR)/filereader.h $(INC_DIR)/combinatorics.h $(INC_DIR)/baseexception.h $(INC_DIR)/shell.h $(INC_DIR)/shelloverlapintegral.h $(INC_DIR)/cluster.h 
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) main_program main.cxx -I$(DEBUG_OBJS) $(ARMAFLAGS)
	rm -rf $(BIN_DIR)/*
	mv main_program* $(BIN_DIR)
	

run:
	make all
	echo Running Program ....
	$(BIN_DIR)/main_program ../../C2H2.txt | tee $(OUTPUT_DIR)/main_program.out


cleanall:
	cd $(BIN_DIR); rm $(EXECUTABLES)


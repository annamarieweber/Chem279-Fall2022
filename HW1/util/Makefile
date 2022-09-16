# MSSE Chem 279- Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: 09/10/2022
#
# This makefile creates objects from utility functions in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -std=c++11 -O3 -c 

# The directory where the include files needed to create the library objects are
INC_DIR = ../include
# The directory where all object files are going to moved afer creation
LIB_DIR     = ../lib

# The next line contains the list of utility objects created by this Makefile.
UTILS_OBJS = filereader.o

filereader.o: $(INC_DIR)/filereader.h
	$(CXX) $(CXXFLAGS) filereader.cxx -I$(INC_DIR)
	mv filereader.o $(LIB_DIR)

all: $(UTILS_OBJS)

cleanall:
	cd $(LIB_DIR); rm $(UTILS_OBJS)


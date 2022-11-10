# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: November 05, 2022
#
# This makefile creates exception objects in C++
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

# targets this makefile builds
SHELL_OBJS = shell.o shelloverlapintegral.o

DEBUG_OBJS = debug_shell.o debug_shelloverlapintegral.o

shell.o: $(INC_DIR)/shell.h 
	$(CXX) $(CXXFLAGS) $(OPTIM) shell.cxx -I$(INC_DIR)
	mv shell.o $(LIB_DIR)

shelloverlapintegral.o: $(INC_DIR)/shell.h $(INC_DIR)/shelloverlapintegral.h 
	$(CXX) $(CXXFLAGS) $(OPTIM) shelloverlapintegral.cxx -I$(INC_DIR)
	mv shelloverlapintegral.o $(LIB_DIR)

debug_shell.o: $(INC_DIR)/shell.h 
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) shell.cxx -I$(INC_DIR)
	mv shell.o $(LIB_DIR)/debug_shell.o

debug_shelloverlapintegral.o: $(INC_DIR)/shell.h $(INC_DIR)/shelloverlapintegral.h 
	$(CXX) $(DEBUG_FLAG) $(CXXFLAGS) shelloverlapintegral.cxx -I$(INC_DIR)
	mv shelloverlapintegral.o $(LIB_DIR)/debug_shelloverlapintegral.o

all: $(SHELL_OBJS)

debug: $(DEBUG_OBJS)

cleanall: 
	cd $(LIB_DIR); rm $(SHELL_OBJS); rm$(DEBUG_OBJS)



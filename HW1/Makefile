# This Makefile generates the Components Needed for HW1 
# for MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  Anna Weber
# Date Created: Sept 04, 2022

# This is the main build makefile

TEST_DIR = tests/
ATOMS_DIR = src/atoms
MOL_DIR = src/molecules
EXCEPTIONS = src/exceptions
MAIN_DIR = src/main
UTIL_DIR = util/

all:
	cd $(EXCEPTIONS); make all
	cd $(UTIL_DIR); make all
	cd $(MOL_DIR); make all
	cd $(ATOMS_DIR); make all
	cd $(MAIN_DIR); make all

testall:
	make all
	cd $(TEST_DIR); make test

run:
	cd $(MAIN_DIR); make run

cleanall:
	cd $(ATOMS_DIR); make cleanall
	cd $(MOL_DIR); make cleanall
	cd $(UTIL_DIR); make cleanall
	cd $(MAIN_DIR); make cleanall
	cd $(TEST_DIR); make cleanall



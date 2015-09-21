# LIBFORBES MAKEFILE

#
# Configurable properties
#

# C++ compiler
CXX = g++

# Additional compiler flags (e.g., -O2 or -O3 optimization flags, etc)
# To create a test coverage report add: -fprofile-arcs -ftest-coverage
CFLAGS_ADDITIONAL = -O3 -fprofile-arcs -ftest-coverage

# Additional link flags
# To create a test coverage report add: -fprofile-arcs
LFLAGS_ADDITIONAL = -fprofile-arcs

# Lapacke Include directory on your system
LAPACKE_INCLUDE = /usr/include/lapacke

include Configuration

#                         #
# Do not modify this file #
#                         #
OBJ_DIR = build/Debug
BIN_DIR = dist/Debug

OBJ_TEST_DIR = ${OBJ_DIR}/Test
BIN_TEST_DIR = ${BIN_DIR}/Test

TEST_DIR = source/tests

ARCHIVE = ${BIN_DIR}/libforbes.a

CFLAGS = -c -DUSE_LIBS $(CFLAGS_ADDITIONAL)

IFLAGS = \
	 -I. \
	 -I./source \
	 -I$(SS_DIR)/CHOLMOD/Include \
	 -I$(SS_DIR)/LDL/Include \
	 -I$(SS_DIR)/SuiteSparse_config \
	 -I$(LAPACKE_INCLUDE)

lFLAGS = \
	-lldl \
	-lcholmod \
	-lamd \
	-lcolamd \
	-lsuitesparseconfig \
	-lccolamd \
	-lcamd \
	-llapacke \
	-lblas \
	-llapack \
	-lopenblas \
	-lm -lrt \
	$(LFLAGS_ADDITIONAL)

LFLAGS = -L$(SS_DIR)/CHOLMOD/Lib \
	 -L$(SS_DIR)/AMD/Lib \
	 -L$(SS_DIR)/COLAMD/Lib \
	 -L$(SS_DIR)/SuiteSparse_config \
	 -L$(SS_DIR)/CCOLAMD/Lib \
	 -L$(SS_DIR)/CAMD/Lib \
	 -L$(SS_DIR)/LDL/Lib

SOURCES = \
    CholeskyFactorization.cpp \
    FactoredSolver.cpp \
	ForBESUtils.cpp \
	Function.cpp \
	IndBox.cpp \
	IndSOC.cpp \
	LDLFactorization.cpp \
	LinearOperator.cpp \
	Matrix.cpp \
	MatrixFactory.cpp \
	MatrixOperator.cpp \
	MatrixWriter.cpp \
	OpAdjoint.cpp \
	OpComposition.cpp \
	OpDCT2.cpp \
	OpDCT3.cpp \
	OpGradient.cpp \
	OpGradient2D.cpp \
	OpLinearCombination.cpp \
	OpReverseVector.cpp \
	OpSum.cpp \
	QuadOverAffine.cpp \
	QuadraticOperator.cpp \
	Quadratic.cpp

OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)

TESTS = \
	TestCholesky.test \
	TestIndBox.test \
	TestIndSOC.test \
	TestLDL.test \
	TestMatrix.test \
	TestMatrixFactory.test \
	TestMatrixOperator.test \
	TestOpAdjoint.test \
	TestOpComposition.test \
	TestOpDCT2.test \
	TestOpDCT3.test \
	TestOpGradient.test \
	TestOpReverseVector.test \
	TestQuadOverAffine.test \
	TestQuadratic.test \
	TestQuadraticOperator.test

TEST_BINS = $(TESTS:%.test=$(BIN_TEST_DIR)/%)
	
$(ARCHIVE): dirs $(OBJECTS)
	@echo "\nArchiving..."
	ar rcs $(ARCHIVE) $(OBJECTS)

all: $(ARCHIVE)
	
build-tests: $(ARCHIVE) $(TEST_BINS)

test: build-tests
	${BIN_TEST_DIR}/TestCholesky
	${BIN_TEST_DIR}/TestIndBox
	${BIN_TEST_DIR}/TestIndSOC
	${BIN_TEST_DIR}/TestLDL
	${BIN_TEST_DIR}/TestMatrix
	${BIN_TEST_DIR}/TestMatrixFactory
	${BIN_TEST_DIR}/TestMatrixOperator
	${BIN_TEST_DIR}/TestOpAdjoint
	${BIN_TEST_DIR}/TestOpComposition
	${BIN_TEST_DIR}/TestOpDCT2
	${BIN_TEST_DIR}/TestOpDCT3
	${BIN_TEST_DIR}/TestOpGradient
	${BIN_TEST_DIR}/TestOpReverseVector
	${BIN_TEST_DIR}/TestQuadOverAffine
	${BIN_TEST_DIR}/TestQuadratic
	${BIN_TEST_DIR}/TestQuadraticOperator
	
	

$(BIN_TEST_DIR)/%: $(OBJECTS) $(TEST_DIR)/%.cpp $(TEST_DIR)/%Runner.cpp $(TEST_DIR)/%.h
	@echo 
	@echo Compiling $*
	$(CXX) $(CFLAGS) $(IFLAGS) -o $(OBJ_TEST_DIR)/$*.o $(TEST_DIR)/$*.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) -o $(OBJ_TEST_DIR)/$*Runner.o $(TEST_DIR)/$*Runner.cpp
	@echo 
	@echo Linking $*
	$(CXX) $(LFLAGS) -o $(BIN_TEST_DIR)/$* $(OBJECTS) $(OBJ_TEST_DIR)/$*.o $(OBJ_TEST_DIR)/$*Runner.o $(lFLAGS) `cppunit-config --libs`

$(OBJ_DIR)/%.o: source/%.cpp
	@echo 
	@echo Compiling $*
	$(CXX) $(CFLAGS) $(IFLAGS) $< -o $@

dirs: $(OBJ_DIR) $(OBJ_TEST_DIR) $(BIN_DIR) $(BIN_TEST_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_TEST_DIR):
	mkdir -p $(OBJ_TEST_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_TEST_DIR):
	mkdir -p $(BIN_TEST_DIR)

clean:				
	@echo Cleaning...
	rm -rf $(OBJ_DIR)/*.o ;
	rm -rf $(OBJ_TEST_DIR)/*.o; 
	rm -rf $(BIN_DIR)/*;

help:
	@echo "Makefile targets for libforbes:\n"
	@echo "make                     - Compiles, links and archives [creates libforbes.a]"
	@echo "make clean               - Cleans all previously built files"
	@echo "make all                 - Same as make (tests are not built)"	
	@echo "make build-tests         - Compiles and links the tests"
	@echo "make test                - Compiles [if necessary] and runs all tests"
	@echo "make help                - This help message\n"
	
	
.SECONDARY: 

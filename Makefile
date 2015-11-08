# LIBFORBES MAKEFILE
	
# Licence note:
#
# ForBES is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#   
# ForBES is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with ForBES. If not, see <http://www.gnu.org/licenses/>.


#						 #
# Do not modify this file #
#						 #


include config.mk

	
# Enable parallel make on N-1 processors	
NPROCS := 1
OS:=$(shell uname -s)
ifeq ($(OS),Linux)
  NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
endif
ifeq ($(OS),Darwin) # Assume Mac OS X
  NPROCS:=$(shell sysctl -n hw.ncpu)
endif
NPROCS:=$$(($(NPROCS)-1))
#MAKEFLAGS += -j $(NPROCS)
#MAKEFLAGS += --no-print-directory

# C++ compiler
CXX = g++

# Enable CCACHE
ifneq (, $(shell which ccache))
	CXX := ccache $(CXX)
endif

# Additional compiler flags (e.g., -O2 or -O3 optimization flags, etc)
# To create a test coverage report add: -fprofile-arcs -ftest-coverage
CFLAGS_ADDITIONAL = -O0
#CFLAGS_ADDITIONAL += -fprofile-arcs
#CFLAGS_ADDITIONAL += -ftest-coverage

CFLAGS_WARNINGS = \
	-pedantic \
	-Wall \
	-Wextra \
	-Wcast-align \
	-Wcast-qual \
	-Wdisabled-optimization \
	-Wformat=2 \
	-Winit-self \
	-Wlogical-op \
	-Wmissing-declarations \
	-Wmissing-include-dirs \
	-Wnoexcept \
	-Wold-style-cast \
	-Woverloaded-virtual \
	-Wredundant-decls \
	-Wshadow \
	-Wsign-promo \
	-Wstrict-null-sentinel \
	-Wstrict-overflow=5 \
	-Wswitch-default \
	-Wundef \
	-Wno-unused \
	-Wno-sign-compare
	
CFLAGS_ADDITIONAL += ${CFLAGS_WARNINGS}

# Additional link flags
# To create a test coverage report add: -fprofile-arcs
LFLAGS_ADDITIONAL = -fprofile-arcs

OBJ_DIR = build/Debug
BIN_DIR = dist/Debug

OBJ_TEST_DIR = ${OBJ_DIR}/Test
BIN_TEST_DIR = ${BIN_DIR}/Test

TEST_DIR = source/tests

ARCHIVE = ${BIN_DIR}/libforbes.a

CFLAGS = -c -DUSE_LIBS $(CFLAGS_ADDITIONAL) $(CEXTRA)

IFLAGS = \
	-I. \
	-I./source \
	-I$(SS_DIR)/CHOLMOD/Include \
	-I$(SS_DIR)/LDL/Include \
	-I$(SS_DIR)/SuiteSparse_config \
	-I$(IEXTRA)


ifeq ($(OS),Linux)
	# Use the real time POSIX library on Linux
	LFLAGS_ADDITIONAL += -lrt
endif

lFLAGS = \
	-lldl \
	-lcholmod \
	-lamd \
	-lcolamd \
	-lsuitesparseconfig \
	-lccolamd \
	-lcamd \
	-llapacke \
	-llapack \
	-lopenblas \
	-lm \
	$(LFLAGS_ADDITIONAL)

LFLAGS = \
	-L$(SS_DIR)/CHOLMOD/Lib \
	-L$(SS_DIR)/AMD/Lib \
	-L$(SS_DIR)/COLAMD/Lib \
	-L$(SS_DIR)/SuiteSparse_config \
	-L$(SS_DIR)/CCOLAMD/Lib \
	-L$(SS_DIR)/CAMD/Lib \
	-L$(SS_DIR)/LDL/Lib \
	-L$(LEXTRA)

SOURCES = \
	LinSysSolver.cpp \
	S_LDLFactorization.cpp \
	CholeskyFactorization.cpp \
	FactoredSolver.cpp \
	ForBESUtils.cpp \
	Function.cpp \
	IndBox.cpp \
	IndSOC.cpp \
	IndPos.cpp \
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
	Quadratic.cpp \
	FunctionOntologicalClass.cpp \
	FunctionOntologyRegistry.cpp \
	DistanceToBox.cpp \
	ElasticNet.cpp \
	LogLogisticLoss.cpp \
	QuadraticLoss.cpp \
	HingeLoss.cpp \
	HuberLoss.cpp \
	Norm.cpp \
	Norm1.cpp \
	Norm2.cpp \
	SeparableSum.cpp \
	IndBall2.cpp \
	QuadraticLossOverAffine.cpp \
	SumOfNorm2.cpp \
	FBProblem.cpp \
	FBCache.cpp \
	ConjugateFunction.cpp

OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)

TESTS = \
	TestSLDL.test \
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
	TestQuadraticOperator.test \
	TestOntRegistry.test \
	TestDistanceToBox.test \
	TestElasticNet.test \
	TestQuadraticLoss.test \
	TestLogLogisticLoss.test \
	TestHingeLoss.test \
	TestHuber.test \
	TestNorm1.test \
	TestNorm2.test \
	TestFunctionOntologicalClass.test \
	TestFunctionOntologyRegistry.test \
	TestIndBall2.test \
	TestSeparableSum.test \
	TestFBCache.test \
	TestConjugateFunction.test

TEST_BINS = $(TESTS:%.test=$(BIN_TEST_DIR)/%)

$(ARCHIVE): dirs $(OBJECTS)
	@echo "\nArchiving..."
	ar rcs $(ARCHIVE) $(OBJECTS)

all: $(ARCHIVE)

build-tests: $(ARCHIVE) $(TEST_BINS)

test: build-tests
	@echo "\n*** FACTORIZERS ***"
	${BIN_TEST_DIR}/TestCholesky
	${BIN_TEST_DIR}/TestLDL
	${BIN_TEST_DIR}/TestSLDL
	@echo "\n*** FUNCTIONS ***"
	${BIN_TEST_DIR}/TestConjugateFunction
	${BIN_TEST_DIR}/TestQuadOverAffine
	${BIN_TEST_DIR}/TestQuadratic
	${BIN_TEST_DIR}/TestQuadraticOperator
	${BIN_TEST_DIR}/TestIndBox
	${BIN_TEST_DIR}/TestIndSOC
	${BIN_TEST_DIR}/TestNorm1
	${BIN_TEST_DIR}/TestNorm2
	${BIN_TEST_DIR}/TestIndBall2
	${BIN_TEST_DIR}/TestDistanceToBox
	${BIN_TEST_DIR}/TestElasticNet
	${BIN_TEST_DIR}/TestQuadraticLoss
	${BIN_TEST_DIR}/TestLogLogisticLoss
	${BIN_TEST_DIR}/TestHingeLoss
	${BIN_TEST_DIR}/TestHuber
	${BIN_TEST_DIR}/TestSeparableSum
	@echo "\n*** UTILITIES ***"
	${BIN_TEST_DIR}/TestMatrixFactory
	${BIN_TEST_DIR}/TestMatrix
	${BIN_TEST_DIR}/TestOntRegistry
	${BIN_TEST_DIR}/TestFunctionOntologicalClass
	${BIN_TEST_DIR}/TestFunctionOntologyRegistry
	@echo "\n*** LINEAR OPERATORS ***"
	${BIN_TEST_DIR}/TestMatrixOperator
	${BIN_TEST_DIR}/TestOpAdjoint
	${BIN_TEST_DIR}/TestOpComposition
	${BIN_TEST_DIR}/TestOpDCT2
	${BIN_TEST_DIR}/TestOpDCT3
	${BIN_TEST_DIR}/TestOpGradient
	${BIN_TEST_DIR}/TestOpReverseVector
	@echo "\n*** ALGORITHMS ***"
	${BIN_TEST_DIR}/TestFBCache

$(BIN_TEST_DIR)/%: $(OBJECTS) $(TEST_DIR)/%.cpp $(TEST_DIR)/%Runner.cpp $(TEST_DIR)/%.h
	@echo
	@echo [Compiling $*]
	$(CXX) $(CFLAGS) $(IFLAGS) -o $(OBJ_TEST_DIR)/$*.o $(TEST_DIR)/$*.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) -o $(OBJ_TEST_DIR)/$*Runner.o $(TEST_DIR)/$*Runner.cpp
	@echo
	@echo [Linking $*]
	$(CXX) $(LFLAGS) -o $(BIN_TEST_DIR)/$* $(OBJECTS) $(OBJ_TEST_DIR)/$*.o $(OBJ_TEST_DIR)/$*Runner.o $(lFLAGS) `cppunit-config --libs`
	@echo "\n\n\n"

main:
	$(CXX) $(CFLAGS) $(IFLAGS) source/main.cpp 
	$(CXX) $(LFLAGS) -L./dist/Debug main.o -lforbes $(lFLAGS)

$(OBJ_DIR)/%.o: source/%.cpp
	@echo 
	@echo Compiling $*
	$(CXX) $(CFLAGS) $(IFLAGS) $< -o $@
	@echo "\n\n\n"

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
	@echo "make					 - Compiles, links and archives [creates libforbes.a]"
	@echo "make clean			   - Cleans all previously built files"
	@echo "make all				 - Same as make (tests are not built)"
	@echo "make build-tests		 - Compiles and links the tests"
	@echo "make test				- Compiles [if necessary] and runs all tests"
	@echo "make docs				- Used doxygen to build documentation"
	@echo "make help				- This help message\n"

docs:
	doxygen forbes.doxygen

.SECONDARY:

.PHONY: clean

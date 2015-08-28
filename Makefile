include Configuration

CXX = g++

OBJ_DIR = build/Debug
BIN_DIR = dist/Debug

OBJ_TEST_DIR = ${OBJ_DIR}/Test
BIN_TEST_DIR = ${BIN_DIR}/Test

TEST_DIR = tests

CFLAGS = -c -DUSE_LIBS

IFLAGS = -I. -I./source -I$(SS_DIR)/CHOLMOD/Include -I$(SS_DIR)/LDL/Include -I$(SS_DIR)/SuiteSparse_config -I/usr/include/lapacke

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
	-lm -lrt

LFLAGS = -L$(SS_DIR)/CHOLMOD/Lib -L$(SS_DIR)/AMD/Lib -L$(SS_DIR)/COLAMD/Lib -L$(SS_DIR)/SuiteSparse_config -L$(SS_DIR)/CCOLAMD/Lib -L$(SS_DIR)/CAMD/Lib -L$(SS_DIR)/LDL/Lib -L/usr/lib64

SOURCES = \
	ForBESUtils.cpp \
	Function.cpp \
	Matrix.cpp \
	MatrixFactory.cpp \
	Quadratic.cpp \
	QuadOverAffine.cpp \
	LinearOperator.cpp \
	QuadraticOperator.cpp \
	FactoredSolver.cpp \
	CholeskyFactorization.cpp \
	LDLFactorization.cpp \
	IndBox.cpp \
	MatrixOperator.cpp

OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)

TESTS = \
	TestCholesky.test \
	TestMatrix.test \
	TestMatrixFactory.test \
	TestMatrixOperator.test \
	TestLDL.test \
	TestIndBox.test \
	TestQuadOverAffine.test \
	TestQuadratic.test \
	TestQuadraticOperator.test

build-tests: dirs $(TESTS)

%.test: $(OBJECTS) $(TEST_DIR)/%.cpp $(TEST_DIR)/%Runner.cpp $(TEST_DIR)/%.h
	@echo Compiling $*
	$(CXX) $(CFLAGS) $(IFLAGS) -o $(OBJ_TEST_DIR)/$*.o $(TEST_DIR)/$*.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) -o $(OBJ_TEST_DIR)/$*Runner.o $(TEST_DIR)/$*Runner.cpp
	@echo Linking $*
	$(CXX) $(LFLAGS) -o $(BIN_TEST_DIR)/$* $(OBJECTS) $(OBJ_TEST_DIR)/$*.o $(OBJ_TEST_DIR)/$*Runner.o $(lFLAGS) `cppunit-config --libs`

$(OBJ_DIR)/%.o: source/%.cpp
	@echo Compiling $*
	$(CXX) $(CFLAGS) $(IFLAGS) $< -o $@
			

dirs:						
	@echo Creating directories
	mkdir -p $(OBJ_DIR); 
	mkdir -p $(OBJ_TEST_DIR); 
	mkdir -p $(BIN_DIR); 
	mkdir -p $(BIN_TEST_DIR)

clean:				
	@echo Cleaning
	rm -rf $(OBJ_DIR)/*.o ;
	rm -rf $(OBJ_TEST_DIR)/*.o; 
	rm -rf $(BIN_DIR)/*;

test:
	${BIN_TEST_DIR}/TestMatrix
	${BIN_TEST_DIR}/TestLDL
	${BIN_TEST_DIR}/TestCholesky
	

.SECONDARY: 

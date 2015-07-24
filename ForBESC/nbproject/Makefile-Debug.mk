#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/529256040/ForBESUtils.o \
	${OBJECTDIR}/_ext/529256040/QuadOverAffine.o \
	${OBJECTDIR}/Function.o \
	${OBJECTDIR}/Matrix.o \
	${OBJECTDIR}/MatrixFactory.o \
	${OBJECTDIR}/Quadratic.o \
	${OBJECTDIR}/main.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f3 \
	${TESTDIR}/TestFiles/f2

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/home/chung/Documents/MATLAB/SuiteSparse/CHOLMOD/Lib -L/home/chung/Documents/MATLAB/SuiteSparse/AMD/Lib -L/home/chung/Documents/MATLAB/SuiteSparse/COLAMD/Lib -L/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -L/home/chung/Documents/MATLAB/SuiteSparse/CCOLAMD/Lib -L/home/chung/Documents/MATLAB/SuiteSparse/CAMD/Lib -Wl,-rpath,. -lcholmod -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd -llapacke -lblas -llapack -lopenblas

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/forbesc

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/forbesc: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/forbesc ${OBJECTFILES} ${LDLIBSOPTIONS} -lm -lrt

${OBJECTDIR}/_ext/529256040/ForBESUtils.o: /home/chung/Documents/MATLAB/ForBES/ForBESC/ForBESUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/529256040
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/529256040/ForBESUtils.o /home/chung/Documents/MATLAB/ForBES/ForBESC/ForBESUtils.cpp

${OBJECTDIR}/_ext/529256040/QuadOverAffine.o: /home/chung/Documents/MATLAB/ForBES/ForBESC/QuadOverAffine.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/529256040
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/529256040/QuadOverAffine.o /home/chung/Documents/MATLAB/ForBES/ForBESC/QuadOverAffine.cpp

${OBJECTDIR}/Function.o: Function.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Function.o Function.cpp

${OBJECTDIR}/Matrix.o: Matrix.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Matrix.o Matrix.cpp

${OBJECTDIR}/MatrixFactory.o: MatrixFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixFactory.o MatrixFactory.cpp

${OBJECTDIR}/Quadratic.o: Quadratic.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Quadratic.o Quadratic.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/TestMatrix.o ${TESTDIR}/tests/TestMatrixRunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} `cppunit-config --libs` `cppunit-config --libs`   

${TESTDIR}/TestFiles/f3: ${TESTDIR}/tests/TestMatrixFactory.o ${TESTDIR}/tests/TestMatrixFactoryRunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f3 $^ ${LDLIBSOPTIONS} `cppunit-config --libs`   

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/TestQuadratic.o ${TESTDIR}/tests/TestQuadraticRunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} `cppunit-config --libs`   


${TESTDIR}/tests/TestMatrix.o: tests/TestMatrix.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -I. -include Matrix.h -include Matrix.h `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TestMatrix.o tests/TestMatrix.cpp


${TESTDIR}/tests/TestMatrixRunner.o: tests/TestMatrixRunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -I. -include Matrix.h -include Matrix.h `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TestMatrixRunner.o tests/TestMatrixRunner.cpp


${TESTDIR}/tests/TestMatrixFactory.o: tests/TestMatrixFactory.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -I. -include Matrix.h `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TestMatrixFactory.o tests/TestMatrixFactory.cpp


${TESTDIR}/tests/TestMatrixFactoryRunner.o: tests/TestMatrixFactoryRunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -I. -include Matrix.h `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TestMatrixFactoryRunner.o tests/TestMatrixFactoryRunner.cpp


${TESTDIR}/tests/TestQuadratic.o: tests/TestQuadratic.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -I. -include Matrix.h `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TestQuadratic.o tests/TestQuadratic.cpp


${TESTDIR}/tests/TestQuadraticRunner.o: tests/TestQuadraticRunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -I. -include Matrix.h `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TestQuadraticRunner.o tests/TestQuadraticRunner.cpp


${OBJECTDIR}/_ext/529256040/ForBESUtils_nomain.o: ${OBJECTDIR}/_ext/529256040/ForBESUtils.o /home/chung/Documents/MATLAB/ForBES/ForBESC/ForBESUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/529256040
	@NMOUTPUT=`${NM} ${OBJECTDIR}/_ext/529256040/ForBESUtils.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/529256040/ForBESUtils_nomain.o /home/chung/Documents/MATLAB/ForBES/ForBESC/ForBESUtils.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/_ext/529256040/ForBESUtils.o ${OBJECTDIR}/_ext/529256040/ForBESUtils_nomain.o;\
	fi

${OBJECTDIR}/_ext/529256040/QuadOverAffine_nomain.o: ${OBJECTDIR}/_ext/529256040/QuadOverAffine.o /home/chung/Documents/MATLAB/ForBES/ForBESC/QuadOverAffine.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/529256040
	@NMOUTPUT=`${NM} ${OBJECTDIR}/_ext/529256040/QuadOverAffine.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/529256040/QuadOverAffine_nomain.o /home/chung/Documents/MATLAB/ForBES/ForBESC/QuadOverAffine.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/_ext/529256040/QuadOverAffine.o ${OBJECTDIR}/_ext/529256040/QuadOverAffine_nomain.o;\
	fi

${OBJECTDIR}/Function_nomain.o: ${OBJECTDIR}/Function.o Function.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/Function.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Function_nomain.o Function.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/Function.o ${OBJECTDIR}/Function_nomain.o;\
	fi

${OBJECTDIR}/Matrix_nomain.o: ${OBJECTDIR}/Matrix.o Matrix.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/Matrix.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Matrix_nomain.o Matrix.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/Matrix.o ${OBJECTDIR}/Matrix_nomain.o;\
	fi

${OBJECTDIR}/MatrixFactory_nomain.o: ${OBJECTDIR}/MatrixFactory.o MatrixFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/MatrixFactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatrixFactory_nomain.o MatrixFactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/MatrixFactory.o ${OBJECTDIR}/MatrixFactory_nomain.o;\
	fi

${OBJECTDIR}/Quadratic_nomain.o: ${OBJECTDIR}/Quadratic.o Quadratic.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/Quadratic.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Quadratic_nomain.o Quadratic.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/Quadratic.o ${OBJECTDIR}/Quadratic_nomain.o;\
	fi

${OBJECTDIR}/main_nomain.o: ${OBJECTDIR}/main.o main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/main.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DUSE_LIBS -I../../SuiteSparse/CHOLMOD/Include -I/home/chung/Documents/MATLAB/SuiteSparse/SuiteSparse_config -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_nomain.o main.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/main.o ${OBJECTDIR}/main_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f3 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/forbesc

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

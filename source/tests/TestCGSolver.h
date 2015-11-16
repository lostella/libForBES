/*
 * File:   TestCGSolver.h
 * Author: chung
 *
 * Created on Nov 16, 2015, 6:41:56 PM
 */

#ifndef TESTCGSOLVER_H
#define	TESTCGSOLVER_H

#define FORBES_TEST_UTILS

#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestCGSolver : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestCGSolver);

    CPPUNIT_TEST(testSolve);
    CPPUNIT_TEST(testSolve2);

    CPPUNIT_TEST_SUITE_END();

public:
    TestCGSolver();
    virtual ~TestCGSolver();
    void setUp();
    void tearDown();

private:
    void testSolve();
    void testSolve2();

};

#endif	/* TESTCGSOLVER_H */


/*
 * File:   TestLDL.h
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 4, 2015, 8:18:33 PM
 */

#ifndef TESTLDL_H
#define	TESTLDL_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "Matrix.h"
#include "MatrixFactory.h"
#include "LDLFactorization.h"

#include <cppunit/extensions/HelperMacros.h>

class TestLDL : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestLDL);

    CPPUNIT_TEST(testSolveDense);
    CPPUNIT_TEST(testSolveSymmetric);
    CPPUNIT_TEST(testSolveSparse);
    CPPUNIT_TEST(testSolveSparse2);

    CPPUNIT_TEST_SUITE_END();

public:
    TestLDL();
    virtual ~TestLDL();
    void setUp();
    void tearDown();

private:
    void testSolveDense();
    void testSolveSymmetric();
    void testSolveSparse();
    void testSolveSparse2();

};

#endif	/* TESTLDL_H */


/*
 * File:   TestMatrixOperator.h
 * Author: Pantelis Sopasakis
 *
 * Created on Jul 24, 2015, 8:39:54 PM
 */

#ifndef TESTMATRIXOPERATOR_H
#define	TESTMATRIXOPERATOR_H

#include <cppunit/extensions/HelperMacros.h>
#include "MatrixOperator.h"
#include "MatrixFactory.h"

#define FORBES_TEST_UTILS
#include "ForBES.h"

class TestMatrixOperator : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestMatrixOperator);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallAdjoint);

    CPPUNIT_TEST_SUITE_END();

public:
    TestMatrixOperator();
    virtual ~TestMatrixOperator();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallAdjoint();
    
};

#endif	/* TESTMATRIXOPERATOR_H */

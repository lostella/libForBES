/*
 * File:   TestQuadraticOperator.h
 * Author: Pantelis Sopasakis
 *
 * Created on Jul 24, 2015, 8:52:27 PM
 */

#ifndef TESTQUADRATICOPERATOR_H
#define	TESTQUADRATICOPERATOR_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS

#include "Matrix.h"
#include "MatrixFactory.h"
#include "LinearOperator.h"
#include "MatrixOperator.h"
#include "QuadraticOperator.h"
#include "ForBES.h"

class TestQuadraticOperator : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadraticOperator);

    CPPUNIT_TEST(testCall);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadraticOperator();
    virtual ~TestQuadraticOperator();
    void setUp();
    void tearDown();

private:
    void testCall();

};

#endif	/* TESTQUADRATICOPERATOR_H */

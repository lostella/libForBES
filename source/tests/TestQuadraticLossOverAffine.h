/*
 * File:   TestQuadraticLossOverAffine.h
 * Author: chung
 *
 * Created on Nov 5, 2015, 12:06:30 AM
 */

#ifndef TESTQUADRATICLOSSOVERAFFINE_H
#define	TESTQUADRATICLOSSOVERAFFINE_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestQuadraticLossOverAffine : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadraticLossOverAffine);

    CPPUNIT_TEST(testMethod);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadraticLossOverAffine();
    virtual ~TestQuadraticLossOverAffine();
    void setUp();
    void tearDown();

private:
    void testMethod();
};

#endif	/* TESTQUADRATICLOSSOVERAFFINE_H */


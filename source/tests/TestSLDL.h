/*
 * File:   TestSLDL.h
 * Author: chung
 *
 * Created on Nov 5, 2015, 4:08:41 AM
 */

#ifndef TESTSLDL_H
#define	TESTSLDL_H

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include <cppunit/extensions/HelperMacros.h>

class TestSLDL : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestSLDL);

    CPPUNIT_TEST(testFactorizeAndSolve);
    CPPUNIT_TEST(testDenseShort);
    CPPUNIT_TEST(testDenseTall);

    CPPUNIT_TEST_SUITE_END();

public:
    TestSLDL();
    virtual ~TestSLDL();
    void setUp();
    void tearDown();

private:
    void testFactorizeAndSolve();
    void testDenseShort();
    void testDenseTall();
};

#endif	/* TESTSLDL_H */


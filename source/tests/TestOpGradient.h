/*
 * File:   TestOpGradient.h
 * Author: chung
 *
 * Created on Sep 16, 2015, 1:57:49 AM
 */

#ifndef TESTOPGRADIENT_H
#define	TESTOPGRADIENT_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpGradient.h"

#include <cppunit/extensions/HelperMacros.h>

class TestOpGradient : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpGradient);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testLinearity);
    CPPUNIT_TEST(testAdjointLinearity);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpGradient();
    virtual ~TestOpGradient();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testLinearity();
    void testAdjointLinearity();

};

#endif	/* TESTOPGRADIENT_H */


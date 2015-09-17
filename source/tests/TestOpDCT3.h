/*
 * File:   TestOpDCT3.h
 * Author: chung
 *
 * Created on Sep 16, 2015, 3:42:02 AM
 */

#ifndef TESTOPDCT3_H
#define	TESTOPDCT3_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpDCT3.h"

#include <cppunit/extensions/HelperMacros.h>

class TestOpDCT3 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpDCT3);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testLinearity);
    CPPUNIT_TEST(testAdjointLinearity);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpDCT3();
    virtual ~TestOpDCT3();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testLinearity();
    void testAdjointLinearity();

};

#endif	/* TESTOPDCT3_H */


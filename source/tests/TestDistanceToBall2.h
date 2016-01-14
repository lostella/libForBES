/*
 * File:   TestDistanceToBall2.h
 * Author: chung
 *
 * Created on Jan 14, 2016, 4:32:12 PM
 */

#ifndef TESTDISTANCETOBALL2_H
#define	TESTDISTANCETOBALL2_H

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include "DistanceToBall2.h"
#include <math.h>

#include <cppunit/extensions/HelperMacros.h>

class TestDistanceToBall2 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestDistanceToBall2);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCall2);
    CPPUNIT_TEST(testCall3);
    CPPUNIT_TEST(testCall4);
    CPPUNIT_TEST(testCall5);
    CPPUNIT_TEST(testGradient);
    CPPUNIT_TEST(testGradient2);

    CPPUNIT_TEST_SUITE_END();

public:
    TestDistanceToBall2();
    virtual ~TestDistanceToBall2();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCall2();
    void testCall3();
    void testCall4();
    void testCall5();
    void testGradient();
    void testGradient2();

};

#endif	/* TESTDISTANCETOBALL2_H */


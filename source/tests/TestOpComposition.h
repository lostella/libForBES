/*
 * File:   TestOpComposition.h
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 2:40:14 AM
 */

#ifndef TESTOPCOMPOSITION_H
#define	TESTOPCOMPOSITION_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpComposition.h"
#include <cppunit/extensions/HelperMacros.h>

class TestOpComposition : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpComposition);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallAdjoint);
    CPPUNIT_TEST(testDimension);    

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpComposition();
    virtual ~TestOpComposition();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallAdjoint();
    void testDimension();

};

#endif	/* TESTOPCOMPOSITION_H */


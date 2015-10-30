/*
 * File:   TestNorm1.h
 * Author: chung
 *
 * Created on Oct 30, 2015, 6:15:56 PM
 */

#ifndef TESTNORM1_H
#define	TESTNORM1_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestNorm1 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestNorm1);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testCallProx2);
    CPPUNIT_TEST(testDualNorm);

    CPPUNIT_TEST_SUITE_END();

public:
    TestNorm1();
    virtual ~TestNorm1();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();
    void testCallProx2();
    void testDualNorm();

};

#endif	/* TESTNORM1_H */


/*
 * File:   TestIndBall2.h
 * Author: chung
 *
 * Created on Nov 3, 2015, 5:21:20 PM
 */

#ifndef TESTINDBALL2_H
#define	TESTINDBALL2_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestIndBall2 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestIndBall2);

    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testCategory);

    CPPUNIT_TEST_SUITE_END();

public:
    TestIndBall2();
    virtual ~TestIndBall2();
    void setUp();
    void tearDown();

private:
    void testCallProx();
    void testCategory();

};

#endif	/* TESTINDBALL2_H */


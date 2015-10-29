/*
 * File:   TestHingeLoss.h
 * Author: chung
 *
 * Created on Oct 29, 2015, 11:36:49 PM
 */

#ifndef TESTHINGELOSS_H
#define	TESTHINGELOSS_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestHingeLoss : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestHingeLoss);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);

    CPPUNIT_TEST_SUITE_END();

public:
    TestHingeLoss();
    virtual ~TestHingeLoss();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();

};

#endif	/* TESTHINGELOSS_H */


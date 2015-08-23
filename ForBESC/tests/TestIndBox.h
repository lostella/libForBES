/*
 * File:   TestIndBox.h
 * Author: chung
 *
 * Created on Jul 26, 2015, 5:41:14 PM
 */

#ifndef TESTINDBOX_H
#define	TESTINDBOX_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include "IndBox.h"

class TestIndBox : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestIndBox);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCall2);
    CPPUNIT_TEST(testCall3);
    CPPUNIT_TEST(testCallProx);

    CPPUNIT_TEST_SUITE_END();

public:
    TestIndBox();
    virtual ~TestIndBox();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCall2();
    void testCall3();
    void testCallProx();

};

#endif	/* TESTINDBOX_H */


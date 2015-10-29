/*
 * File:   TestElasticNet.h
 * Author: chung
 *
 * Created on Oct 29, 2015, 6:51:55 PM
 */

#ifndef TESTELASTICNET_H
#define	TESTELASTICNET_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestElasticNet : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestElasticNet);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testOther);

    CPPUNIT_TEST_SUITE_END();

public:
    TestElasticNet();
    virtual ~TestElasticNet();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();
    void testOther();

};

#endif	/* TESTELASTICNET_H */


/*
 * File:   TestLogLogisticLoss.h
 * Author: chung
 *
 * Created on Oct 29, 2015, 7:38:05 PM
 */

#ifndef TESTLOGLOGISTICLOSS_H
#define	TESTLOGLOGISTICLOSS_H

#define FORBES_TEST_UTILS
#include "ForBES.h"


#include <cppunit/extensions/HelperMacros.h>

class TestLogLogisticLoss : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestLogLogisticLoss);

    CPPUNIT_TEST(testCall);

    CPPUNIT_TEST_SUITE_END();

public:
    TestLogLogisticLoss();
    virtual ~TestLogLogisticLoss();
    void setUp();
    void tearDown();

private:
    void testCall();

};

#endif	/* TESTLOGLOGISTICLOSS_H */


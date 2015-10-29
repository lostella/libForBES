/*
 * File:   TestQuadraticLoss.h
 * Author: chung
 *
 * Created on Oct 29, 2015, 7:32:54 PM
 */

#ifndef TESTQUADRATICLOSS_H
#define	TESTQUADRATICLOSS_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestQuadraticLoss : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadraticLoss);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallConj);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadraticLoss();
    virtual ~TestQuadraticLoss();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallConj();

};

#endif	/* TESTQUADRATICLOSS_H */


/*
 * File:   TestOpReverseVector.h
 * Author: chung
 *
 * Created on Sep 16, 2015, 4:00:46 AM
 */

#ifndef TESTOPREVERSEVECTOR_H
#define	TESTOPREVERSEVECTOR_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpReverseVector.h"
#include "MatrixFactory.h"

#include <cppunit/extensions/HelperMacros.h>

class TestOpReverseVector : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpReverseVector);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallNotFixedSize);
    CPPUNIT_TEST(testCallAdjoint);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpReverseVector();
    virtual ~TestOpReverseVector();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallNotFixedSize();
    void testCallAdjoint();

};

#endif	/* TESTOPREVERSEVECTOR_H */


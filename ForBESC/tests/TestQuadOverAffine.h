/*
 * File:   TestQuadOverAffine.h
 * Author: chung
 *
 * Created on Aug 27, 2015, 2:19:01 PM
 */

#ifndef TESTQUADOVERAFFINE_H
#define	TESTQUADOVERAFFINE_H

#include <cppunit/extensions/HelperMacros.h>

class TestQuadOverAffine : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadOverAffine);

    CPPUNIT_TEST(testMethod);
    CPPUNIT_TEST(testFailedMethod);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadOverAffine();
    virtual ~TestQuadOverAffine();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();
};

#endif	/* TESTQUADOVERAFFINE_H */


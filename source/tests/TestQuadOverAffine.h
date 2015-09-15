/*
 * File:   TestQuadOverAffine.h
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 27, 2015, 2:19:01 PM
 */

#ifndef TESTQUADOVERAFFINE_H
#define	TESTQUADOVERAFFINE_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"
#include <cppunit/extensions/HelperMacros.h>

class TestQuadOverAffine : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadOverAffine);

    CPPUNIT_TEST(testQuadOverAffine);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadOverAffine();
    virtual ~TestQuadOverAffine();
    void setUp();
    void tearDown();

private:
    void testQuadOverAffine();
};

#endif	/* TESTQUADOVERAFFINE_H */

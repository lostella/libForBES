/*
 * File:   TestCholesky.h
 * Author: chung
 *
 * Created on Aug 5, 2015, 12:28:54 AM
 */

#ifndef TESTCHOLESKY_H
#define	TESTCHOLESKY_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"

#include <cppunit/extensions/HelperMacros.h>

class TestCholesky : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestCholesky);

    CPPUNIT_TEST(testCholeskyDense);
    CPPUNIT_TEST(testCholeskySymmetric);
    CPPUNIT_TEST(testCholeskySymmetric2);
    CPPUNIT_TEST(testCholeskySparse);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestCholesky();
    virtual ~TestCholesky();
    void setUp();
    void tearDown();

private:
    void testCholeskyDense();
    void testCholeskySymmetric();
    void testCholeskySymmetric2();
    void testCholeskySparse();
    
};

#endif	/* TESTCHOLESKY_H */


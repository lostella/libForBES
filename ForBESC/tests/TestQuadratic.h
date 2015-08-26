/*
 * File:   TestQuadratic.h
 * Author: Pantelis Sopasakis
 *
 * Created on Jul 9, 2015, 4:14:39 AM
 */

#ifndef TESTQUADRATIC_H
#define	TESTQUADRATIC_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include <cmath>

class TestQuadratic : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadratic);

    CPPUNIT_TEST(testQuadratic);
    CPPUNIT_TEST(testQuadratic2);
    CPPUNIT_TEST(testQuadratic3);
    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallWithGradient);
    CPPUNIT_TEST(testCallConj);
    CPPUNIT_TEST(testCallConj2);
    CPPUNIT_TEST(testCategory);
    CPPUNIT_TEST(testCallDiagonalMatrix);
    CPPUNIT_TEST(testCallSparse);
    CPPUNIT_TEST(testCallSparse2);
    CPPUNIT_TEST(testCallSparse3);
    CPPUNIT_TEST(testCallConjSparse);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadratic();
    virtual ~TestQuadratic();
    void setUp();
    void tearDown();

private:
    void testQuadratic();
    void testQuadratic2();
    void testQuadratic3();
    void testCallProx();
    void testCall();
    void testCallWithGradient();
    void testCallConj();
    void testCallConj2();
    void testCategory();
    
    void testCallSparse();
    void testCallSparse2();
    void testCallSparse3();
    void testCallConjSparse();
    
    void testCallDiagonalMatrix();

};

#endif	/* TESTQUADRATIC_H */


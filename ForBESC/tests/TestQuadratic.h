/*
 * File:   TestQuadratic.h
 * Author: chung
 *
 * Created on Jul 9, 2015, 4:14:39 AM
 */

#ifndef TESTQUADRATIC_H
#define	TESTQUADRATIC_H

#include <cppunit/extensions/HelperMacros.h>

class TestQuadratic : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadratic);

    CPPUNIT_TEST(testQuadratic);
    CPPUNIT_TEST(testQuadratic2);
    CPPUNIT_TEST(testQuadratic3);
    CPPUNIT_TEST(testQuadratic4);
    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallWithGradient);
    CPPUNIT_TEST(testCallConj);
    CPPUNIT_TEST(testCategory);
    CPPUNIT_TEST(testCallDiagonalMatrix);

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
    void testQuadratic4();
    void testCall();
    void testCallWithGradient();
    void testCallConj();
    void testCategory();
    
    void testCallDiagonalMatrix();

};

#endif	/* TESTQUADRATIC_H */


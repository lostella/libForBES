/*
 * File:   TestMatrixExtras.h
 * Author: chung
 *
 * Created on Nov 8, 2015, 4:33:45 PM
 */

#ifndef TESTMATRIXEXTRAS_H
#define	TESTMATRIXEXTRAS_H

#define FORBES_TEST_UTILS

#include <cppunit/extensions/HelperMacros.h>

class TestMatrixExtras : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestMatrixExtras);

    CPPUNIT_TEST(test_add_DD);
    CPPUNIT_TEST(test_add_DS);
    CPPUNIT_TEST(test_add_SS);
    
    CPPUNIT_TEST(test_mult_DD);
    CPPUNIT_TEST(test_mult_DDT);
    CPPUNIT_TEST(test_mult_DX);
    CPPUNIT_TEST(test_mult_DH);
    
    CPPUNIT_TEST(test_mult_SS);
    CPPUNIT_TEST(test_mult_SS2);
    CPPUNIT_TEST(test_mult_SS3);

    CPPUNIT_TEST_SUITE_END();

public:
    TestMatrixExtras();
    virtual ~TestMatrixExtras();
    void setUp();
    void tearDown();

private:
    void test_add_DD(); // DENSE  := gamma * DENSE  + alpha * DENSE
    void test_add_DS(); // DENSE  := gamma * DENSE  + alpha * SPARSE
    void test_add_SS(); // SPARSE := gamma * SPARSE + alpha * SPARSE
    
    void test_mult_DD(); // DENSE := gamma * DENSE + alpha * DENSE * DENSE
    void test_mult_DDT(); // DENSE := gamma * DENSE + alpha * DENSE * DENSE'
    void test_mult_DX(); // DENSE := gamma * DENSE + alpha * DENSE * DIAGONAL
    void test_mult_DH(); // DENSE := gamma * DENSE + alpha * DENSE * SYMMETRIC
    void test_mult_SS(); // SPARSE := gamma * SPARSE + alpha * SPARSE * SPARSE
    void test_mult_SS2(); // SPARSE := gamma * SPARSE + alpha * SPARSE * SPARSE
    void test_mult_SS3(); // SPARSE := gamma * SPARSE + alpha * SPARSE * SPARSE
};

#endif	/* TESTMATRIXEXTRAS_H */


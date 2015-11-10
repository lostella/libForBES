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
    CPPUNIT_TEST(test_add_DDT);
    CPPUNIT_TEST(test_add_DTD);
    CPPUNIT_TEST(test_add_DTDT);
    
    CPPUNIT_TEST(test_add_DS);
    
    CPPUNIT_TEST(test_add_SS);
    CPPUNIT_TEST(test_add_SST);
    
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
    void test_add_DD(); 
    void test_add_DDT();
    void test_add_DTD();
    void test_add_DTDT();
    
    void test_add_DS();
    void test_add_SS();
    void test_add_SST();
    
    
    
    void test_mult_DD();
    void test_mult_DDT();
    void test_mult_DX();
    void test_mult_DH();
    void test_mult_SS();
    void test_mult_SS2();
    void test_mult_SS3();
};

#endif	/* TESTMATRIXEXTRAS_H */


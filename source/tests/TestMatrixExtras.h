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
    CPPUNIT_TEST(test_add_SS);

    CPPUNIT_TEST_SUITE_END();

public:
    TestMatrixExtras();
    virtual ~TestMatrixExtras();
    void setUp();
    void tearDown();

private:
    void test_add_DD(); // DENSE  := gamma * DENSE  + alpha * DENSE
    void test_add_SS(); // SPARSE := gamma * SPARSE + alpha * SPARSE
};

#endif	/* TESTMATRIXEXTRAS_H */


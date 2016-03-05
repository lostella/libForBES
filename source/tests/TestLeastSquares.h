/*
 * File:   TestLeastSquares.h
 * Author: chung
 *
 * Created on Mar 3, 2016, 2:47:16 PM
 */

#ifndef TESTLEASTSQUARES_H
#define	TESTLEASTSQUARES_H

#include <cppunit/extensions/HelperMacros.h>

class TestLeastSquares : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestLeastSquares);

    CPPUNIT_TEST(testSolveLS);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestLeastSquares();
    virtual ~TestLeastSquares();
    void setUp();
    void tearDown();

private:
    void testSolveLS();
    
};

#endif	/* TESTLEASTSQUARES_H */


/*
 * File:   TestOpDCT2.h
 * Author: chung
 *
 * Created on Sep 15, 2015, 4:24:44 PM
 */

#ifndef TESTOPDCT2_H
#define	TESTOPDCT2_H

#include <cppunit/extensions/HelperMacros.h>

class TestOpDCT2 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpDCT2);

    CPPUNIT_TEST(testCall);   

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpDCT2();
    virtual ~TestOpDCT2();
    void setUp();
    void tearDown();

private:
    void testCall();
    
};

#endif	/* TESTOPDCT2_H */


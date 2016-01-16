#ifndef TESTFBSPLITTING_H
#define	TESTFBSPLITTING_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"

#include <cppunit/extensions/HelperMacros.h>

class TestFBSplitting : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestFBSplitting);

    CPPUNIT_TEST(testBoxQP_small);
    CPPUNIT_TEST(testLasso_small);
    CPPUNIT_TEST(testSparseLogReg_small);
    
    CPPUNIT_TEST_SUITE_END();

public:
    TestFBSplitting();
    virtual ~TestFBSplitting();
    void setUp();
    void tearDown();

private:
    void testBoxQP_small();
    void testLasso_small();
    void testSparseLogReg_small();
};

#endif	/* TESTFBSPLITTING_H */

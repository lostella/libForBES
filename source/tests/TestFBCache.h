#ifndef TESTFBCACHE_H
#define	TESTFBCACHE_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"

#include <cppunit/extensions/HelperMacros.h>

class TestFBCache : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestFBCache);

    CPPUNIT_TEST(testBoxQP_small);
    CPPUNIT_TEST(testLasso_small);
    CPPUNIT_TEST(testSparseLogReg_small);
    
    CPPUNIT_TEST_SUITE_END();

public:
    TestFBCache();
    virtual ~TestFBCache();
    void setUp();
    void tearDown();

private:
    void testBoxQP_small();
    void testLasso_small();
    void testSparseLogReg_small();
};

#endif	/* TESTFBCACHE_H */


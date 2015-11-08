#ifndef TESTFBCACHE_H
#define	TESTFBCACHE_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"

#include <cppunit/extensions/HelperMacros.h>

class TestFBCache : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestFBCache);

    CPPUNIT_TEST(test_boxqp_small);
    
    CPPUNIT_TEST_SUITE_END();

public:
    TestFBCache();
    virtual ~TestFBCache();
    void setUp();
    void tearDown();

private:
    void test_boxqp_small();
    void test_l1logreg_small();
};

#endif	/* TESTFBCACHE_H */

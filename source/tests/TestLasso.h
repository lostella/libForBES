#ifndef TESTLASSO_H
#define	TESTLASSO_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"

#include <cppunit/extensions/HelperMacros.h>

class TestLasso : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestLasso);

    CPPUNIT_TEST(runTest);
    
    CPPUNIT_TEST_SUITE_END();

public:
    TestLasso();
    virtual ~TestLasso();
    void setUp();
    void tearDown();

private:
    void runTest();
};

#endif	/* TESTLASSO_H */


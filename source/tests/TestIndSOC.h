/*
 * File:   TestIndSOC.h
 * Author: Lorenzo Stella
 *
 * Created on Sept 21, 2015, 11:59 AM
 */

#ifndef TESTINDSOC_H
#define	TESTINDSOC_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include "IndSOC.h"

class TestIndSOC : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestIndSOC);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);

    CPPUNIT_TEST_SUITE_END();

public:
    TestIndSOC();
    virtual ~TestIndSOC();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();

};

#endif	/* TESTINDSOC_H */

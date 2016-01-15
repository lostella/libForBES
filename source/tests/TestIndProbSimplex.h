/*
 * File:   TestIndProbSimplex.h
 * Author: chung
 *
 * Created on Jan 15, 2016, 5:38:56 PM
 */

#ifndef TESTINDPROBSIMPLEX_H
#define	TESTINDPROBSIMPLEX_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include <cppunit/extensions/HelperMacros.h>

class TestIndProbSimplex : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestIndProbSimplex);

    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testCallProxLarge);

    CPPUNIT_TEST_SUITE_END();

public:
    TestIndProbSimplex();
    virtual ~TestIndProbSimplex();
    void setUp();
    void tearDown();

private:
    void testCallProx();
    void testCallProxLarge();

};

#endif	/* TESTINDPROBSIMPLEX_H */


/*
 * File:   TestConjugateFunction.h
 * Author: chung
 *
 * Created on Nov 7, 2015, 4:10:56 PM
 */

#ifndef TESTCONJUGATEFUNCTION_H
#define	TESTCONJUGATEFUNCTION_H

#include <cppunit/extensions/HelperMacros.h>

class TestConjugateFunction : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestConjugateFunction);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCall2);
    CPPUNIT_TEST(testCallConj);
    CPPUNIT_TEST(testCallConj2);
    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testCategory);

    CPPUNIT_TEST_SUITE_END();

public:
    TestConjugateFunction();
    virtual ~TestConjugateFunction();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCall2();
    void testCallConj();
    void testCallConj2();
    void testCallProx();
    void testCategory();

};

#endif	/* TESTCONJUGATEFUNCTION_H */


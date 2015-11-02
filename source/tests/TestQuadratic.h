/*
 * File:   TestQuadratic.h
 * Author: Pantelis Sopasakis
 *
 * Created on Jul 9, 2015, 4:14:39 AM
 * 
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TESTQUADRATIC_H
#define	TESTQUADRATIC_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include <cmath>

class TestQuadratic : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestQuadratic);

    CPPUNIT_TEST(testQuadratic);
    CPPUNIT_TEST(testQuadratic2);
    CPPUNIT_TEST(testQuadratic3);
    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallWithGradient);
    CPPUNIT_TEST(testCallConj);
    CPPUNIT_TEST(testCallConj2);
    CPPUNIT_TEST(testCallDiagonalMatrix);
    CPPUNIT_TEST(testCallSparse);
    CPPUNIT_TEST(testCallSparse2);
    CPPUNIT_TEST(testCallSparse3);
    CPPUNIT_TEST(testCallConjSparse);

    CPPUNIT_TEST_SUITE_END();

public:
    TestQuadratic();
    virtual ~TestQuadratic();
    void setUp();
    void tearDown();

private:
    void testQuadratic();
    void testQuadratic2();
    void testQuadratic3();
    void testCallProx();
    void testCall();
    void testCallWithGradient();
    void testCallConj();
    void testCallConj2();
    
    void testCallSparse();
    void testCallSparse2();
    void testCallSparse3();
    void testCallConjSparse();
    
    void testCallDiagonalMatrix();

};

#endif	/* TESTQUADRATIC_H */

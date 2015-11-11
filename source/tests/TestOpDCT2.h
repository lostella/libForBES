/*
 * File:   TestOpDCT2.h
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 4:24:44 PM
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

#ifndef TESTOPDCT2_H
#define	TESTOPDCT2_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "LinearOperator.h"
#include "OpDCT2.h"
#include "MatrixFactory.h"
#include <cppunit/extensions/HelperMacros.h>

class TestOpDCT2 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpDCT2);

    CPPUNIT_TEST(testCall0);   
    CPPUNIT_TEST(testCall1);   
    CPPUNIT_TEST(testAdj0);   
    CPPUNIT_TEST(testCall);   
    CPPUNIT_TEST(testLinearity);   
    CPPUNIT_TEST(testAdjointLinearity);   

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpDCT2();
    virtual ~TestOpDCT2();
    void setUp();
    void tearDown();

private:
    void testCall0();
    void testCall1();
    void testAdj0();
    void testCall();
    void testLinearity();
    void testAdjointLinearity();
    
};

#endif	/* TESTOPDCT2_H */


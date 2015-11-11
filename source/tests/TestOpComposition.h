/*
 * File:   TestOpComposition.h
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 2:40:14 AM
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

#ifndef TESTOPCOMPOSITION_H
#define	TESTOPCOMPOSITION_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpComposition.h"
#include <cppunit/extensions/HelperMacros.h>

class TestOpComposition : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpComposition);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCall2);
    CPPUNIT_TEST(testCallAdjoint);
    CPPUNIT_TEST(testDimension);    

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpComposition();
    virtual ~TestOpComposition();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCall2();
    void testCallAdjoint();
    void testDimension();

};

#endif	/* TESTOPCOMPOSITION_H */


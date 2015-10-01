/*
 * File:   TestOpDCT3.h
 * Author: chung
 *
 * Created on Sep 16, 2015, 3:42:02 AM
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

#ifndef TESTOPDCT3_H
#define	TESTOPDCT3_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpDCT3.h"

#include <cppunit/extensions/HelperMacros.h>

class TestOpDCT3 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpDCT3);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testLinearity);
    CPPUNIT_TEST(testAdjointLinearity);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpDCT3();
    virtual ~TestOpDCT3();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testLinearity();
    void testAdjointLinearity();

};

#endif	/* TESTOPDCT3_H */


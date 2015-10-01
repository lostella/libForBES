/*
 * File:   TestOpGradient.h
 * Author: chung
 *
 * Created on Sep 16, 2015, 1:57:49 AM
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

#ifndef TESTOPGRADIENT_H
#define	TESTOPGRADIENT_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "OpGradient.h"

#include <cppunit/extensions/HelperMacros.h>

class TestOpGradient : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpGradient);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testLinearity);
    CPPUNIT_TEST(testAdjointLinearity);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpGradient();
    virtual ~TestOpGradient();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testLinearity();
    void testAdjointLinearity();

};

#endif	/* TESTOPGRADIENT_H */


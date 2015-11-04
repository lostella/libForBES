/*
 * File:   TestSeparableSum.h
 * Author: Pantelis Sopasakis
 *
 * Created on Nov 4, 2015, 12:10:44 AM
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

#ifndef TESTSEPARABLESUM_H
#define	TESTSEPARABLESUM_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <vector>

#include <cppunit/extensions/HelperMacros.h>

class TestSeparableSum : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestSeparableSum);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);

    CPPUNIT_TEST_SUITE_END();

public:
    TestSeparableSum();
    virtual ~TestSeparableSum();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();
    
};

#endif	/* TESTSEPARABLESUM_H */


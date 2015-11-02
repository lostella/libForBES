/*
 * File:   TestHuber.h
 * Author: Pantelis Sopasakis
 *
 * Created on Oct 30, 2015, 2:19:35 AM
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

#ifndef TESTHUBER_H
#define	TESTHUBER_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestHuber : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestHuber);

    CPPUNIT_TEST(testCall);

    CPPUNIT_TEST_SUITE_END();

public:
    TestHuber();
    virtual ~TestHuber();
    void setUp();
    void tearDown();

private:
    void testCall();

};

#endif	/* TESTHUBER_H */


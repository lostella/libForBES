/*
 * File:   TestProperties.h
 * Author: Pantelis Sopasakis
 *
 * Created on Mar 5, 2016, 4:53:48 PM
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

#ifndef TESTFBSTATS_H
#define	TESTFBSTATS_H

#define FORBES_TEST_UTILS
#include "Properties.h"
#include "ForBES.h"


#include <cppunit/extensions/HelperMacros.h>

class TestProperties : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestProperties);

    CPPUNIT_TEST(testDouble);
    CPPUNIT_TEST(testInteger);
    CPPUNIT_TEST(testVector);
    CPPUNIT_TEST(testNotExisting);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestProperties();
    virtual ~TestProperties();
    void setUp();
    void tearDown();

private:
    void testDouble();
    void testInteger();
    void testVector();
    void testNotExisting();
    
};

#endif	/* TESTFBSTATS_H */


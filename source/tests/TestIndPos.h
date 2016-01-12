/*
 * File:   TestIndPos.h
 * Author: Pantelis Sopasakis
 *
 * Created on Jan 12, 2016, 4:06:34 PM
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

#ifndef TESTINDPOS_H
#define	TESTINDPOS_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include "IndBox.h"

class TestIndPos : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestIndPos);

    CPPUNIT_TEST(testCall1);
    CPPUNIT_TEST(testCall2);
    CPPUNIT_TEST(testCall3);
    CPPUNIT_TEST(testConjugate1);
    CPPUNIT_TEST(testConjugate2);
    CPPUNIT_TEST(testConjugate3);
    CPPUNIT_TEST(testProx1);
    CPPUNIT_TEST(testProx2);
    CPPUNIT_TEST(testProx3);
    CPPUNIT_TEST(testCategory);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestIndPos();
    virtual ~TestIndPos();
    void setUp();
    void tearDown();

private:
    void testCall1();
    void testCall2();
    void testCall3();
    void testConjugate1();
    void testConjugate2();
    void testConjugate3();
    void testProx1();
    void testProx2();
    void testProx3();
    void testCategory();
    
};

#endif	/* TESTINDPOS_H */


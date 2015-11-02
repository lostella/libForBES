/*
 * File:   TestCholesky.h
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 5, 2015, 12:28:54 AM
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

#ifndef TESTCHOLESKY_H
#define	TESTCHOLESKY_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include "ForBESUtils.h"

#include <cppunit/extensions/HelperMacros.h>

class TestCholesky : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestCholesky);

    CPPUNIT_TEST(testCholeskyDense);
    CPPUNIT_TEST(testCholeskySymmetric);
    CPPUNIT_TEST(testCholeskySymmetric2);
    CPPUNIT_TEST(testCholeskySparse);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestCholesky();
    virtual ~TestCholesky();
    void setUp();
    void tearDown();

private:
    void testCholeskyDense();
    void testCholeskySymmetric();
    void testCholeskySymmetric2();
    void testCholeskySparse();
    
};

#endif	/* TESTCHOLESKY_H */

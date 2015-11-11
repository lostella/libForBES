/*
 * File:   TestMatrixExtras.h
 * Author: Pantelis Sopasakis
 *
 * Created on Nov 8, 2015, 4:33:45 PM
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

#ifndef TESTMATRIXEXTRAS_H
#define	TESTMATRIXEXTRAS_H

#define FORBES_TEST_UTILS

#include <cppunit/extensions/HelperMacros.h>

class TestMatrixExtras : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestMatrixExtras);
    
    CPPUNIT_TEST(test_add_wrong_args);
    CPPUNIT_TEST(test_mult_wrong_args);

    /* ADDITION */
    CPPUNIT_TEST(test_add_DD);
    CPPUNIT_TEST(test_add_DDT);
    CPPUNIT_TEST(test_add_DTD);
    CPPUNIT_TEST(test_add_DTDT);
    
    CPPUNIT_TEST(test_add_DS);
    CPPUNIT_TEST(test_add_DST);
    
    CPPUNIT_TEST(test_add_SS);
    CPPUNIT_TEST(test_add_SST);
    CPPUNIT_TEST(test_add_STS);
    CPPUNIT_TEST(test_add_STST);
    
    
    
    /* MULTIPLICATION */    
    CPPUNIT_TEST(test_mult_DD);
    CPPUNIT_TEST(test_mult_DDT);
    CPPUNIT_TEST(test_mult_DX);
    CPPUNIT_TEST(test_mult_DH);
    
    CPPUNIT_TEST(test_mult_SS);
    CPPUNIT_TEST(test_mult_SST);
    CPPUNIT_TEST(test_mult_SS2);
    CPPUNIT_TEST(test_mult_SS3);
    
    CPPUNIT_TEST(test_mult_SX);

    CPPUNIT_TEST_SUITE_END();

public:
    TestMatrixExtras();
    virtual ~TestMatrixExtras();
    void setUp();
    void tearDown();

private:
    void test_add_wrong_args();
    void test_mult_wrong_args();
    
    // dense with dense: fully tested!
    void test_add_DD(); 
    void test_add_DDT();
    void test_add_DTD();
    void test_add_DTDT();
    
    void test_add_DS();
    void test_add_DST();
    
    // sparse with sparse: fully tested!
    void test_add_SS();
    void test_add_SST();
    void test_add_STS();
    void test_add_STST();
    
    
    
    void test_mult_DD();
    void test_mult_DDT();
    void test_mult_DX();
    void test_mult_DH();
    
    void test_mult_SS();
    void test_mult_SST();
    void test_mult_SS2();
    void test_mult_SS3();
    
    void test_mult_SX();
};

#endif	/* TESTMATRIXEXTRAS_H */


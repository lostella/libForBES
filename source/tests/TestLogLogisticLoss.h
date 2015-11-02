/*
 * File:   TestLogLogisticLoss.h
 * Author: chung
 *
 * Created on Oct 29, 2015, 7:38:05 PM
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

#ifndef TESTLOGLOGISTICLOSS_H
#define	TESTLOGLOGISTICLOSS_H

#define FORBES_TEST_UTILS
#include "ForBES.h"


#include <cppunit/extensions/HelperMacros.h>

class TestLogLogisticLoss : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestLogLogisticLoss);

    CPPUNIT_TEST(testCall);

    CPPUNIT_TEST_SUITE_END();

public:
    TestLogLogisticLoss();
    virtual ~TestLogLogisticLoss();
    void setUp();
    void tearDown();

private:
    void testCall();

};

#endif	/* TESTLOGLOGISTICLOSS_H */


/*
 * File:   TestSeparableSum.cpp
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

#include "TestSeparableSum.h"
#include "ElasticNet.h"
#include "HingeLoss.h"
#include "MatrixFactory.h"
#include "SeparableSum.h"
#include "ForBES.h"



CPPUNIT_TEST_SUITE_REGISTRATION(TestSeparableSum);

TestSeparableSum::TestSeparableSum() {
}

TestSeparableSum::~TestSeparableSum() {
}

void TestSeparableSum::setUp() {
}

void TestSeparableSum::tearDown() {
}

void TestSeparableSum::testCall() {
    size_t num_functs = 2; /* 2 functions */
    size_t num_idx_1 = 3; /* Function #1 points to 3 indices */
    size_t num_idx_2 = 2; /* Function #2 points to 2 indices */

    std::vector<size_t> idx1(num_idx_1);
    std::vector<size_t> idx2(num_idx_2);

    idx1[0] = 0;
    idx1[1] = 3;
    idx1[2] = 4; /* f1(x_0, x_3, x_4) */

    idx2[0] = 1;
    idx2[1] = 2; /* f2(x_1, x_2) */

    Function * f1 = new ElasticNet(2.5, 1.3);

    Matrix b(2, 1);
    b[0] = 0.8;
    b[1] = 1.5;
    Function * f2 = new HingeLoss(b, 1.4);

    std::map<Function*, std::vector<size_t>* > fun_idx_map;
    fun_idx_map[f1] = &idx1;
    fun_idx_map[f2] = &idx2;

    Function * sep_sum = new SeparableSum(fun_idx_map);

    Matrix x(5, 1);
    for (size_t i = 0; i < 5; i++) {
        x[i] = 0.1 * (i + 1);
    }


    double f_val;
    _ASSERT(sep_sum->category().defines_f());
    int status = sep_sum->call(x, f_val);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    Matrix x1(num_idx_1, 1);
    x1[0] = x[0];
    x1[1] = x[3];
    x1[2] = x[4];

    Matrix x2(num_idx_2, 1);
    x2[0] = x[1];
    x2[1] = x[2];

    double f1_val;
    status = f1->call(x1, f1_val);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    double f2_val;
    status = f2->call(x2, f2_val);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    _ASSERT_NUM_EQ(f1_val + f2_val, f_val, 1e-8);

    delete f1;
    delete f2;
    delete sep_sum;

}

void TestSeparableSum::testCallProx() {
    size_t num_functs = 2; /* 2 functions */
    size_t num_idx_1 = 3; /* Function #1 points to 3 indices */
    size_t num_idx_2 = 2; /* Function #2 points to 2 indices */

    std::vector<size_t> idx1(num_idx_1);
    std::vector<size_t> idx2(num_idx_2);

    idx1[0] = 0;
    idx1[1] = 3;
    idx1[2] = 4; /* f1(x_0, x_3, x_4) */

    idx2[0] = 1;
    idx2[1] = 2; /* f2(x_1, x_2) */

    Function * f1 = new ElasticNet(2.5, 1.3);

    Matrix b(2, 1);
    b[0] = 0.8;
    b[1] = 1.5;
    Function * f2 = new HingeLoss(b, 1.4);

    std::map<Function*, std::vector<size_t>* > fun_idx_map;
    fun_idx_map[f1] = &idx1;
    fun_idx_map[f2] = &idx2;

    Function * sep_sum = new SeparableSum(fun_idx_map);

    Matrix x(5, 1);
    for (size_t i = 0; i < 5; i++) {
        x[i] = 10.9 * (i + 1);
    }


    double gamma = 0.87;
    Matrix prox(num_idx_1 + num_idx_2, 1);
    _ASSERT(sep_sum->category().defines_prox());
    int status = sep_sum->callProx(x, gamma, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);


    Matrix x1(num_idx_1, 1);
    x1[0] = x[0];
    x1[1] = x[3];
    x1[2] = x[4];

    Matrix x2(num_idx_2, 1);
    x2[0] = x[1];
    x2[1] = x[2];

    Matrix prox1(num_idx_1, 1);
    Matrix prox2(num_idx_2, 1);

    status = f1->callProx(x1, gamma, prox1);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    status = f2->callProx(x2, gamma, prox2);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    double tol = 1e-8;
    _ASSERT_NUM_EQ(prox1[0], prox[0], tol);
    _ASSERT_NUM_EQ(prox2[0], prox[1], tol);
    _ASSERT_NUM_EQ(prox2[1], prox[2], tol);
    _ASSERT_NUM_EQ(prox1[1], prox[3], tol);
    _ASSERT_NUM_EQ(prox1[2], prox[4], tol);
   

    delete f1;
    delete f2;
    delete sep_sum;
}




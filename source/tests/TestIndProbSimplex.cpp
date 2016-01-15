/*
 * File:   TestIndProbSimplex.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Jan 15, 2016, 5:38:56 PM
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

#include "TestIndProbSimplex.h"
#include "IndProbSimplex.h"
#include "MatrixFactory.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestIndProbSimplex);

TestIndProbSimplex::TestIndProbSimplex() {
}

TestIndProbSimplex::~TestIndProbSimplex() {
}

void TestIndProbSimplex::setUp() {
}

void TestIndProbSimplex::tearDown() {
}

void TestIndProbSimplex::testCallProx() {
    const size_t n = 10;

    Function * F = new IndProbSimplex();

    double x_vals[n] = {0.910569988523029, -0.300000000000000,
        -0.500000000000000, 0.813112813610761, 0.383306318625529, 0.617279232316449,
        0.575494859702814, 0.530051704765016, 0.275069755821935, 0.248628959661970};

    Matrix x(n, 1, x_vals);
    Matrix x_orig(x);
    Matrix prox(n, 1);
    int status = F->callProx(x, 1.0, prox);
    _ASSERT(ForBESUtils::is_status_ok(status));

    double prox_expected_vals[n] = {0.421268268739415, 0, 0, 0.323811093827147, 0,
        0.127977512532836, 0.086193139919200, 0.040749984981402, 0, 0};
    Matrix prox_expected = MatrixFactory::ShallowVector(prox_expected_vals, n, 0);

    double f_at_prox = 2.0;
    status = F->callProx(x, 1.0, prox, f_at_prox);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f_at_prox, 1e-14);   
    _ASSERT_EQ(prox_expected, prox);
    _ASSERT_EQ(x_orig, x);
    delete F;

}

void TestIndProbSimplex::testCallProxLarge() {
    std::srand(10);
    const size_t n = 1e6;
    Function * F = new IndProbSimplex();
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);
    Matrix prox(n, 1);

//    clock_t start = clock();
    int status = F->callProx(x, 1.0, prox);
//    clock_t end = clock();
//    float elapsed_time_secs = static_cast<float>(end - start)*1000.0 / CLOCKS_PER_SEC;
    _ASSERT(ForBESUtils::is_status_ok(status));
    
//    std::cout << "\ntime = " << elapsed_time_secs << std::endl;
    delete F;
}



/*
 * File:   TestLogLogisticLoss.cpp
 * Author: Pantelis Sopasakis
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

#include "TestLogLogisticLoss.h"
#include "LogLogisticLoss.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestLogLogisticLoss);

TestLogLogisticLoss::TestLogLogisticLoss() {
}

TestLogLogisticLoss::~TestLogLogisticLoss() {
}

void TestLogLogisticLoss::setUp() {
}

void TestLogLogisticLoss::tearDown() {
}

void TestLogLogisticLoss::testCall() {
    const double mu = 1.5;
    Function *logLogisticLoss = new LogLogisticLoss(mu);
    const size_t n = 10;
    const double xdata[n] = {
        0.537667139546100,
        1.833885014595086,
        -2.258846861003648,
        0.862173320368121,
        0.318765239858981,
        -1.307688296305273,
        -0.433592022305684,
        0.342624466538650,
        3.578396939725760,
        2.769437029884877
    };
    const double ddata[n] = {
        -0.013569096140251,
        -0.333434615170239,
        -0.431122888787024,
        0.407923343445659,
        0.404836067818807,
        -0.257685527804887,
         -0.057441551301617,
         -0.101942255039910,
         -0.226500118101661,
          0.490999700556483
    };

    Matrix x(n, 1, xdata);
    Matrix d(n, 1, ddata);

    double f;
    Matrix grad(n, 1);
    Matrix Hd(n, 1);

    const double tol = 1e-12;
    const double f_expected = 10.455339600285200;

    const double grad_expected_data[n] = {
        -0.553095647613005,
        -0.206664161288562,
        -1.358116380141062,
        -0.445328215242575,
        -0.631465046697517,
        -1.180689100643326,
        -0.910096624790372,
        -0.622758155942074,
        -0.040743067205105,
        -0.088497390570449
    };

    const double Hd_expected_data[n] = {
        -0.004737683170800,
        -0.059414973349659,
        -0.055383330523630,
         0.127727625721715,
         0.148021416496942,
        -0.064766147765657,
        -0.020559061512420,
        -0.037128015811522,
        -0.008977649776764,
         0.040888588516473,
    };

    Matrix grad_expected(n, 1, grad_expected_data);
    Matrix Hd_expected(n, 1, Hd_expected_data);

    _ASSERT(logLogisticLoss->category().defines_f());
    int status = logLogisticLoss->call(x, f, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_expected, f, tol);
    _ASSERT_EQ(grad_expected, grad);

    status = logLogisticLoss->hessianProduct(x, d, Hd);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(Hd_expected, Hd);

    double f2;
    status = logLogisticLoss->call(x, f2);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f, f2, tol);

    x.set(1, 0, 40);
    x.set(2, 0, 34);
    
    double f3;
    status = logLogisticLoss->call(x, f3);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(6.695659275314881, f3, tol);
    
    double f4;
    x.set(0, 0, -70);
    status = logLogisticLoss->call(x, f4);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(1.110056258258362e+02, f4, tol);

    delete logLogisticLoss;

}


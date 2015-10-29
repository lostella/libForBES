/*
 * File:   TestLogLogisticLoss.cpp
 * Author: chung
 *
 * Created on Oct 29, 2015, 7:38:05 PM
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
    Matrix x(n, 1, xdata);
    double f;
    Matrix grad(n, 1);

    const double tol = 1e-10;
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

    Matrix grad_expected(n, 1, grad_expected_data);

    int status = logLogisticLoss->call(x, f, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_expected, f, tol);
    _ASSERT_EQ(grad_expected, grad);


    double f2;
    status = logLogisticLoss->call(x, f2);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f, f2, tol);
    
    delete logLogisticLoss;

}


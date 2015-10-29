/*
 * File:   TestQuadraticLoss.cpp
 * Author: chung
 *
 * Created on Oct 29, 2015, 7:32:54 PM
 */

#include "TestQuadraticLoss.h"
#include "QuadraticLoss.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestQuadraticLoss);

TestQuadraticLoss::TestQuadraticLoss() {
}

TestQuadraticLoss::~TestQuadraticLoss() {
}

void TestQuadraticLoss::setUp() {
}

void TestQuadraticLoss::tearDown() {
}

void TestQuadraticLoss::testCall() {
    const size_t n = 10;
    const double weights_data[n] = {
        0.157613081677548,
        0.970592781760616,
        0.957166948242946,
        0.485375648722841,
        0.800280468888800,
        0.141886338627215,
        0.421761282626275,
        0.915735525189067,
        0.792207329559554,
        0.959492426392903
    };
    const double p_data[n] = {
        2.665186895272925,
        -3.441210320907452,
        -3.206611374504095,
        -2.428496083274627,
        -8.832852485984688,
        4.315140878445295,
        0.975571618368594,
        -2.264784957509110,
        4.110895620285683,
        -5.134549256561094,
    };
    Matrix w(n, 1, weights_data);
    Matrix p(n, 1, p_data);

    Function *quadLoss = new QuadraticLoss(&w, &p);

    const double xdata[n] = {
        -0.511212230427454,
        -1.207235208036790,
        1.596033695827509,
        1.564292983187142,
        -4.324399586622282,
        -0.150256480981343,
        -0.824395096045192,
        3.138536437643633,
        5.466328345197420,
        5.546366488071987
    };
    Matrix x(n, 1, xdata);

    double f;
    const double f_expected = 97.181889336760406;
    const double tol = 1e-8;
    int status = quadLoss->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_expected, f, tol);

    delete quadLoss;

}

void TestQuadraticLoss::testCallConj() {
//    if (true /*check result*/) {
//        CPPUNIT_ASSERT(false);
//    }
}


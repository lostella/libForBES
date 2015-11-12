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

    Function *quadLoss = new QuadraticLoss(w, p);

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
    Matrix grad(n, 1);

    double f0;
    double f;
    const double f_expected = 97.181889336760406;
    const double tol = 1e-8;
    _ASSERT(quadLoss->category().defines_f());
    int status = quadLoss->call(x, f0, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    status = quadLoss->call(x, f, grad);
    _ASSERT_EQ(f, f0);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_expected, f, tol);

    const double grad_expected_data[n] = {
        -0.500642054839506,
        2.168280119185122,
        4.596933125463330,
        1.938002583347348,
        3.608026800264816,
        -0.633578881844680,
        -0.759156270155760,
        4.948013355555519,
        1.073783739399957,
        10.248257763916156
    };
    Matrix grad_expected(n, 1, grad_expected_data);

    _ASSERT_EQ(grad_expected, grad);

    double f_2;
    status = quadLoss->call(x, f_2);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f, f_2, tol);

    delete quadLoss;
}

void TestQuadraticLoss::testCallConj() {
    const size_t n = 10;
    const double weights_data[n] = {
        8.407172559836626,
        2.542821789715310,
        8.142848260688163,
        2.435249687249893,
        9.292636231872278,
        3.499837659848088,
        1.965952504312082,
        2.510838579760311,
        6.160446761466392,
        4.732888489027292
    };
    const double p_data[n] = {
        -11.658439314820487,
        -11.479527788985941,
        1.048747160164940,
        7.222540322250016,
        25.854912526162416,
        -6.668906707013855,
        1.873310245789398,
        -0.824944253709554,
        -19.330229178509867,
        -4.389661539347733
    };
    Matrix w(n, 1, weights_data);
    Matrix p(n, 1, p_data);

    Function *quadLoss = new QuadraticLoss(w, p);

    const double xdata[n] = {
        -3.589357682910247,
        1.680751059507811,
        -1.776064164658021,
        0.200185666278645,
        -1.089057859981095,
        0.607041589298709,
        -1.200653124267468,
        0.979930642347896,
        1.478726247208948,
        3.423775565963109
    };
    Matrix x(n, 1, xdata);

    Matrix grad(n, 1);
    double f_star;
    double f_star_expected = -53.127646012575340;
    const double tol = 1e-8;
    _ASSERT(quadLoss->category().defines_conjugate());
    _ASSERT(quadLoss->category().defines_conjugate_grad());
    int status = quadLoss->callConj(x, f_star, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_star_expected, f_star, tol);

    const double grad_expected_data[n] = {
        -12.085379247046644,
        -10.818549081667255,
        0.830633779220282,
        7.304743665638888,
        25.737716734483215,
        -6.495458207884025,
        1.262586883127290,
        -0.434664030103367,
        -19.090193625170070,
        -3.666260687225234
    };
    Matrix grad_expected(n, 1, grad_expected_data);

    _ASSERT_EQ(grad_expected, grad);

    double f_star2;
    status = quadLoss->callConj(x, f_star2);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_star, f_star2, tol);

    delete quadLoss;
}


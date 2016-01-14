/*
 * File:   TestDistanceToBall2.cpp
 * Author: chung
 *
 * Created on Jan 14, 2016, 4:32:13 PM
 */

#include "TestDistanceToBall2.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestDistanceToBall2);

TestDistanceToBall2::TestDistanceToBall2() {
}

TestDistanceToBall2::~TestDistanceToBall2() {
}

void TestDistanceToBall2::setUp() {
}

void TestDistanceToBall2::tearDown() {
}

void TestDistanceToBall2::testCall() {
    Function * d2b = new DistanceToBall2();
    Matrix x(2, 1);
    x[0] = 1.0;
    x[1] = 1.0;

    double b = std::pow(M_SQRT2 - 1, 2) / 2.0;
    double f;
    int status = d2b -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(b, f, 1e-10);

    double y_values[10] = {
        0.706046088019609,
        0.031832846377421,
        0.276922984960890,
        0.046171390631154,
        0.097131781235848,
        0.823457828327293,
        0.694828622975817,
        0.317099480060861,
        0.950222048838355,
        0.034446080502909
    };
    Matrix y(10, 1, y_values);
    b = 0.217342378005960;
    status = d2b -> call(y, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(b, f, 1e-8);

    delete d2b;
}

void TestDistanceToBall2::testCall2() {
    Function * d2b = new DistanceToBall2(5.4312);
    double x_values[10] = {
        0.438744359656398,
        0.381558457093008,
        0.765516788149002,
        0.795199901137063,
        0.186872604554379,
        0.489764395788231,
        0.445586200710899,
        0.646313010111265,
        0.709364830858073,
        0.754686681982361
    };
    Matrix x(10, 1, x_values);
    double expected = 2.084995711340361;
    double f;
    int status = d2b->call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(expected, f, 1e-8);
    delete d2b;
}

void TestDistanceToBall2::testCall3() {
    Function * d2b = new DistanceToBall2(0.998, 2.54);
    double x_values[10] = {
        0.877488719312796,
        0.763116914186016,
        1.531033576298004,
        1.590399802274126,
        0.373745209108758,
        0.979528791576462,
        0.891172401421798,
        1.292626020222530,
        1.418729661716146,
        1.509373363964722
    };
    Matrix x(10, 1, x_values);
    double expected = 0.733567206618426;
    double f;
    int status = d2b->call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(expected, f, 1e-8);
    delete d2b;
}

void TestDistanceToBall2::testCall4() {
    double x_values[10] = {
        0.877488719312796,
        0.763116914186016,
        1.531033576298004,
        1.590399802274126,
        0.373745209108758,
        0.979528791576462,
        0.891172401421798,
        1.292626020222530,
        1.418729661716146,
        1.509373363964722
    };
    Matrix x(10, 1, x_values);

    double c_values[10] = {
        -0.863652821988714,
        0.077359091130425,
        -1.214117043615409,
        -1.113500741486764,
        -0.006849328103348,
        1.532630308284750,
        -0.769665913753682,
        0.371378812760058,
        -0.225584402271252,
        1.117356138814467
    };


    Matrix c_shallow = MatrixFactory::ShallowVector(c_values, 10, 0);
    Function * d2b = new DistanceToBall2(0.998, 2.54, c_shallow);

    double expected = 3.084754993825916;
    double f;
    int status = d2b->call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(expected, f, 1e-8);
    delete d2b;
}

void TestDistanceToBall2::testCall5() {
    double x_values[10] = {
        0.877488719312796,
        0.763116914186016,
        1.531033576298004,
        1.590399802274126,
        0.373745209108758,
        0.979528791576462,
        0.891172401421798,
        1.292626020222530,
        1.418729661716146,
        1.509373363964722
    };
    Matrix x(10, 1, x_values);
    Matrix c(x);
    double rho = 3.23;
    double w = 0.887;
    for (size_t i = 0; i < 10; i++) {
        c[i] += rho / 11.0;
    }

    Function * d2b = new DistanceToBall2(w, rho, c);
    double expected = 0.0;
    double f;
    int status = d2b->call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(expected, f, 1e-8);
    delete d2b;
}

void TestDistanceToBall2::testGradient() {
    Function * d2b = new DistanceToBall2(5.4312);
    double x_values[10] = {
        0.438744359656398,
        0.381558457093008,
        0.765516788149002,
        0.795199901137063,
        0.186872604554379,
        0.489764395788231,
        0.445586200710899,
        0.646313010111265,
        0.709364830858073,
        0.754686681982361
    };
    Matrix x(10, 1, x_values);
    Matrix grad(10, 1);

    double grad_expected_values[10] = {
        1.112859194860052,
        0.967809221946716,
        1.941705637375089,
        2.016995779558804,
        0.473995600555761,
        1.242269670645696,
        1.130213277162638,
        1.639349567076561,
        1.799278229238064,
        1.914235464907696

    };
    Matrix grad_expected(10, 1, grad_expected_values);
    double expected = 2.084995711340361;
    double f;
    int status = d2b->call(x, f, grad);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(expected, f, 1e-8);
    _ASSERT_EQ(grad_expected, grad);

    delete d2b;
}

void TestDistanceToBall2::testGradient2() {
    double x_values[10] = {
        0.438744359656398,
        0.381558457093008,
        0.765516788149002,
        0.795199901137063,
        0.186872604554379,
        0.489764395788231,
        0.445586200710899,
        0.646313010111265,
        0.709364830858073,
        0.754686681982361
    };
    Matrix x(10, 1, x_values);

    double c_values[10] = {
        -1.089064295052236,
        0.032557464164973,
        0.552527021112224,
        1.100610217880866,
        1.544211895503951,
        0.085931133175425,
        -1.491590310637609,
        -0.742301837259857,
        -1.061581733319986,
        2.350457224002042
    };
    Matrix c_shallow = MatrixFactory::ShallowVector(c_values, 10, 0);

    double grad_expected_values[10] = {
        0.555435626832365,
        0.126879491534933,
        0.077432539996649,
        -0.111032078656484,
        -0.493461401445295,
        0.146813791546350,
        0.704261522900583,
        0.504831646162784,
        0.643828539607569,
        -0.580143206124434

    };
    Matrix grad_expected(10, 1, grad_expected_values);

    Function * d2b = new DistanceToBall2(0.998, 2.54, c_shallow);
    Matrix grad(10, 1);
    double expected = 1.057069197096058;
    double f;
    int status = d2b->call(x, f, grad);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(expected, f, 1e-8);
    _ASSERT_EQ(grad_expected, grad);
    delete d2b;
}

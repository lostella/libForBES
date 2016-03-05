/*
 * File:   TestFBStats.cpp
 * Author: chung
 *
 * Created on Mar 5, 2016, 4:53:48 PM
 */

#include "TestFBStats.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestFBStats);

TestFBStats::TestFBStats() {
}

TestFBStats::~TestFBStats() {
}

void TestFBStats::setUp() {
}

void TestFBStats::tearDown() {
}

void TestFBStats::testDouble() {
    FBStats stats;

    const double first_value = 34.59210;
    const double second_value = -52.132;
    double *a = new double;
    *a = first_value;

    /* Add a value to FBStats uder `key` */
    std::string key = "key-1";
    bool status = stats.add_property<double>(key, FBStats::DOUBLE_TYPE, a);
    _ASSERT(status);

    /* Retrieve the value under `key` */
    int type; // type
    void * value; // actual value
    status = stats.get_typed_property(key, type, value);
    _ASSERT(status); // make sure the key-value-pair exists
    _ASSERT_EQ(FBStats::DOUBLE_TYPE, type); // make sure the retrieved value is of type double

    double * dbl_val = static_cast<double*> (value); // cast as double*
    _ASSERT_EQ(a, dbl_val);
    _ASSERT_EQ(first_value, *dbl_val);
    _ASSERT_EQ(first_value, *a);

    /* Now add again the same key (update) */
    *a = second_value;
    status = stats.add_property<double>(key, FBStats::DOUBLE_TYPE, a);
    _ASSERT_NOT(status);

    status = stats.get_typed_property(key, type, value);
    dbl_val = static_cast<double*> (value); // cast as double*
    _ASSERT(status); // make sure the key-value-pair exists
    _ASSERT_EQ(FBStats::DOUBLE_TYPE, type); // make sure the retrieved value is of type double
    _ASSERT_EQ(second_value, *dbl_val);

    delete a;
}

void TestFBStats::testInteger() {
    FBStats * stats = new FBStats();
    int i = 365;
    bool how = stats->add_property<int>("i", FBStats::INT_TYPE, &i);
    _ASSERT(how);
    how = stats->add_property<int>("i", FBStats::INT_TYPE, &i);
    _ASSERT_NOT(how);

    int * i_retrieved;
    how = stats->get_property<>("i", i_retrieved);
    _ASSERT(how);
    _ASSERT_NEQ(NULL, i_retrieved);
    _ASSERT_EQ(i_retrieved, &i);
    _ASSERT_EQ(*i_retrieved, i);

    how = stats->get_property<>("asdf-qwerty", i_retrieved);
    _ASSERT_NOT(how);

    delete stats;
}

void TestFBStats::testNotExisting() {
    FBStats * stats = new FBStats();

    double a = 345.678;
    stats->add_property<>("a", FBStats::DOUBLE_TYPE, &a);

    size_t q = 122;
    stats->add_property("q", FBStats::SIZE_T_TYPE, &q);

    unsigned long s = 4564564;
    stats->add_property("s", FBStats::CUSTOM_TYPE, &s);

    void * val = NULL;
    bool found = stats->get_property("chung", val);
    _ASSERT_NOT(found);
    _ASSERT_EQ(static_cast<void*> (NULL), val);

    _ASSERT_EQ(static_cast<size_t> (3), stats->size());


    delete stats;
}

void TestFBStats::testVector() {
    FBStats * stats = new FBStats();

    std::vector<double> vec;
    vec.push_back(1.23);
    vec.push_back(2.34);
    vec.push_back(3.45);

    string key = "key-vec-1";

    stats->add_property<std::vector<double> >(key, FBStats::VECTOR_TYPE, &vec);

    std::vector<double> * found = NULL;
    bool found_status = stats->get_property<std::vector<double> >(key, found);
    _ASSERT(found_status);
    _ASSERT_NEQ(NULL, found);
    for (size_t i = 0; i < vec.size(); i++) {
        _ASSERT_EQ(vec[i], found->at(i));
    }

    delete stats;
}





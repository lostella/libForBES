/*
 * File:   TestProperties.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Mar 5, 2016, 4:53:48 PM
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

#include "TestProperties.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestProperties);

TestProperties::TestProperties() {
}

TestProperties::~TestProperties() {
}

void TestProperties::setUp() {
}

void TestProperties::tearDown() {
}

void TestProperties::testDouble() {
    Properties stats;

    const double first_value = 34.59210;
    const double second_value = -52.132;
    double *a = new double;
    *a = first_value;

    /* Add a value to Properties under `key` */
    std::string key = "key-1";
    bool status = stats.add_property<double>(key, Properties::DOUBLE_TYPE, a);
    _ASSERT(status);

    /* Retrieve the value under `key` */
    int type; // type
    void * value; // actual value
    status = stats.get_typed_property(key, type, value);
    _ASSERT(status); // make sure the key-value-pair exists
    _ASSERT_EQ(Properties::DOUBLE_TYPE, type); // make sure the retrieved value is of type double

    double * dbl_val = static_cast<double*> (value); // cast as double*
    _ASSERT_EQ(a, dbl_val);
    _ASSERT_EQ(first_value, *dbl_val);
    _ASSERT_EQ(first_value, *a);

    /* Now add again the same key (update) */
    *a = second_value;
    status = stats.add_property<double>(key, Properties::DOUBLE_TYPE, a);
    _ASSERT_NOT(status);

    status = stats.get_typed_property(key, type, value);
    dbl_val = static_cast<double*> (value); // cast as double*
    _ASSERT(status); // make sure the key-value-pair exists
    _ASSERT_EQ(Properties::DOUBLE_TYPE, type); // make sure the retrieved value is of type double
    _ASSERT_EQ(second_value, *dbl_val);

    delete a;
}

void TestProperties::testInteger() {
    Properties * stats = new Properties();
    int i = 365;
    bool how = stats->add_property<int>("i", Properties::INT_TYPE, &i);
    _ASSERT(how);
    how = stats->add_property<int>("i", Properties::INT_TYPE, &i);
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

void TestProperties::testNotExisting() {
    Properties * stats = new Properties();

    double a = 345.678;
    stats->add_property<>("a", Properties::DOUBLE_TYPE, &a);

    size_t q = 122;
    stats->add_property("q", Properties::SIZE_T_TYPE, &q);

    unsigned long s = 4564564;
    stats->add_property("s", Properties::CUSTOM_TYPE, &s);

    void * val = NULL;
    bool found = stats->get_property("chung", val);
    _ASSERT_NOT(found);
    _ASSERT_EQ(static_cast<void*> (NULL), val);

    _ASSERT_EQ(static_cast<size_t> (3), stats->size());


    delete stats;
}

void TestProperties::testVector() {
    Properties * stats = new Properties();

    std::vector<double> vec;
    vec.push_back(1.23);
    vec.push_back(2.34);
    vec.push_back(3.45);

    string key = "key-vec-1";

    stats->add_property<std::vector<double> >(key, Properties::VECTOR_TYPE, &vec);

    std::vector<double> * found = NULL;
    bool found_status = stats->get_property<std::vector<double> >(key, found);
    _ASSERT(found_status);
    _ASSERT_NEQ(NULL, found);
    for (size_t i = 0; i < vec.size(); i++) {
        _ASSERT_EQ(vec[i], found->at(i));
    }
    
    
    /* here, let's assume we don't know the type of the value of `key` */
    int type_found = 0;
    void * obj;
    found_status = stats->get_typed_property(key, type_found, obj);
    _ASSERT(found_status);
    _ASSERT_EQ(Properties::VECTOR_TYPE, type_found);
    _ASSERT_NEQ(NULL, found);
    
    std::vector<double>* obj_as_vec = static_cast<std::vector<double>*>(obj); 
    obj_as_vec->push_back(10.34);
    obj_as_vec->push_back(11.01);
    obj_as_vec->push_back(-8.75);
    
    _ASSERT_EQ(static_cast<size_t>(6), vec.size());
    _ASSERT_EQ(-8.75, vec[5]);

    delete stats;
}





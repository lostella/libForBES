/* 
 * File:   Properties.h
 * Author: Pantelis Sopasakis
 *
 * Created on March 4, 2016, 4:07 PM
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

#ifndef FORBES_PROPERTIES_H
#define	FORBES_PROPERTIES_H

#include <stdlib.h>
#include <utility>
#include <map>
#include <vector>
#include <string>

/**
 * \class Properties
 * \brief Generic properties class
 * \version version 0.1
 * \date Created on March 4, 2016, 4:07 PM
 * \author Pantelis Sopasakis
 * 
 * 
 * This class can be used to provide options to solvers as well as to 
 * get generic results and statistics from the solver's output.
 */
class Properties {
public:
    Properties();
    virtual ~Properties();


    const static int DOUBLE_TYPE; /**< <code>double</code> type */
    const static int INT_TYPE; /**< <code>size_t</code> type  */
    const static int SIZE_T_TYPE; /**< <code>size_t</code> type  */
    const static int STRING_TYPE; /**< <code>std::string</code> type  */
    const static int VECTOR_TYPE; /**< type <code>std::vector</code> type  */
    const static int CUSTOM_TYPE; /**< any custom type not listed here */

    /**
     * 
     * Adds a new key-value pair.
     * 
     * Example of use:
     * 
     * \code{.cpp}
     *  Properties stats;
     *  double a = 34.11;
     *  std::string key = "key-1";
     *  bool status = stats.add_property<double>(key, Properties::DOUBLE_TYPE, &a);
     * \endcode
     * 
     * @param key key as a string 
     * @param type data-type
     * @param value pointer to where the data are stored 
     * @return returns <code>true</code> iff the key-value pair did not pre-exist
     * and it was created now for the first time.
     */
    template<typename T> bool add_property(std::string key, int type, T* value) {
        TypedPair pair = std::make_pair<int, void*>(type, value);
        std::pair < PropertiesMap::iterator, bool> insertion_status =
                m_map->insert(std::make_pair<std::string, TypedPair>(key, pair));
        return insertion_status.second;
    }

    /**
     * 
     * Returns the value stored under a key (if found) as well as its type.
     * 
     * This method is suitable in case you are unsure about the type of the stored
     * value.
     * 
     * Here is an example of use:
     * 
     * \code{.cpp}
     *  Properties stats;
     * 
     *  // Add data in Properties 
     *  double a = 34.11;
     *  std::string key = "key-1";
     *  bool status = stats.add_property<double>(key, Properties::DOUBLE_TYPE, &a);
     * 
     *  // Look-up data 
     *  int retrieved_type;
     *  void * retrieved_value;
     *  bool found = stats.get_typed_property(key, retrieved_type, retrieved_value);
     * 
     *  // Make sure that the data-type is indeed Properties::DOUBLE_TYPE
     *  if (found && retrieved_type == Properties::DOUBLE_TYPE) {
     *    double * dbl_data = static_cast<double*>(retrieved_value);
     *  }
     * \endcode
     * 
     * @param key property name
     * @param type property type (output value)
     * @param value pointer to where the data can be found 
     * @return returns <code>true</code> iff the value is found
     */
    bool get_typed_property(std::string key, int& type, void*& value);

    /**
     * 
     * In case the data-type is known, this method can be used to easily retrieve 
     * the stored value under a given key.
     * 
     * Example of use:
     * 
     * \code{.cpp}
     *  Properties stats;
     * 
     *  // Add data in Properties 
     *  std::vector<double> x;
     *  x.push_back(1.23);
     *  x.push_back(2.34);
     * 
     *  std::string key = "key-1";
     *  bool status = stats.add_property<>(key, Properties::VECTOR_TYPE, &x);
     * 
     *  // Look-up data 
     *  std::vector<double> * retrieved_x;
     *  bool found = stats.get_property<>(retrieved_x);
     * \endcode
     * 
     * 
     * @param key property name
     * @param value value (if found)
     * @return returns <code>true</code> iff the value is found
     */
    template<typename T> bool get_property(std::string key, T*& value) {
        PropertiesMap::iterator value_it = m_map->find(key);
        if (value_it != m_map->end()) {
            TypedPair value_found = value_it->second;
            /* 
             * here, static_cast may result in erroneous interpretation of (void*) 
             * if the method is not used properly. The use of:
             *   bool get_typed_property(std::string key, int& type, void*& value);
             * is recommended when the return type is unknown.
             */
            value = static_cast<T*> (value_found.second);
            return true;
        }
        return false;
    }

    /**
     * Returns the number of key-value pairs stored.
     * @return size of internal map.
     */
    size_t size();

private:

    /**
     * A typed pair (type, data)
     */
    typedef std::pair<int, void*> TypedPair;
    /**
     * A map used to store the data
     */
    typedef std::map<std::string, TypedPair> PropertiesMap;
    /**
     * An entry of the above map
     */
    typedef std::pair<std::string, TypedPair> PropertiesPair;

    void init();

    PropertiesMap * m_map;



};

#endif	/* Properties_H */


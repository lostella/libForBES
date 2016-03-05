/* 
 * File:   Properties.cpp
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


#include "Properties.h"

const int Properties::DOUBLE_TYPE = 1;
const int Properties::INT_TYPE = 2;
const int Properties::SIZE_T_TYPE = 3;
const int Properties::STRING_TYPE = 4;
const int Properties::VECTOR_TYPE = 100;
const int Properties::CUSTOM_TYPE = 666;

Properties::Properties() {
    init();
}

Properties::~Properties() {
    if (m_map != NULL) {
        delete m_map;
        m_map = NULL;
    }
}

void Properties::init() {
    m_map = new PropertiesMap;
}


bool Properties::get_typed_property(std::string key, int& type, void*& value) {
    PropertiesMap::iterator value_it = m_map->find(key);
    if (value_it != m_map->end()) {
        TypedPair value_found = value_it->second;
        type = value_found.first;
        value = value_found.second;
        return true;
    }    
    return false;
}

size_t Properties::size() {
    return m_map->size();
}






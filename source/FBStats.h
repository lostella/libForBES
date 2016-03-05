/* 
 * File:   FBStatscpp
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

#ifndef FBSTATS_H
#define	FBSTATS_H

#include <stdlib.h>
#include <utility>
#include <map>
#include <vector>
#include <string>

class FBStats {
public:
    FBStats();
    virtual ~FBStats();

    const static int DOUBLE_TYPE;
    const static int INT_TYPE;
    const static int SIZE_T_TYPE;
    const static int STRING_TYPE;
    const static int CUSTOM_TYPE;
    const static int VECTOR_TYPE;


    typedef std::pair<int, void*> TypedPair;
    typedef std::map<std::string, TypedPair> FBStatsMap;
    typedef std::pair<std::string, TypedPair> FBStatsPair;

    
    /**
     * 
     * Adds a new key-value pair.
     * 
     * @param key
     * @param type
     * @param value
     * @return 
     */
    template<typename T> bool add_property(std::string key, int type, T* value) {
        TypedPair pair = std::make_pair<int, void*>(type, value);
        std::pair < FBStatsMap::iterator, bool> insertion_status =
                m_map->insert(std::make_pair<std::string, TypedPair>(key, pair));
        return insertion_status.second;
    }    

    /**
     * 
     * @param key
     * @param type
     * @param value
     * @return 
     */
    bool get_typed_property(std::string key, int& type, void*& value);

    /**
     * 
     * @param key
     * @param value
     * @return 
     */
    template<typename T> bool get_property(std::string key, T*& value) {
        FBStatsMap::iterator value_it = m_map->find(key);
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
    
    size_t size();;

private:

    void init();

    FBStatsMap * m_map;

};

#endif	/* FBSTATS_H */


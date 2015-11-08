/* 
 * File:   ForBESUtils.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 24, 2015, 5:17 PM
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

#ifndef FORBESUTILS_H
#define	FORBESUTILS_H

/**
 * \brief ForBES utilities.
 */
class ForBESUtils {
public:
    
     /**
     * Method has succeeded.
     */
    const static int STATUS_OK;
    /**
     * Method is undefined.
     */
    const static int STATUS_UNDEFINED_FUNCTION;
    /**
     * The result is unreliable, or could not be computed because
     * of numerical errors.
     */
    const static int STATUS_NUMERICAL_PROBLEMS;
    
    const static int STATUS_HAD_TO_REALLOC;
    
    
private:
    
    ForBESUtils();
    virtual ~ForBESUtils();
    
};

#endif	/* FORBESUTILS_H */


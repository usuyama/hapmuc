/*    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * OutputData.hpp
 *
 *  Created on: Sep 7, 2009
 *      Author: caa
 */

#ifndef STRINGHASH_HPP_
#define STRINGHASH_HPP_
#include <vector>
#include <ext/hash_map>
#include <string>
#include <iostream>
#include "foreach.hpp"
using namespace std;
using namespace boost;
using __gnu_cxx::hash;
  
struct my_hash_funct : public unary_function<string, size_t>
{	
	size_t operator()(const string & x) const { return hash<const char*>() (x.c_str()); }
};

template<class T>  class string_hash : public hash_map<string, T, my_hash_funct> {};

#endif /* OUTPUTDATA_HPP_ */

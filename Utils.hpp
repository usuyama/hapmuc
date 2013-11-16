/*
 * Dindel
 * http://www.sanger.ac.uk/resources/software/dindel/
 *
 * Copyright 2010, Kees Albers
 * Released under the GPL Version 3.
 */
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
 * Utils.hpp
 *
 *  Created on: Mar 11, 2009
 *      Author: caa
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
using namespace std;

inline double addLogs(const double l1, const double l2)
{
	if (l1>l2) {
		double diff=l2-l1;
		return l1+log(1.0+exp(diff));
	} else {
		double diff=l1-l2;
		return l2+log(1.0+exp(diff));
	}
}

double lbeta(double a, double b);

template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}



#endif /* UTILS_HPP_ */

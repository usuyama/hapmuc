//
//  Utils.cpp
//  hapmuc2
//
//  Created by 直人 臼山 on 10/12/12.
//  Copyright 2012 Univ. of Tokyo. All rights reserved.
//

#include <iostream>
#include "Utils.hpp"
#include <boost/math/special_functions/digamma.hpp>
double lbeta(double a, double b) {
    return lgamma(a) + lgamma(b) - lgamma(a+b);
}


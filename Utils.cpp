/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#include <iostream>
#include "Utils.hpp"
#include <boost/math/special_functions/digamma.hpp>
double lbeta(double a, double b) {
    return lgamma(a) + lgamma(b) - lgamma(a+b);
}


/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#include <iostream>
#include "Utils.hpp"
double lbeta(double a, double b) {
    return lgamma(a) + lgamma(b) - lgamma(a+b);
}


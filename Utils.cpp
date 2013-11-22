/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */
#include <iostream>
#include "Utils.hpp"
#include <string>
#include <sstream>
#include <vector>
double lbeta(double a, double b) {
  return lgamma(a) + lgamma(b) - lgamma(a+b);
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

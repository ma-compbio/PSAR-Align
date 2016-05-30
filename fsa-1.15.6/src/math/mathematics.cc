
/**
 * \file mathematics.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "math/mathematics.h"

using namespace fsa;

const double Mathematics::double_tiny = 0.001;
const double Mathematics::double_very_tiny = 0.000001;

unsigned Mathematics::factorial (unsigned x) {

  unsigned f = 1;

  while (x > 0) {
    f *= x--;
  }

  return f;
}

double Mathematics::log_factorial (unsigned x) {

  double logf = 0;

  while (x > 1) {
    logf += std::log (static_cast<double> (x--));
  }

  return logf;
}


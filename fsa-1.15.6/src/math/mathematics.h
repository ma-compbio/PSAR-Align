
/**
 * \file mathematics.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef MATH_MATHEMATICS_INCLUDED
#define MATH_MATHEMATICS_INCLUDED

#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>

#include "config.h"
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#include "util/misc.h"

namespace fsa {

  struct Mathematics {

  public:

    /**
     * \brief Is a number a power of 2?
     */
    static bool is_power_of_2 (const uintmax_t num);

    /**
     * \brief Calculates power of 2.
     */
    static uintmax_t power_of_2 (const unsigned n);

    /**
     * \brief Calculates power of k^n, where both k and n are integers.
     */
    static uintmax_t power_of_integer (const unsigned k, const unsigned n);

    /**
     * \brief Calculate factorial x!.
     */
    static unsigned factorial (unsigned x);

    /**
     * \brief Calculate log-factorial log (x!).
     */
    static double log_factorial (unsigned x);

    /**
     * \brief Compute the average value.
     */
    template<typename T>
    static double mean (const std::vector<T>& values);

    /**
     * \brief Compute the median value.
     * \param is_sorted true if the data are already sorted low-to-high
     */
    template<typename T>
    static T median (std::vector<T>& values,
		     const bool is_sorted = false);

    /**
     * \brief Compute the value at a particular percentile (in [0, 1]).
     * \param is_sorted true if the data are already sorted low-to-high
     */
    template<typename T>
    static T percentile_value (std::vector<T>& values, const float percentile,
			       const bool is_sorted = false);

    /**
     * \brief Return original value if inside range, or bound if outside.
     */
    template<typename N>
    static N bounded_value (N n, N min_bound, N max_bound) {
      return n < min_bound
	? min_bound
	: (n > max_bound ? max_bound : n);
    }

    /**
     * \brief Function object to compute the normalized difference function (x - y) / (x + y).
     */
    template<typename T>
    struct normalized_difference {
      double operator() (const T x, const T y) const {
	return (x - y) / (x + y);
      }
    };



    static const double double_tiny;             ///< small double value
    static const double double_very_tiny;        ///< very small double value

  };

  inline bool Mathematics::is_power_of_2 (const uintmax_t num) {
    return (num & (num - 1)) == 0;
  }

  inline uintmax_t Mathematics::power_of_2 (const unsigned n) {
#ifndef NDEBUG
    if (std::pow (static_cast<double> (2), static_cast<int> (n)) > std::numeric_limits<uintmax_t>::max()) {
      cerr << "ERROR: Overflow of uintmax_t type: attempting to calculate 2^" << n << endl;
      exit (1);
    }      
#endif
    return static_cast<uintmax_t> (1) << n;
  }

  inline uintmax_t Mathematics::power_of_integer (const unsigned k, const unsigned n) {
#ifndef NDEBUG
    if (std::pow (static_cast<double> (k), static_cast<int> (n)) > std::numeric_limits<uintmax_t>::max()) {
      cerr << "ERROR: Overflow of uintmax_t type: attempting to calculate " << k << "^" << n << endl;
      exit (1);
    }      
#endif
    uintmax_t sum = 1;
    for (size_t i = 0; i < n; ++i)
      sum *= k;
    return sum;
  }

  template<typename T>
    double Mathematics::mean (const std::vector<T>& values) {

    if (!values.size()) {
      cerr << "ERROR: Tried to take the mean of an empty vector." << endl;
      exit (1);
    }

    double total = 0.0;
    for (typename std::vector<T>::const_iterator value = values.begin(); value != values.end(); ++value)
      total += *value;

    return (total / values.size());

  }
  
  template<typename T>
    T Mathematics::median (std::vector<T>& values,
			   const bool is_sorted /* = false */) {

    return percentile_value (values, 0.5,
			     is_sorted);

  }

  template<typename T>
    T Mathematics::percentile_value (std::vector<T>& values, const float percentile,
				     const bool is_sorted /* = false */) {

    // check sane
    if (percentile < 0.0 || percentile > 1.0) {
      cerr << "ERROR: A percentile of " << percentile << " is not meaningful." << endl;
      exit (1);
    }

    if (!values.size()) {
      cerr << "ERROR: Tried to find a percentile within an empty vector." << endl;
      exit (1);
    }

    // sort values if necessary
    if (!is_sorted)
      std::sort (values.begin(), values.end());

    // pull out the index into the vector for the appropriate percentile
    size_t index = static_cast<size_t> (std::floor (values.size() * percentile));

    // catch the case of percentile_value = 1
    if (index == values.size())
      index = values.size() - 1;
    else if (index > values.size()) {
      cerr << "ERROR: Array index indicated by percentile value overflows the vector indices!" << endl
	   << "index = " << index << "; array size = " << values.size() << endl;
      exit (1);
    }

    // take array slice to get appropriate percentile
    return values[index];

  }

}

#endif /* MATH_MATHEMATICS_INCLUDED */

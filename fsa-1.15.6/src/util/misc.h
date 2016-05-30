
/**
 * \file misc.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef UTIL_MISC_INCLUDED
#define UTIL_MISC_INCLUDED

#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <map>

#include "config.h"
#include "util/regexp.h"

using std::cerr;
using std::cout;
using std::endl;

namespace fsa {

  struct Util {

  public:

    /**
     * \brief Convert to a string.
     */
    template<typename T>
    static std::string to_string (const T& t);

    /**
     * \brief Remove newline character, if it exists, from the end of a string.
     */
    static void chomp (std::string& str);

    /**
     * \brief Convert a string to upper-case.
     */
    static void toupper (std::string& str);

    /**
     * \brief Convert a string to lower-case.
     */
    static void tolower (std::string& str);

    /**
     * \brief Tokenize a string according to the specified delimiters.
     *
     * Code from:
     * http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
     */
    static std::vector<std::string> split (const std::string& str,
					   const std::string& delimiters = " ");

    /**
     * \brief Join strings into a single string separated by delimiter.
     */
    template<typename T>
    static std::string join (const std::vector<T>& items,
			     const std::string& delimiter = "");

    /**
     * \brief Search for an instance of a specified item in a vector of items.
     * \param item list of items to search through
     * \param query item to search for
     */
    template<typename T>
    static bool contains (const std::vector<T>& items,
			  const T query);

    /**
     * \brief Add prefix to all strings in vector.
     */
    static std::vector<std::string> addprefix (const std::vector<std::string>& items,
					       const std::string prefix);

    /**
     * \brief Add suffix to all strings in vector.
     */
    static std::vector<std::string> addsuffix (const std::vector<std::string>& items,
					       const std::string suffix);

    /**
     * \brief Get keys in a map or multimap.
     *
     * Returns possibly-duplicate key values if T is a multimap.
     * \param container map or multimap
     */
    template<typename T>
    static std::vector<typename T::key_type> keys (const T& container);

    /**
     * \brief Get values in a map or multimap.
     * \param container map or multimap
     */
    template<typename T>
    static std::vector<typename T::mapped_type> values (const T& container);

    /**
     * \brief Strip leading 'chr' (if present) from the passed chromosome name.
     */
    static void strip_leading_chr (std::string& chromosome);

    /**
     * \brief Return a random number in [0, max].
     */
    static unsigned rand (const unsigned max);

    /**
     * \brief Return a random probability (double in [0, 1]).
     */
    static double prob();

    /**
     * \brief Seed random number generator.
     */
    static void seed_rand (const unsigned seed);

    /**
     * \brief Seed random number generate on current time.
     */
    static void seed_rand_on_time();

    /**
     * \brief Choose from a probability vector.
     * 
     * Assumes that the distribution is properly normalized.
     */
    static size_t choose_from_distribution (const std::vector<double>& dist);

    /**
     * \brief Does a file exist?
     *
     * Uses stat to see if we can get the file attributes.
     * From http://www.techbytes.ca/techbyte103.html.
     */
    static bool exists_file (const std::string& filename);

    /**
     * \brief Counts the number of newlines in a file.
     * \param filename name of file
     */
    static size_t count_newlines (const std::string& filename);

    /**
     * \brief Remove directory information from filename.
     *
     * Converts, e.g., '/tmp/rob.data' to 'rob.data' or 'rob'.
     * \param strip_extension remove extension from filename as well
     */
    static void basename (std::string& filename,
			  const bool strip_extension = false);

    /**
     * \brief Case-insensitive character comparison.
     * \see String_equal_ci
     */
    static bool char_less_ci (const char c1, const char c2) {
      return std::tolower (static_cast<unsigned char> (c1)) < std::tolower (static_cast<unsigned char> (c2));
    }

    /**
     * \brief Case-insensitive string comparison.
     * 
     * Taken from Meyers, "Effective C++" 3rd ed.
     */
    struct String_equal_ci {
    public:
      bool operator() (const std::string& s1, const std::string& s2) const;
    };

    /**
     * \brief Function object for comparing duples.
     *
     * Orders by first coordinate, then second coordinate.
     */
    template<typename T1, typename T2>
    struct Duple_less : std::binary_function<T1, T2, bool> {
    public:
      bool operator() (const std::pair<T1, T2> l, const std::pair<T1, T2> r) const {
	if (l.first == r.first)
	  return (l.second < r.second);
	else
	  return (l.first < r.first);
      }
    };

    /**
     * \brief Function object for comparing map entries based on their values.
     *
     * T might be, e.g., std::map<char, size_t>
     */
    template<typename T>
    struct Map_value_less : std::binary_function<T, T, bool> {
    public:
      bool operator() (const typename T::value_type& l, const typename T::value_type& r) const {
	return l.second < r.second;
      }
    };


  private:

    static Regexp re_chr_stripper;                  ///< \see strip_leading_chr
    static Regexp re_basename;                      ///< \see basename
    static Regexp re_basename_extension;            ///< \see basename
    
  };


  template<typename T>
    std::string Util::to_string (const T& t) {

    std::stringstream ss;
    ss << t;

    return ss.str();

  }

  template<typename T>
    std::string Util::join (const std::vector<T>& items,
			    const std::string& delimiter /* = "" */) {
    
    std::stringstream joined;

    // case the catch of no items
    if (!items.size())
      return joined.str();

    // awkward loop structure is to avoid adding the delimiter onto the end of the string
    typename std::vector<T>::const_iterator item;
    for (item = items.begin(); (item + 1) != items.end(); ++item)
      joined << *item << delimiter;
    joined << *item;

    return joined.str();

  }

  inline std::vector<std::string> Util::addprefix (const std::vector<std::string>& items,
						   const std::string prefix) {

    std::vector<std::string> modified;
    if (!items.size())
      return modified;

    for (std::vector<std::string>::const_iterator item = items.begin(); item != items.end(); ++item)
      modified.push_back (prefix + *item);

    return modified;

  }

  inline std::vector<std::string> Util::addsuffix (const std::vector<std::string>& items,
						   const std::string suffix) {

    std::vector<std::string> modified;
    if (!items.size())
      return modified;

    for (std::vector<std::string>::const_iterator item = items.begin(); item != items.end(); ++item)
      modified.push_back (*item + suffix);

    return modified;

  }

  template<typename T>
    bool Util::contains (const std::vector<T>& items,
			 const T query) {

    for (typename std::vector<T>::const_iterator item = items.begin(); item != items.end(); ++item) {
      if (*item == query)
	return true;
    }

    return false;

  }

  template<typename T>
    std::vector<typename T::key_type> Util::keys (const T& container) {

    std::vector<typename T::key_type> keys;
    for (typename T::const_iterator item = container.begin(); item != container.end(); ++item)
      keys.push_back (item->first);

    return keys;

  }

  template<typename T>
    std::vector<typename T::mapped_type> Util::values (const T& container) {

    std::vector<typename T::mapped_type> values;
    for (typename T::const_iterator item = container.begin(); item != container.end(); ++item)
      values.push_back (item->second);

    return values;

  }

  inline void Util::chomp (std::string& str) {

    const int end = str.length() - 1;
    if (end >= 0 && str[end] == '\n')
      str.erase (end);

  }

  inline void Util::toupper (std::string& str) {

    for (std::string::iterator s = str.begin(); s != str.end(); ++s)
      *s = std::toupper (*s);

  }

  inline void Util::tolower (std::string& str) {

    for (std::string::iterator s = str.begin(); s != str.end(); ++s)
      *s = std::tolower (*s);

  }

  inline std::vector<std::string> Util::split (const std::string& str,
					       const std::string& delimiters /* = " " */) {

    std::vector<std::string> tokens;

    // Skip delimiters at beginning.
    std::string::size_type lastpos = str.find_first_not_of (delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of (delimiters, lastpos);

    while (std::string::npos != pos || std::string::npos != lastpos) {
      // Found a token, add it to the vector.
      tokens.push_back (str.substr (lastpos, pos - lastpos));
      // Skip delimiters.  Note the "not_of"
      lastpos = str.find_first_not_of (delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of (delimiters, lastpos);
    }

    return tokens;

  }

  inline unsigned Util::rand (const unsigned max) {
    assert (max < RAND_MAX);
    return std::rand() % (max + 1);
  }

  inline double Util::prob() {
    return static_cast<double> (std::rand()) / RAND_MAX;
  }

  inline void Util::seed_rand (const unsigned seed) {
    std::srand (seed);
  }

  inline void Util::seed_rand_on_time() {
    std::srand (std::time (NULL));
  }

  inline size_t Util::choose_from_distribution (const std::vector<double>& dist) {
    double prob = Util::prob();
    for (size_t i = 0; i < dist.size() - 1; ++i) {
      prob -= dist[i];
      if (prob <= 0.)
	return i;
    }
    
    // handle last case separately to accomodate floating-point error
    return dist.size() - 1;

  }

}

#endif /* UTIL_MISC_INCLUDED */

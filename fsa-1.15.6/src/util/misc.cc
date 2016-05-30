
/**
 * \file misc.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "util/misc.h"
#include <algorithm>
#include <fstream>

using namespace fsa;


Regexp Util::re_chr_stripper = Regexp ("^(chr)?(.+)");
Regexp Util::re_basename = Regexp ("^.*\\/([^\\/]+)$");
Regexp Util::re_basename_extension = Regexp ("^(.+)\\.[^\\.]+$");

bool Util::exists_file (const std::string& filename) {

  struct stat fileinfo;
  bool gotattr;

  // Attempt to get the file attributes
  int statresults = stat (filename.c_str(), &fileinfo);

  if (statresults == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    gotattr = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    gotattr = false;

  }
  
  return gotattr; 

}

size_t Util::count_newlines (const std::string& filename) {

  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << filename << "' for reading." << endl;
    exit (1);
  }

  size_t num = std::count (std::istreambuf_iterator<char> (filestream), std::istreambuf_iterator<char>(),
			   '\n');

  return num;

}

void Util::strip_leading_chr (std::string& chromosome) {

  if (re_chr_stripper.Match (chromosome.c_str())) // strip off leading 'chr' if present
    chromosome = std::string (re_chr_stripper[2]);

}

void Util::basename (std::string& filename,
		     const bool strip_extension /* = false */) {

  // remove directory information
  if (re_basename.Match (filename.c_str()))
    filename = std::string (re_basename[1]);

  // remove filename extension
  if (strip_extension && re_basename_extension.Match (filename.c_str()))
    filename = re_basename_extension[1];

}

bool Util::String_equal_ci::operator() (const std::string& s1, const std::string& s2) const {

    return std::lexicographical_compare (s1.begin(),
					 s1.end(),
					 s2.begin(),
					 s2.end(),
					 Util::char_less_ci);

}


/**
 * \file gff.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "util/regexp.h"
#include "math/mathematics.h"
#include "seq/gff.h"

using namespace fsa;

const std::string GFF::comment_char = "#";
const std::string GFF::undef_char = ".";

const char GFF::unknown_strand = '.';

const std::string GFF::attributes_split_char = ";";
const std::string GFF::attributes_assign_char = "=";
const std::string GFF::attributes_list_char = ",";
    
const std::string GFF::key_id = "ID";
const std::string GFF::key_name = "Name";

const size_t GFF::default_flanking = 10000;

Regexp GFF::re_key_value = "^([^=]+)=(.*)$";

void GFF::from_string (const std::string& str, const bool strip_leading_chr /* = true */) {

  std::vector<std::string> tokens = Util::split (str, "\t"); // hold tab-separated tokens

  // now parse tokens

  // sanity checks
  if (tokens.size() != 9) {
    cerr << "ERROR: Not a GFF-format line: " << str << endl;
    exit (1);
  }

  seqid = tokens[0];
  if (strip_leading_chr)
    Util::strip_leading_chr (seqid);
  source = tokens[1];
  type = tokens[2];
  set_start (atoi (tokens[3].c_str()));
  end = static_cast<unsigned> (atoi (tokens[4].c_str()));
  score = (tokens[5] != "" && tokens[5] != GFF::undef_char) ? atof (tokens[5].c_str()) : -1.;
  strand = (tokens[6])[0];
  phase = (tokens[7] != "" && tokens[7] != GFF::undef_char) ? static_cast<unsigned> (atoi (tokens[7].c_str())) : 3;

  parse_attributes_string (tokens[8]);
  
}

void GFF::set_start (const int s) {

  if (s < 1) {
    cerr << "WARNING: setting GFF feature start coordinate " << s << " to 1." << endl;
    start = 1;
  }
  
  else {
    start = static_cast<unsigned> (s);
  }

}

void GFF::parse_attributes_string (const std::string& str) {

  std::vector<std::string> key_value_pairs = Util::split (str, GFF::attributes_split_char);
  for (std::vector<std::string>::const_iterator key_value = key_value_pairs.begin(); key_value != key_value_pairs.end(); ++key_value) {
    if (re_key_value.Match (key_value->c_str())) {
      std::string key = re_key_value[1];
      std::vector<std::string> values = Util::split (re_key_value[2], GFF::attributes_list_char);
      if (values.size())
	attributes_ordering.push_back (key);
      for (std::vector<std::string>::const_iterator value = values.begin(); value != values.end(); ++value)
	attributes_map[key].push_back (*value);
    }
  }

}

std::string GFF::get_attributes_string() const {

  std::string s;

  // iterate over keys
  for (std::vector<std::string>::const_iterator key = attributes_ordering.begin(); key != attributes_ordering.end(); ++key) {
    std::map<std::string, std::vector<std::string> >::const_iterator key_value = attributes_map.find (*key);
    // iterate over values
    if (key_value->second.size())
      s += *key + GFF::attributes_assign_char;
    for (std::vector<std::string>::const_iterator value = key_value->second.begin(); value != key_value->second.end(); ++value)
      s += *value + GFF::attributes_list_char;
    if (s.length())
      s.erase (s.end() - 1, s.end()); // strip off final comma
    s += GFF::attributes_split_char;
  }

  if (!s.length())
    s += GFF::undef_char;

  return s;

}

std::string GFF::to_string() const {

  std::stringstream ss;

  ss << (seqid != "" ? seqid : GFF::undef_char) << "\t"
     << (source != "" ? source : GFF::undef_char) << "\t"
     << (type != "" ? type : GFF::undef_char) << "\t"
     << start << "\t" << end << "\t";
  if (score != -1.) { ss << score; }
  else { ss << GFF::undef_char; }
  ss << "\t";
  ss << strand << "\t";
  if (phase != 3) { ss << phase; }
  else { ss << GFF::undef_char; }
  ss << "\t";
  ss << get_attributes_string();

  return ss.str();

}

void GFF_database::from_file (const std::string& filename, const bool strip_leading_chr /* = true */) {

  Regexp re_gff ("^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*$");

  // open file
  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << filename << "' for reading." << endl;
    exit (1);
  }

  cerr << "Reading GFF file '" << filename << "'...";

  // parse file
  std::string line;
  while (!filestream.eof()) {
    getline (filestream, line);
    Util::chomp (line);
    // skip lines which don't seem to be in GFF format
    if (!line.length() || !re_gff.Match (line.c_str()))
      continue;
    GFF gff;
    gff.from_string (line, strip_leading_chr);
    store_entry (gff);
    __maxlen = std::max (__maxlen, gff.length());
  }

  cerr << "done." << endl;

  // clean up
  filestream.close();

}

void GFF_database::write (std::ostream& o) const {

  for (std::vector<GFF>::const_iterator gff = this->begin(); gff != this->end(); ++gff)
    o << *gff;

}

void GFF_database::append (const GFF_database& gff_db) {

  for (std::vector<GFF>::const_iterator gff = gff_db.begin(); gff != gff_db.end(); ++gff)
    this->store_entry (*gff);
  __maxlen = (gff_db.maxlen() > __maxlen) ? gff_db.maxlen() : __maxlen;

}

void GFF_database::create_unique_ids() {

  if (!size())
    return;

  size_t width = static_cast<size_t> (log10 (size()));
  std::string format = "%" + Util::to_string (width) + "d";
  for (size_t i = 0; i < size(); ++i) {
    std::string id (width - static_cast<size_t> (log10 (i)), '0');  // hack to format the IDs nicely
    id += Util::to_string (i);
    __entries[i].set_id (id);
  }

}

void GFF_database::sort_entries() {

  std::sort (__entries.begin(), __entries.end(), GFF::GFF_less());
  __is_sorted = true;

}

GFF GFF_database::find_closest_feature_five_prime (const std::string& chromosome,
						   unsigned start, unsigned end,
						   const size_t flanking /* = default_flanking */) const {

  return find_closest_feature (chromosome,
			       start, end,
			       flanking,
			       true); // use_feature_five_prime

}

GFF GFF_database::find_closest_feature_three_prime (const std::string& chromosome,
						    unsigned start, unsigned end,
						    const size_t flanking /* = default_flanking */) const {

  return find_closest_feature (chromosome,
			       start, end,
			       flanking,
			       false); // use_feature_five_prime

}

GFF GFF_database::find_closest_feature (const std::string& chromosome,
					unsigned start, unsigned end,
					const size_t flanking /* = default_flanking */,
					const bool use_feature_five_prime /* = true */) const {

  const unsigned centroid = start + (end - start + 1) / 2;
  start = (start < flanking) ? 0 : start - flanking;
  end += flanking;

  const GFF_database isects = intersect_genomic_interval (chromosome,
							  start, end);

  if (!isects.size())
    return GFF();

  size_t closest_i = 0;
  size_t closest_distance = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < isects.size(); ++i) {
    const GFF& gff = isects[i];
    size_t distance = 0;
    unsigned gff_boundary; // either the 5' or 3' end of the feature, as specified by use_feature_five_prime (0-based)
    if (gff.strand == '-')
      gff_boundary = use_feature_five_prime ? gff.end - 1 : gff.start - 1;
    else
      gff_boundary = use_feature_five_prime ? gff.start - 1 : gff.end - 1;
    distance = std::max (centroid, gff_boundary) - std::min (centroid, gff_boundary);
    if (distance < closest_distance) {
      closest_i = i;
      closest_distance = distance;
    }
  }

  return isects[closest_i];

}

GFF_database GFF_database::intersect_genomic_interval (const std::string& chromosome,
						       const unsigned start, const unsigned end) const {

  if (!__is_sorted) {
    cerr << "ERROR: You must call sort_entries before attepting use this function." << endl;
    exit (1);
  }

  GFF_database intersections;

  // check for empty interval or no features
  if (end < start || !size())
    return intersections;

  // initialize dummy object for comparison
  GFF g;
  g.seqid = chromosome;
  Util::strip_leading_chr (g.seqid);

  // two intervals 1 and 2 intersect iff
  // 1.start <= 2.end && 1.end >= 2.start
  // features in database are 1; query interval is 2

  // check first condition (1.start <= 2.end)
  g.start = end + 1; // convert to 1-based coordinates
  g.end = end + 1;
  std::vector<GFF>::const_iterator gff_upper = std::upper_bound (__entries.begin(),
								 __entries.end(),
								 g,
								 GFF::GFF_less());

  // decrement element to try to get to the upper-bound entry
  // which (possibly) intersects the query interval
  if (gff_upper != __entries.begin())
    --gff_upper;

  // now check that we're on the correct chromosome
  if (gff_upper->seqid != chromosome)
    return intersections;

  // check that gff_upper doesn't fall completely to the right of the query interval
  // if gff_upper is __entries.begin(), then it may!
  if (gff_upper->start > end)
    return intersections;

  // now iterate backwards through the sorted entries,
  // checking the second condition (1.end >= 2.start)
  // as we go
  while (gff_upper != __entries.begin() - 1) {

    // if we're no longer on the correct chromosome, then we're done
    if (gff_upper->seqid != chromosome)
      break;

    // does it intersect the query interval?
    if (gff_upper->end >= start)
      intersections.store_entry (*gff_upper);

    // check whether we're done
    // (i.e., whether given __maxlen, there cannot be any more features close
    // enough to the query interval to intersect it)
    if (gff_upper->start + __maxlen < start)
      break;

    --gff_upper;
  }

  // sort nicely
  intersections.sort_entries();

  return intersections;

}

size_t GFF_database::meanlen() const {

  std::vector<size_t> lengths (size(), 0);
  for (size_t i = 0; i < size(); ++i)
    lengths[i] = (*this)[i].length();

  assert (static_cast<size_t> (Mathematics::mean (lengths)) <= maxlen());

  return static_cast<size_t> (Mathematics::mean (lengths));

}


size_t GFF_database::medianlen() const {

  std::vector<size_t> lengths (size(), 0);
  for (size_t i = 0; i < size(); ++i)
    lengths[i] = (*this)[i].length();
  std::sort (lengths.begin(), lengths.end());
  
  assert (Mathematics::median (lengths) <= maxlen());

  return Mathematics::median (lengths);

}


/**
 * \file constraints.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "util/logfile.h"
#include "fsa/constraints.h"

using namespace fsa;

void Constraints::write_mercator (std::ostream& o) const {

  for (std::vector<Constraint>::const_iterator constraint = begin(); constraint != end(); ++constraint)
    constraint->write_mercator (__xname, __yname, o);

}


Constraints_set::Constraints_set (const Sequence_database& seq_db)
  : seq_db (seq_db)
{ }

void Constraints_set::read_mercator (const std::string& filename) {

  std::ifstream filestream (filename.c_str(), std::ifstream::in);
  if (!filestream.is_open())
    THROWEXPR ("ERROR: Couldn't open file '" << filename << "'.");

  std::string line;
  while (!filestream.eof()) {
    getline (filestream, line);
    // log
    if (CTAGGING(-1,CONSTRAINTS_VERBOSE))
      CL << line;
    std::string buffer;
    // skip comments and sequence names
    std::stringstream ss (line);
    std::vector<std::string> tokens; // vector to hold whitespace-separated tokens (words)
    while (ss >> buffer)
      tokens.push_back (buffer);

    // now parse tokens into constraints
    // sample line:
    //   drosophila_melanogaster-5.0 30 1344 drosophila_yakuba-2.0 413 1727
    // Note that Mercator constraints are 0-based half-open [start, end).
    if (tokens.size() == 6) {
      std::string xname = tokens[0];
      unsigned xstart = atoi (tokens[1].c_str()); // [start, end)
      unsigned xend = atoi (tokens[2].c_str()) - 1;
      std::string yname = tokens[3];
      unsigned ystart = atoi (tokens[4].c_str());
      unsigned yend = atoi (tokens[5].c_str()) - 1;

      // sanity checks on sequence names
      if (!seq_db.exists_seq (xname) || !seq_db.exists_seq (yname))
	THROWEXPR ("ERROR: I don't recognize the sequence names '" << xname << "' and '" << yname << "' in constraint file '" << filename << "'.");

      // map sequence names to indices
      size_t i = seq_db.get_seq_index (xname);
      size_t j = seq_db.get_seq_index (yname);

      // be consistent
      if (i > j) {
	std::swap (i, j);
	std::swap (xname, yname);
	std::swap (xstart, ystart);
	std::swap (xend, yend);
      }

      // if no constraint information for this sequence pair,
      // create a new Constraints object
      if (constraints_list.find (std::make_pair (i, j)) == constraints_list.end())
	constraints_list.insert (std::make_pair (std::make_pair (i, j), Constraints (xname, yname)));

      // store constraint
      (constraints_list.find (std::make_pair (i, j)))->second.store (Constraint (Interval (xstart, xend), Interval (ystart, yend)));

    } else if (tokens.size() == 0) {
      continue;
    } else {
      if (CTAGGING(4,CONSTRAINTS))
	CL << "WARNING: I don't know how to parse Mercator output line " << line << "; skipping it." << endl;
      continue;
    }
  }

  filestream.close();

  // log
  CTAG(8,CONSTRAINTS) << "Read Mercator constraints file '" << filename << "'." << endl;

}

void Constraints_set::write_mercator (std::ostream& o) const {

  for (std::map<std::pair<size_t, size_t>, Constraints>::const_iterator iter = constraints_list.begin(); iter != constraints_list.end(); ++iter)
    iter->second.write_mercator (o);

}

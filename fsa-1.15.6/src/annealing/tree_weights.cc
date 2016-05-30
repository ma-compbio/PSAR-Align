
/**
 * \file tree_weights.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "util/logfile.h"
#include "annealing/tree_weights.h"

using namespace fsa;

const double Tree_weights::min_sequence_weight = 0.000001;

void Tree_weights::from_file (const Sequence_database& seq_db, const std::string& filename) {

  // mark populated
  __populated = true;

  // initialize all weights to 0.0
  __weights.assign (seq_db.size(), std::vector<float> (seq_db.size(), 0.0));

  // read in weights from file
  std::ifstream filestream (filename.c_str(), std::ifstream::in);
  if (!filestream.is_open())
    THROWEXPR ("ERROR: Couldn't open file '" << filename << "'.");

  std::string line;
  while (!filestream.eof()) {
    getline (filestream, line);
    Util::chomp (line);
    if (!line.length())
      continue;

    // log
    if (CTAGGING(-1,TREE_WEIGHTS))
      CL << line;
    std::string buffer;
    // skip comments and sequence names
    std::stringstream ss (line);
    std::vector<std::string> tokens; // vector to hold whitespace-separated tokens (words)
    while (ss >> buffer)
      tokens.push_back (buffer);

    // skip lines which we don't understand
    if (tokens.size() != 3) {
      CTAG(8,ANNEALING) << "WARNING: I'm skipping unparseable line '" << line << "'" << endl;
      continue;
    }
    else if (tokens[0] == tokens[1])
      THROWEXPR ("ERROR: A sequence can't have a weight with itself!" << endl << line << endl);

    // format is:
    //   seqX seqY weight
    // symmetrize as we go
    if (!seq_db.exists_seq (tokens[0]) || !seq_db.exists_seq (tokens[1]))
      THROWEXPR ("ERROR: No sequence pair '" << tokens[0] << "' and '" << tokens[1] << "' in sequence file.");
    const size_t i = seq_db.get_seq_index (tokens[0]);
    const size_t j = seq_db.get_seq_index (tokens[1]);
    __weights[i][j] = __weights[j][i] = static_cast<float> (atof (tokens[2].c_str()));

  }

  // normalize weights
  normalize();

  CTAG(9,ANNEALING) << "Read sequence pair weights from file '" << filename << "'." << endl;

}

void Tree_weights::normalize() {

  const size_t num_seqs = __weights.size();
  const size_t num_seq_pairs = num_seqs * (num_seqs - 1) / 2;

  std::vector<double> seqtotal (num_seqs, 0.0); // total for single sequences
  double total = 0.0;                           // total for all sequence pairs
  for (size_t i = 0; i < num_seqs; ++i) {
    for (size_t j = i + 1; j < num_seqs; ++j) {
      seqtotal[i] += __weights[i][j];
      seqtotal[j] += __weights[i][j];
      total += __weights[i][j];
    }
  }
  for (size_t i = 0; i < num_seqs; ++i) {
    if (seqtotal[i] < Tree_weights::min_sequence_weight)
      THROWEXPR ("ERROR: Weights for sequence " << i << " are too small; the sequence will be left unaligned.");
  }

  if (total < num_seqs * Tree_weights::min_sequence_weight)
    THROWEXPR ("ERROR: Sequence pair weights are too small for reliable use.");

  const double scaling = num_seq_pairs / total;
  for (size_t i = 0; i < num_seqs; ++i) {
    for (size_t j = i + 1; j < num_seqs; ++j)
      __weights[i][j] *= scaling;
  }
  
}


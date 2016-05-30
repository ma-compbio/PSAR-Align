
/**
 * \file model.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "fsa/model.h"

using namespace fsa;

Params::Params() {
  
  // initialize state names
  states.resize (NUM_HMM_STATES);
  states[0] = "Start"; states[1] = "Match"; states[2] = "Insert1"; states[3] = "Delete1"; states[4] = "Insert2"; states[5] = "Delete2"; states[6] = "End";

  // initialize transition matrix
  transition_matrix.assign (NUM_HMM_STATES, vector<double> (NUM_HMM_STATES, 0.));

  // initialize stationary distribution
  stat_dist.assign (NUM_HMM_STATES, 0.);

}

void Params::init_transition_matrix (bool left_match /* = false */, bool right_match /* = false */, bool ragged_ends /* = false */) {
    
  // these must manually be set to 0 because they are re-estimated
  // by train_params even if the HMM only has 1 set of indel states
  // (sanity checks against divide-by-zero errors give it a non-zero value)
  if (!is_indel2)
    gap_open2 = gap_extend2 = 0.;

  // calculate stationary distribution of chain
  stat_dist[1] = (1 - gap_extend1) * (1 -gap_extend2) / ((1 - gap_extend1) * (1 - gap_extend2) + (gap_open1 / 2) * (1 - gap_extend2) + (gap_open2 / 2) * (1 - gap_extend1));
  stat_dist[2] = stat_dist[3] = ((gap_open1 / 2) / (2 * (1 - gap_extend1))) * stat_dist[1];
  stat_dist[4] = stat_dist[5] = ((gap_open2 / 2) / (2 * (1 - gap_extend2))) * stat_dist[1];

  // Start
  transition_matrix[0][1] = left_match ? 1.0 : 1.0 - gap_open1 - gap_open2;             // Start -> Match
  transition_matrix[0][2] = transition_matrix[0][3] = left_match ? 0. : gap_open1 / 2;  //       -> Insert1 | Delete1
  transition_matrix[0][4] = transition_matrix[0][5] = left_match ? 0. : gap_open2 / 2;  //       -> Insert2 | Delete2

  // Match
  transition_matrix[1][1] = 1.0 - gap_open1 - gap_open2 - to_end;      // Match -> Match
  transition_matrix[1][2] = transition_matrix[1][3] = gap_open1 / 2;   //       -> Insert1 | Delete1
  transition_matrix[1][4] = transition_matrix[1][5] = gap_open2 / 2;   //       -> Insert2 | Delete2
  transition_matrix[1][6] = to_end;                                    //       -> End

  // Insert1 and Delete1
  transition_matrix[2][2] = transition_matrix[3][3] = gap_extend1;                                                    // Insert1 -> Insert1
  transition_matrix[2][1] = transition_matrix[3][1] = right_match ? 1.0 - gap_extend1 : 1.0 - gap_extend1 - to_end;   //         -> Match
  transition_matrix[2][6] = transition_matrix[3][6] = right_match ? 0. : to_end;                                      //         -> End

  // Insert2 and Delete2
  transition_matrix[4][4] = transition_matrix[5][5] = gap_extend2;                                                    // Insert2 -> Insert2
  transition_matrix[4][1] = transition_matrix[5][1] = right_match ? 1.0 - gap_extend2 : 1.0 - gap_extend2 - to_end;   //         -> Match
  transition_matrix[4][6] = transition_matrix[5][6] = right_match ? 0. : to_end;                                      //         -> End

  assert_transition_matrix_valid();

}

void Params::normalize_transition_matrix() {
    
  for (size_t i = 0; i < transition_matrix.size(); ++i) {
    if (i == 6) // don't try to normalize transitions from End
      continue;
    double total = 0;
    for (size_t j = 0; j < transition_matrix[i].size(); ++j) {
      assert (transition_matrix[i][j] >= 0.);
      total += transition_matrix[i][j];
    }
    for (size_t j = 0; j < transition_matrix[i].size(); ++j)
      transition_matrix[i][j] /= total;
  }

}

void Params::assert_transition_matrix_valid() {

#ifndef NDEBUG
  for (size_t i = 0; i < transition_matrix.size(); ++i)
    for (size_t j = 0; j < transition_matrix.size(); ++j)
      if (transition_matrix[i][j] < 0 || transition_matrix[i][j] > 1) {
	show_transition_matrix (CL);
	THROWEXPR ("ERROR: Invalid probabilities in transition matrix.");
      }
#endif

}

void Params::copy_all (const Params& from) {
  copy_indel (from);
  copy_subst (from);
  bandwidth = from.bandwidth;
}

void Params::copy_indel (const Params& from) {
  is_indel2 = from.is_indel2;
  gap_open1 = from.gap_open1;
  gap_open2 = from.gap_open2;
  gap_extend1 = from.gap_extend1;
  gap_extend2 = from.gap_extend2;
  to_end = from.to_end;
  transition_matrix.assign (from.transition_matrix.begin(), from.transition_matrix.end());
}

void Params::copy_subst (const Params& from) {
  time = from.time;
  single_dist.assign (from.single_dist.begin(), from.single_dist.end());
  pair_dist.assign (from.pair_dist.begin(), from.pair_dist.end());
  alphabet_string = from.alphabet_string;
}

void Params::show_transition_matrix (std::ostream& o) const {

  // show raw params
  o << "gap_open1 = " << gap_open1 << endl;
  o << "gap_extend1 = " << gap_extend1 << endl;
  o << "gap_open2 = " << gap_open2 << endl;
  o << "gap_extend2 = " << gap_extend2 << endl;
  o << "to_end = " << to_end << endl;
  o << endl;

  // now show transition matrix
  const size_t width = 15;
  // state label row
  o.width (width); o << "";
  for (size_t i = 0; i < states.size(); ++i) {
    if (!is_indel2 && ((i == 4) || (i == 5)))
      continue;
    o.width (width);
    o << states[i];
  }
  o << endl;
  // transition matrix
  for (size_t i = 0; i < states.size(); ++i) {
    if (!is_indel2 && ((i == 4) || (i == 5)))
      continue;
    o.width (width);
    o << states[i];
    for (size_t j = 0; j < states.size(); ++j) {
      if (!is_indel2 && ((j == 4) || (j == 5)))
	continue;
      o.width (width);
      o << std::setprecision (PRECISION_DEFAULT) << transition_matrix[i][j];
    }
    o << endl;
  }
    
}

void Params::show_emission (std::ostream& o) const {

  const size_t width = 10;

  // single distribution
  // character label row
  o.width (width); o << "";
  for (size_t i = 0; i < single_dist.size(); ++i) {
    o.width (width);
    o << alphabet_string[i];
  }
  o << endl;
  // actual entries
  o.width (width); o << "";
  for (size_t i = 0; i < single_dist.size(); ++i) {
    o.width (width);
    const double p = (single_dist[i] > DOUBLE_VERY_TINY) ? single_dist[i] : DOUBLE_VERY_TINY;
    o << std::setprecision (PRECISION_DEFAULT) << p;
  }
  o << endl; o << endl;

  // pair distribution
  // character label row
  o.width (width); o << "";
  for (size_t i = 0; i < single_dist.size(); ++i) {
    o.width (width);
    o << alphabet_string[i];
  }
  o << endl;
  // actual entries
  for (size_t i = 0; i < pair_dist.size(); ++i) {
    o.width (width);
    o << alphabet_string[i];
    for (size_t j = 0; j < pair_dist.size(); ++j) {
      o.width (width);
      const double p = (pair_dist[i][j] > DOUBLE_VERY_TINY) ? pair_dist[i][j] : DOUBLE_VERY_TINY;
      o << std::setprecision (PRECISION_DEFAULT) << p;
    }
    o << endl;
  }
    
}

void Params::show (std::ostream& o) const {
  show_transition_matrix (o);
  show_emission (o);
}

void Params::write_transition (const std::string& filename) const {
  std::ofstream file (filename.c_str());
  if (!file)
    THROWEXPR ("ERROR: Couldn't create file with name '" << filename << "'.");
  show_transition_matrix (file);
  file.close();
}

void Params::write_emission (const std::string& filename) const {
  std::ofstream file (filename.c_str());
  if (!file)
    THROWEXPR ("ERROR: Couldn't create file with name '" << filename << "'.");
  show_emission (file);
  file.close();
}

void Params::write_emission_dat (const std::string& filename) const {
  std::ofstream file (filename.c_str());
  if (!file)
    THROWEXPR ("ERROR: Couldn't create file with name '" << filename << "'.");
  // print single_dist
  for (size_t i = 0; i < single_dist.size(); ++i)
    file << i+1 << ' ' << single_dist[i] << endl; // viewParams.pl wants 1-based indexing
  // print pair_dist
  for (size_t i = 0; i < single_dist.size(); ++i)
    for (size_t j = 0; j < single_dist.size(); ++j)
      file << i+1 << ' ' << j+1 << ' ' << pair_dist[i][j] << endl;
    
  file.close();
}

bool Params::is_biological() const {
  for (size_t i = 0; i < pair_dist.size(); ++i) {
    const double match_odds = pair_dist[i][i] / (single_dist[i] * single_dist[i]);
    for (size_t j = i + 1; j < pair_dist.size(); ++j) {
      const double mismatch_odds = max (pair_dist[i][j] / (single_dist[i] * single_dist[j]), pair_dist[j][i] / (single_dist[i] * single_dist[j]));
      if (match_odds < mismatch_odds)
	return false;
    }
  }
  return true;
}

double Params::branch_length (const unsigned which) const {
  vector<double> freq (pair_dist.size(), 0.);
  vector<vector<double> > submat (pair_dist.begin(), pair_dist.end());    

  assert ((which == 0) || (which == 1));
  if (which == 0) {
    // assemble frequencies implied by pair_dist P(x,y):
    // sum over second index to get frequency of first index P(x)
    for (size_t i = 0; i < pair_dist.size(); ++i)
      for (size_t j = 0; j < pair_dist.size(); ++j)
	freq[i] += pair_dist[i][j];
    // convert to conditional substitution matrix:
    // divide each entry by frequency of first index => P(x,y) / P(x) = P(y|x)
    // giving a distribution over the /second/ index 
    // (yes, this is counter-intuitive, but it makes sense in terms of the Ps)
    for (size_t i = 0; i < pair_dist.size(); ++i)
      for (size_t j = 0; j < pair_dist.size(); ++j)
	submat[i][j] /= freq[i];
  } else {
    for (size_t i = 0; i < pair_dist.size(); ++i)
      for (size_t j = 0; j < pair_dist.size(); ++j)
	freq[i] += pair_dist[j][i];
    for (size_t i = 0; i < pair_dist.size(); ++i)
      for (size_t j = 0; j < pair_dist.size(); ++j)
	submat[i][j] /= freq[j];
  }

  // check that the probabilities sum appropriately, ie
  // that we actually do have a stochastic matrix
  if (which == 0) {
    for (size_t i = 0; i < submat.size(); ++i) {
      double sum = 0.;
      for (size_t j = 0; j < submat.size(); ++j)
	sum += submat[i][j];
      assert (std::abs (1.0 - sum) < DOUBLE_VERY_TINY);
    }
  } else {
    for (size_t i = 0; i < pair_dist.size(); ++i) {
      double sum = 0;
      for (size_t j = 0; j < pair_dist.size(); ++j)
	sum += submat[j][i];
      assert (std::abs (1.0 - sum) < DOUBLE_VERY_TINY);
    }
  }

  return logdet (submat);
}

double Params::logdet (const vector<vector<double> >& matrix) {
  const double det = determinant (matrix);
  return (det > 0.) ? (-std::log (det)) : -1;
}

double Params::determinant (const vector<vector<double> >& matrix) {
  double det = 0;
  const size_t dim = matrix.size();
  assert (dim >= 2);

  if (dim == 2)
    det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  else {
    for (size_t j1 = 0; j1 < dim; ++j1) {
      vector<vector<double> > minor_mat (dim-1, vector<double> (dim-1));
      for (size_t i = 1; i < dim; ++i) {
	size_t j2 = 0;
	for (size_t j = 0; j < dim; ++j) {
	  if (j == j1)
	    continue;
	  minor_mat[i-1][j2] = matrix[i][j];
	  ++j2;
	}
      }
      det += ((j1 % 2 == 0) ? 1 : -1) * matrix[0][j1] * determinant (minor_mat);
    }
  }

  return det;
}

void Params::assert_normalized() const {
  double single_total = 0;
  double pair_total = 0;
  for (size_t i = 0; i < pair_dist.size(); ++i) {
    assert (!(single_dist[i] < 0));
    single_total += single_dist[i];
    for (size_t j = 0; j < pair_dist.size(); ++j) {
      assert (!(pair_dist[i][j] < 0));
      pair_total += pair_dist[i][j];
    }
  }
  assert (std::abs (1. - single_total ) < DOUBLE_TINY);
  assert (std::abs (1. - pair_total ) < DOUBLE_TINY);
}

void Params::normalize() {
  double single_total = 0;
  double pair_total = 0;
  for (size_t i = 0; i < pair_dist.size(); ++i) {
    assert (!(single_dist[i] < 0));
    single_total += single_dist[i];
    for (size_t j = 0; j < pair_dist.size(); ++j) {
      assert (!(pair_dist[i][j] < 0));
      pair_total += pair_dist[i][j];
    }
  }
  for (size_t i = 0; i < pair_dist.size(); ++i) {
    single_dist[i] /= single_total;
    for (size_t j = 0; j < pair_dist.size(); ++j) {
      pair_dist[i][j] /= pair_total;
    }
  }
  assert_normalized();
}

vector<double> Model::count_char_freq (const Sequence_database& seq_db, const std::string& alphabet_string, bool is_dna) {

  if (CTAGGING(5,FSA))
    CL << "Counting empirical character frequencies." << endl;

  map<char, size_t> char_cnt;
  for (size_t c = 0; c < alphabet_string.length(); ++c)
    char_cnt[alphabet_string[c]] = 0;

  double total = DOUBLE_TINY * alphabet_string.length(); // account for tiny values (ensure proper normalization)
  for (size_t i = 0; i < seq_db.size(); ++i) {
    const std::string& seq = seq_db.get_seq (i).seq;
    for (size_t c = 0; c < seq.length(); ++c) {
      char ch = toupper (seq[c]);   // ensure upper-case (alphabet_string is upper-case)
      if (is_dna && (ch == 'U'))    // handle case of 'U'
	ch = 'T';
      ++char_cnt[ch];
    }
    total += seq.length();
  }
  vector<double> char_freq (alphabet_string.length(), 0.);
  for (size_t c = 0; c < alphabet_string.length(); ++c)
    char_freq[c] = (DOUBLE_TINY + char_cnt[alphabet_string[c]]) / total;

  // log if requested
  if (CTAGGING(-1,FSAEM)) {
    CL << "Character frequencies:" << endl;
    // character label row
    for (size_t c = 0; c < alphabet_string.length(); ++c)
      CL << "char_freq[" << alphabet_string[c] << "] = " << char_freq[c] << endl;
  }

  return char_freq;
}

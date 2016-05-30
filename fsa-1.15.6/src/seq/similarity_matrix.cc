
/**
 * \file similarity_matrix.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <iostream>
#include <iomanip>

#include "math/mathematics.h"
#include "seq/similarity_matrix.h"

using namespace fsa;

Sequence_kmer_counts::Sequence_kmer_counts (const Alphabet& alphabet,
			  const size_t k)
  : __alphabet (alphabet),
    __k (k),
    __num_words (fsa::Mathematics::power_of_integer (__alphabet.size(), __k)) {

  // pre-compute necessary powers of alphabet size
  __powers_of_alphabet_size.resize (__k);
  for (size_t i = 0; i < k; ++i)
    __powers_of_alphabet_size[i] = Mathematics::power_of_integer (__alphabet.size(), i);

}

Sequence_kmer_counts::Kmer_counts Sequence_kmer_counts::compute_counts (const Sequence& sequence) const {

  Sequence_kmer_counts::Kmer_counts kmer_counts;

  // if the sequence is shorter than k, then return
  if (sequence.length() < k())
    return kmer_counts;

  // iterate over positions in the sequence, storing counts as we go
  for (size_t i = 0; i <= sequence.length() - k() + 1; ++i) {

    const Sequence_kmer_counts::Kmer_index index = kmer_to_index (sequence.seq, i);

    // if this word contains a character which isn't in the alphabet,
    // then don't store any information
    if (index == num_words())
      continue;

    // else store count
    else
      ++kmer_counts[index];

  }
  
  return kmer_counts;

}

size_t Sequence_kmer_counts::choose_minimum_word_length (const Sequence_database& seq_db,
							 const Alphabet& alphabet) {
  
  const size_t medianlen = seq_db.median_length();
  const size_t bestk = 1 + static_cast<size_t> (std::log (static_cast<double> (medianlen))
						/ std::log (static_cast<double> (alphabet.size())));

  if (bestk > 0)
    return bestk;
  else
    return 1;

}

Sequence_similarity_matrix::Sequence_similarity_matrix (const Sequence_database& seq_db,
							const Alphabet& alphabet,
							const size_t k)
  : __seq_db (seq_db),
    __similarities (seq_db.size(), std::vector<double> (seq_db.size(), 0.)) {

  compute_kmer_similarity_matrix (seq_db,
				  alphabet,
				  k);

}

size_t Sequence_similarity_matrix::compute_kmer_similarity (const Sequence_kmer_counts::Kmer_counts& counts_x, const Sequence_kmer_counts::Kmer_counts& counts_y) {

  size_t count = 0;

  // for speed, loop over the smaller set of k-mers
  if (counts_x.size() < counts_y.size()) {

    // for each k-mer in X
    for (Sequence_kmer_counts::Kmer_counts::const_iterator kmer_count_x = counts_x.begin(); kmer_count_x != counts_x.end(); ++kmer_count_x) {
      Sequence_kmer_counts::Kmer_counts::const_iterator kmer_count_y = counts_y.find (kmer_count_x->first);
      if (kmer_count_y != counts_y.end())
	count += std::min (kmer_count_x->second, kmer_count_y->second);
    }

  }

  else {

    // for each k-mer in Y
    for (Sequence_kmer_counts::Kmer_counts::const_iterator kmer_count_y = counts_y.begin(); kmer_count_y != counts_y.end(); ++kmer_count_y) {
      Sequence_kmer_counts::Kmer_counts::const_iterator kmer_count_x = counts_x.find (kmer_count_y->first);
      if (kmer_count_x != counts_x.end())
	count += std::min (kmer_count_x->second, kmer_count_y->second);
    }

  }

  return count;

}

void Sequence_similarity_matrix::write (std::ostream& o) const {

  // figure out an appropriate width based on the longest sequence name
  size_t width = 8;
  for (size_t i = 0; i < __seq_db.size(); ++i)
    width = std::max (width, __seq_db.get_seq (i).name.length() + 2);

  // precision of floating-point output
  const size_t precision = 3;

  // write header
  o.width (width); o << "";
  for (size_t i = 0; i < __similarities.size(); ++i) {
    o.width (width);
    o << __seq_db.get_seq (i).name;
  }
  o << endl;

  // matrix itself
  for (size_t i = 0; i < __similarities.size(); ++i) {
    o.width (width);
    o << __seq_db.get_seq (i).name;
    for (size_t j = 0; j < __similarities.size(); ++j) {
      o.width (width);
      if (i == j)
	o << "";
      else
	o << std::setprecision (precision) << get_similarity (i, j);
    }
    o << endl;
  }
  
}

void Sequence_similarity_matrix::compute_kmer_similarity_matrix (const Sequence_database& seq_db,
								 const Alphabet& alphabet,
								 const size_t k) {

  // check sane
  assert (seq_db.size() > 1);
  assert (seq_db.matches_alphabet (alphabet));

  // initialize appropriate Sequence_kmer_counts object
  const Sequence_kmer_counts counter (alphabet, k);

  // accumulate k-mer counts for each sequence
  std::vector<Sequence_kmer_counts::Kmer_counts> kmer_counts_db (seq_db.size());
  for (size_t i = 0; i < seq_db.size(); ++i)
    kmer_counts_db[i] = counter.compute_counts (seq_db.get_seq (i));

  // compute similarities
  assert (__similarities.size() == seq_db.size());
  for (size_t i = 0; i < seq_db.size(); ++i) {

    // check sane
    assert (__similarities[i].size() == seq_db.size());

    // set diagonal entries to one
    __similarities[i][i] = 1;

    // catch the case of a 0-length sequence
    if (seq_db.get_seq (i).length() == 0)
      continue;

    for (size_t j = i + 1; j < seq_db.size(); ++j) {

      // catch the case of a 0-length sequence
      if (seq_db.get_seq (j).length() == 0)
	continue;

      // compute k-mer similarity and normalize with the length of the shorter sequence
      __similarities[i][j] = static_cast<double> (compute_kmer_similarity (kmer_counts_db[i], kmer_counts_db[j]));
      __similarities[i][j] /=
	(std::min (seq_db.get_seq (i).length(), seq_db.get_seq (j).length()) >= k)
	? std::min (seq_db.get_seq (i).length(), seq_db.get_seq (j).length()) - (k - 1)
	: std::min (seq_db.get_seq (i).length(), seq_db.get_seq (j).length());

      // matrix is symmetric
      __similarities[j][i] = __similarities[i][j];

      // check sane
      if (__similarities[i][j] < 0.)
	cerr << __similarities[i][j] << endl;
      assert (__similarities[i][j] >= 0.);
      assert (__similarities[i][j] <= 1. + Mathematics::double_very_tiny);

    }

  }

}


/**
 * \file similarity_matrix.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef SEQ_SIMILARITY_MATRIX_INCLUDED
#define SEQ_SIMILARITY_MATRIX_INCLUDED

#include "util/misc.h"
#include "seq/alphabet.h"
#include "seq/sequence.h"

namespace fsa {

  /**
   * \brief Compute k-mer statistics of sequences.
   */
  struct Sequence_kmer_counts {
    
  public:

    typedef uintmax_t Kmer_index;                      ///< numerical index of a k-mer
    typedef std::map<Kmer_index, size_t> Kmer_counts;  ///< k-mer counts for a sequence

    /**
     * Constructor.
     * \param k word length
     */
    Sequence_kmer_counts (const Alphabet& alphabet,
			  const size_t k);

    /**
     * \brief Compute the k-mer counts of a sequence.
     */
    Kmer_counts compute_counts (const Sequence& sequence) const;
    
    /**
     * \brief Compute the minimum word length (>= 1) appropriate for computing distances.
     *
     * Given a uniform distribution over characters of the alphabet, compute
     *  k = log_|A| (L) + 1,  |A| = alphabet size
     * the minimum word length such that the expected number of occurrences
     * of a particular k-mer is < 1.
     * The length L is the median length of the sequences.
     */
    static size_t choose_minimum_word_length (const Sequence_database& seq_db,
					      const Alphabet& alphabet);

    /**
     * \brief Get number of possible words.
     * 
     * Note that this corresponds to 1 + the maximum index, e.g., a nonsense k-mer.
     * \see kmer_to_index
     */
    Kmer_index num_words() const { return __num_words; }

    /**
     * \brief Length of words.
     */
    size_t k() const { return __k; }

  private:

    /**
     * \brief Convert a substring of the passed string to the corresponding k-mer integral index.
     *
     * Considers the k-mer str.substr (start, __k).
     * It is the responsibility of the calling function to ensure that the
     * string boundaries are not overrun.
     *
     * \param str string to take a substring of
     * \param start position (0-based) of the k-mer
     * \return index or num_words() if the k-mer contains a character which is either degenerate or not in the alphabet
     */
    Kmer_index kmer_to_index (const std::string& str, const unsigned start) const;

    Alphabet __alphabet;              ///< alphabet which the k-mer are defined over
    unsigned short __k;               ///< word length
    Kmer_index __num_words;           ///< \see num_words

    std::vector<size_t> __powers_of_alphabet_size;  ///< \see kmer_to_index

  };

  /**
   * \brief Represent a matrix of similarities between pairs of sequences based on k-mer counts.
   */
  struct Sequence_similarity_matrix {

  public:
    
    /**
     * \brief Constructor.
     */
    Sequence_similarity_matrix (const Sequence_database& seq_db,
				const Alphabet& alphabet,
				const size_t k);

    /**
     * \brief Get the similarity between two sequences.
     */
    double get_similarity (const size_t i, const size_t j) const;

    /**
     * \brief Output method.
     */
    void write (std::ostream& o) const;

  private:

    /**
     * \brief Compute the k-mer similarity between two sequences.
     *
     * Computed as:
     * counts <- 0
     * For each k-mer w in X
     *   counts += min (# instances of w in X, # instances of w in Y)
     * (the outer loop is over k-mers of X or Y, depending upon which is the smaller set)
     *
     * \param counts_x k-mer counts for sequence X
     * \param counts_y k-mer counts for sequence Y
     */
    static size_t compute_kmer_similarity (const Sequence_kmer_counts::Kmer_counts& counts_x, const Sequence_kmer_counts::Kmer_counts& counts_y);

    /**
     * \brief Compute the k-mer similarity matrix between all pairs of sequences.
     *
     * The similarity between two sequences is defined as the
     * total number of shared k-mers divided by the length of the
     * shortest sequence.  We divide by the length of the shortest
     * sequence to avoid a systematic bias towards longer sequences.
     * Each entry in the similarity matrix therefore falls in the interval [0, 1].
     * Sets diagonal entries (the similarities of sequences with themselves) to 1.
     *
     * \param seq_db sequences to consider
     * \param alphabet alphabet over which sequences are defined
     * \param k word length to use
     * \see compute_kmer_similarity
     */
    void compute_kmer_similarity_matrix (const Sequence_database& seq_db,
					 const Alphabet& alphabet,
					 const size_t k);

    const Sequence_database& __seq_db;                  ///< sequence data
    std::vector<std::vector<double> > __similarities;   ///< sequence similarity matrix

  };

  inline Sequence_kmer_counts::Kmer_index Sequence_kmer_counts::kmer_to_index (const std::string& str, const unsigned start) const {

    // catch case of k = 0
    if (k() == 0)
      return 0;

    // assert sane
    assert (start + k() - 1 <= str.length());

    Sequence_kmer_counts::Kmer_index index = 0;
    for (size_t i = 0; i < k(); ++i) {

      const char ch = str[start + i];

      if (!__alphabet.is_nondegen_char (ch))
	return num_words();

      // exploit 2-bit representation of a DNA alphabet if possible (the hare)
      if (__alphabet.size() == 4) {
	index <<= 2; // because 2 bits per character
	index |= __alphabet.get_char_index (ch);
      }

      // else do it with multiplication (the tortoise)
      // (but not too slow, since we've hardcoded in powers of the alphabet size)
      // note that the "backwards" vector indexing is so that the leftmost character,
      // which we see first, corresponds to the largest power of the alphabet size
      // (which is how the above bit-shifting works)
      else
	index += __alphabet.get_char_index (ch) * __powers_of_alphabet_size[k() - i - 1];

    }

    assert (index < num_words());

    return index;

  }

  inline double Sequence_similarity_matrix::get_similarity (const size_t i, const size_t j) const {
    assert (i < __similarities.size());
    assert (j < __similarities.size());
    return __similarities[i][j];
  }

}

#endif /* SEQ_SIMILARITY_MATRIX_INCLUDED */

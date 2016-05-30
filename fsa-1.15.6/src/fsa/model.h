
/**
 * \file model.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef FSA_MODEL_INCLUDED
#define FSA_MODEL_INCLUDED

#define NUM_HMM_STATES 7

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "util/misc.h"
#include "util/opts_list.h"
#include "util/memcheck.h"
#include "seq/sequence.h"
#include "seq/alignment.h"

#include "fsa/mybanding.h"
#include "fsa/nucleotidedp.h"
#include "fsa/aminoaciddp.h"
#include "fsa/nucleotide_indel2dp.h"
#include "fsa/aminoacid_indel2dp.h"
#include "fsa/sequence_pair_selector.h"

#include "annealing/dotplot.h"
#include "annealing/SparseMatrix.h"
#include "annealing/alignment_DAG.h"

namespace fsa {

  /**
   * \param A struct to hold the alignment parameters and the DP bandwidth.
   *
   * Distributions aren't required to be normalized, so they can be used to hold pseudocounts as well.
   */
  struct Params {

  public:
    size_t bandwidth;     ///< banding width (of DP matrix)

    double gap_open1;     ///< gap-open probability (set 1 of indel states)
    double gap_extend1;   ///< gap-extend probability (set 1 of indel states)

    double gap_open2;     ///< gap-open probability (set 2 of indel states)
    double gap_extend2;   ///< gap-extend probability (set 2 of indel states)

    double to_end;        ///< probability to transition to End state of HMM (basically irrelevant)

    double time;                                 ///< evolutionary time parameter (for Jukes-Cantor & Tamura-Nei)
    std::vector<double> single_dist;               ///< distribution over single nucleotides (for Insert and Delete states)
    std::vector<std::vector<double> > pair_dist;   ///< distribution over pairs of aligned nucleotides (for Match state)

    /**
     * \brief Alphabet of the HMM.
     *
     * We must store the alphabet like this for sorting later on.
     * (HMMoC-generated DP code demands alphabetical sorting.)
     */
    std::string alphabet_string;

    bool is_indel2;  ///< one or two sets of indel states?

    std::vector<std::string> states; ///< states of Pair HMM

    /**
     * \brief Transition probabilities of Pair HMM.
     *
     *     0     1       2       3       4       5   6
     * Start Match Insert1 Delete1 Insert2 Delete2 End
     */
    std::vector<std::vector<double> > transition_matrix;

    /**
     * \brief Stationary distribution over states of the Pair HMM.
     */
    std::vector<double> stat_dist;

    /**
     * \brief Constructor.
     */
    Params();

    /**
     * \brief Initialize transition probabilities of Pair HMM.
     *
     * Initialize the transition matrix with the desired HMM structure.
     * If left_match is set, then a transition Start -> Match is forced.
     * This is appropriate if you know there is homology to the "left."
     * If right_match is set, then only the {Insert, Delete} -> End transitions
     * are disallowed.  This is appropriate if you know
     * there is homology to the "right."
     * Note that this transition matrix is NOT symmetric, but rather
     * favors aligments which "begin" with homologous regions and end
     * with possibly non-homologous regions.  This is intentional and reflects
     * the implicit assumptions most people make when constructing alignments.
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     */
    void init_transition_matrix (bool left_match = false, bool right_match = false, bool ragged_ends = false);

    /**
     * \brief Normalize transition matrix to be a (left) stochastic matrix.
     */
    void normalize_transition_matrix();

    /**
     * \brief Assert that entries in the transition matrix are valid probabilities.
     */
    void assert_transition_matrix_valid();

    /**
     * \brief Copy all parameters.
     *
     * \param from Params template
     */
    void copy_all (const Params& from);

    /**
     * \brief Copy indel parameters (for transition probabilities).
     *
     * \param from Params template
     */
    void copy_indel (const Params& from);

    /**
     * \brief Copy emission distributions.
     *
     * \param from Params template
     */
    void copy_subst (const Params& from);

    /**
     * \brief Show indel parameters and the corresponding transition matrix.
     */
    void show_transition_matrix (std::ostream& o) const;

    /**
     * \brief Show emit distributions in human-readable format.
     */
    void show_emission (std::ostream& o) const;

    /**
     * \brief Show all parameters.
     */
    void show (std::ostream& o) const;

    /**
     * \brief Write transition probs to disk.
     * \param filename filename
     */
    void write_transition (const std::string& filename) const;

    /**
     * \brief Write emission probs to disk.
     * \param filename filename
     */
    void write_emission (const std::string& filename) const;

    /**
     * \brief Write emission probs in .dat-format file to disk for viewParams.pl to use.
     * \param filename filename
     */
    void write_emission_dat (const std::string& filename) const;

    /**
     * \brief Are the parameters biologically meaningful?
     *
     * Returns false if odds (match) < odds (mismatch),
     * where odds (mismatch) = P (a ~ c) / (P (a) * P (c)), etc.
     */
    bool is_biological() const;

    /**
     * \brief Compute branch-length based on log-det transform.
     *
     * which specifies which index (0 = first => columns, 1 = second => rows) of pair_dist
     * is the ancestral sequence.  It follows that
     *  which = 0 gives a right stochastic matrix,
     *     where summing over the second index (rows) of the conditional matrix submat gives 1.
     *  which = 1 gives a left stochastic matrix,
     *     where summing over the first index (columns) of the conditional matrix submat gives 1.
     * \param which specifies which index is the ancestral sequence
     */
    double branch_length (unsigned which) const;

    /**
     * \brief Compute -log det.
     *
     * \return -1 if undefined (if determinant < 0).
     */
    static double logdet (const std::vector<std::vector<double> >& matrix);

    /**
     * \brief Compute the determinant of the substitution matrix.
     * \param matrix matrix of interest
     */
    static double determinant (const std::vector<std::vector<double> >& matrix);

    /**
     * \brief Assert distributions normalized.
     */
    void assert_normalized() const;

    /**
     * \brief Normalize distributions.
     */
    void normalize();

  };

  /**
   * \brief The statistical model.
   */
  struct Model {

    /**
     * \brief Create alignment dotplot from DP matrices.
     *
     * \param params parameters
     * \param xseq first sequence X
     * \param yseq second sequence Y
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     * \return Dotplot
     */
    template<class T>
    static Dotplot get_pairwise_dotplot_engine (Params& params, const Sequence& xseq, const Sequence& yseq,
						bool left_match, bool right_match, bool ragged_ends);

    /**
     * \brief Create an alignment posterior probability Post_probs from DP matrices.
     *
     * The construction of the Post_probs object guarantees proper lexical sorting.
     * \param params parameters
     * \param xseq first sequence X
     * \param yseq second sequence Y
     * \param minprob minimum posterior probability to store
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     * \return Post_probs (format expected by SparseMatrix constructor).
     */
    template<class T>
    static Post_probs get_pairwise_post_probs_engine (Params& params, const Sequence& xseq, const Sequence& yseq,
						      double minprob,
						      bool left_match, bool right_match, bool ragged_ends);

    /**
     * \brief Use Baum-Welch training to estimate transition (indel) and emission parameters.
     *
     * Unsupervised training on the passed list of sequences.
     * Returns the total forward likelihood.
     * Dirichlet regularization according to the passed pseudocounts (held in a Params object, but not normalized).
     * \param params parameters
     * \param seq_db sequence data
     * \param seq_pairs sequence pairs to consider
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     * \param learn_gap learn indel parameters
     * \param learn_emit learn emit parameters
     * \param regularize regularize parameters
     * \param pseudocounts pseudocounts for emission regularization with a Dirichlet
     * \param iterations maximum number of rounds of EM
     * \param min_inc minimum fractional increase in log-likelihood before termination
     */
    template<class T>
    static bool train_params_engine (Params& params, const Sequence_database& seq_db, const Sequence_pairs& seq_pairs,
				     bool left_match, bool right_match, bool ragged_ends,
				     bool learn_gap, bool learn_emit, bool regularize, const Params& pseudocounts, size_t iterations, double min_inc);

    /**
     * \brief Get character frequencies.
     *
     * Force tiny value for all characters, even if absent.
     * NB: We CANNOT use seq_db.get_null_model (alphabet.size())
     * because it will return an answer sorted non-alphabetically =>
     * will give us WRONG results!
     * Assumes that the alphabet is not case-sensitive and that the
     * alphabet_string is upper-case.
     * \param seq_db sequence data
     * \param alphabet_string (ordered) characters in the alphabet
     * \param is_dna is the alphabet DNA or protein
     */
    static std::vector<double> count_char_freq (const Sequence_database& seq_db, const std::string& alphabet_string, bool is_dna);

  };

  template<class T>
    Dotplot Model::get_pairwise_dotplot_engine (Params& params, const Sequence& xseq, const Sequence& yseq,
						bool left_match, bool right_match, bool ragged_ends) {

    // print log message
    if (CTAGGING(6,FSA))
      CTAG(6,FSA) << "Summing over alignments of '" << xseq.name << "' and '" << yseq.name << "' to get posterior probs." << endl;

    const size_t xlen = xseq.length();
    const size_t ylen = yseq.length();

    // initialize transition matrix
    params.init_transition_matrix (left_match, right_match, ragged_ends);

    // ensure that left_match and right_match are reasonable
    assert (xlen >= static_cast<size_t> (left_match) + static_cast<size_t> (right_match));
    assert (ylen >= static_cast<size_t> (left_match) + static_cast<size_t> (right_match));

    // calculate the Forward and Backward tables, to compute posteriors
    typename T::DPTable *pFW, *pBW;
    const bfloat fw = T::Forward (&pFW, params, xseq.seq, yseq.seq);
    T::Backward (&pBW, params, xseq.seq, yseq.seq);

    if (CTAGGING(5,FSA))
      CTAG(5,FSA) << "Per-site log-likelihood = " << log (fw) / (xseq.length() + yseq.length()) << "." << endl;

    Dotplot dotplot (xseq, yseq);
    const std::string inner_seq = "sequence2";
    const int match_id = pBW->getId ("match");

    // the first coordinate in getProb refers to the fastest coordinate (inner loop), sequence2;
    // therefore switch x and y when calling getProb if necessary
    if (pFW->outputId[0] == inner_seq) {

      const bfloat total = pBW->getProb ("start", 0, 0);
      for (unsigned x = 0; x < xlen; ++x) {
	for (unsigned y = 0; y < ylen; ++y) {
	  double& p = dotplot (x, y);
	  p += pFW->getProb (match_id, y + 1, x + 1) * pBW->getProb (match_id, y + 1, x + 1) / total; // HMMoC DP accessor code is 1-based
	}
      }

    }
    // else don't
    else {

      const bfloat total = pBW->getProb ("start", 0, 0);
      for (unsigned x = 0; x < xlen; ++x) {
	for (unsigned y = 0; y < ylen; ++y) {
	  double& p = dotplot (x, y);
	  p += pFW->getProb (match_id, x + 1, y + 1) * pBW->getProb (match_id, x + 1, y + 1) / total; // HMMoC DP accessor code is 1-based
	}
      }

    }

    // finished, so clean up
    delete pFW;
    delete pBW;

    return dotplot;
  }

  template<class T>
    Post_probs Model::get_pairwise_post_probs_engine (Params& params, const Sequence& xseq, const Sequence& yseq,
						      double minprob,
						      bool left_match, bool right_match, bool ragged_ends) {
  
    // print log message
    if (CTAGGING(6,FSA))
      CTAG(6,FSA) << "Summing over alignments of '" << xseq.name << "' and '" << yseq.name << "' to get posterior probs." << endl;

    const size_t xlen = xseq.length();
    const size_t ylen = yseq.length();

    // initialize transition matrix
    params.init_transition_matrix (left_match, right_match, ragged_ends);

    // ensure that left_match and right_match are reasonable
    assert (xlen >= static_cast<size_t> (left_match) + static_cast<size_t> (right_match));
    assert (ylen >= static_cast<size_t> (left_match) + static_cast<size_t> (right_match));

    // calculate the Forward and Backward tables, to compute posteriors
    typename T::DPTable *pFW, *pBW;
    const bfloat fw = T::Forward (&pFW, params, xseq.seq, yseq.seq);
    T::Backward (&pBW, params, xseq.seq, yseq.seq);

    if (CTAGGING(5,FSA))
      CTAG (5,FSA) << "Per-site log-likelihood = " << log (fw) / (xlen + ylen) << "." << endl;

    // first assemble posterior vector from post prob information
    Post_probs post_probs;

    const std::string inner_seq = "sequence2";
    const int match_id = pBW->getId ("match");

    // the first coordinate in getProb refers to the fastest coordinate (inner loop), sequence2;
    // therefore switch x and y when calling getProb if necessary
    if (pFW->outputId[0] == inner_seq) {

      const bfloat total = pBW->getProb ("start", 0, 0);
      for (unsigned x = 0; x < xlen; ++x) {
	for (unsigned y = 0; y < ylen; ++y) {
	  double p = pFW->getProb (match_id, y + 1, x + 1) * pBW->getProb (match_id, y + 1, x + 1) / total; // HMMoC DP accessor code is 1-based
	  // sparse storage
	  if (p >= minprob)
	    post_probs.push_back (Post_prob (x, y, p));
	}
      }

    }
    // else don't
    else {

      const bfloat total = pBW->getProb ("start", 0, 0);
      for (unsigned x = 0; x < xlen; ++x) {
	for (unsigned y = 0; y < ylen; ++y) {
	  double p = pFW->getProb (match_id, x + 1, y + 1) * pBW->getProb (match_id, x + 1, y + 1) / total; // HMMoC DP accessor code is 1-based
	  // sparse storage
	  if (p >= minprob)
	    post_probs.push_back (Post_prob (x, y, p));
	}
      }

    }

    // finished, so clean up
    delete pFW;
    delete pBW;

    return post_probs;
  }

  template<class T>
    bool Model::train_params_engine (Params& params, const Sequence_database& seq_db, const Sequence_pairs& seq_pairs,
				     bool left_match, bool right_match, bool ragged_ends,
				     bool learn_gap, bool learn_emit, bool regularize, const Params& pseudocounts, size_t iterations, double min_inc) {

    // get number of seqs
    const size_t num_seqs = seq_db.size();
    const size_t num_seq_pairs = seq_pairs.size();

    // log
    assert (num_seqs >= 2);
    std::string what_training;
    if (learn_gap && learn_emit)
      what_training = "gap and emit";
    else if (learn_gap && !learn_emit)
      what_training = "gap";
    else if (!learn_gap && learn_emit)
      what_training = "emit";
    else
      what_training = "nothing; why am I training?";
    if ((num_seqs == 2) && CTAGGING(6,FSAEM))
      CTAG(6,FSAEM) << "Training parameters (" << what_training << ") for '" << seq_db.get_seq (0).name << "' and '" << seq_db.get_seq (1).name << "'." << endl;
    else if ((num_seqs != 2) && CTAGGING(9,FSAEM))
      CTAG(9,FSAEM) << "Training parameters (" << what_training << ") for " << seq_db.size() << " sequences (using " << num_seq_pairs << "/" << num_seqs * (num_seqs - 1) / 2 << " possible sequence pairs)." << endl;

    // get alphabet size and alphabet
    const size_t alphabet_size = params.alphabet_string.length();

    // init DP tables
    typename T::DPTable* pFW;
    typename T::BaumWelchCounters baum_welch;

    bfloat curr_fw_total = 0;
    bfloat prev_fw_total = 0;

    const std::string emit1 = "emit1";
    const std::string emit2 = "emit2";
    const std::string emit12 = "emit12";

    for (size_t iter = 1; iter <= iterations; ++iter) { // 1-based for easy logging
    
      // initialize transition matrix
      params.init_transition_matrix (left_match, right_match, ragged_ends);

      // reset counts and begin round of EM
      baum_welch.resetCounts();
      curr_fw_total = 0;
      size_t sites = 0;

      if ((num_seqs == 2) && CTAGGING(6,FSAEM))
	CTAG(6,FSAEM) << "Beginning EM round " << iter << "/" << iterations << endl;
      else if ((num_seqs != 2) && CTAGGING(7,FSAEM))
	CTAG(7,FSAEM) << "Beginning EM round " << iter << "/" << iterations << endl;

      // accumulate counts for all sequence pairs
      for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt) {
	const size_t i = seq_pairs[cnt].first;
	const size_t j = seq_pairs[cnt].second;

	sites += seq_db.get_seq (i).length() + seq_db.get_seq (j).length();

	// ensure that left_match and right_match are reasonable
	assert (seq_db.get_seq (i).length() >= static_cast<size_t> (left_match) + static_cast<size_t> (right_match));
	assert (seq_db.get_seq (j).length() >= static_cast<size_t> (left_match) + static_cast<size_t> (right_match));

	// calculate the forward DP table
	const bfloat fw = T::Forward (&pFW, params, seq_db.get_seq (i).seq, seq_db.get_seq (j).seq);
	curr_fw_total += fw;

	if (CTAGGING(4,FSAEM))
	  CTAG(4,FSAEM) << "Adding counts for '" << seq_db.get_seq (i).name << "' and '" << seq_db.get_seq (j).name << "' (per-site log-likelihood = " << log (fw) / (seq_db.get_seq (i).length() + seq_db.get_seq (j).length()) << ")" << endl;

	// calculate the Baum-Welch estimated transition counts
	const bfloat bw = T::BackwardBaumWelch (baum_welch, pFW, params, seq_db.get_seq (i).seq, seq_db.get_seq (j).seq);

	// log progress
	if (num_seqs > 2) {
	  const unsigned percent_done = static_cast<unsigned> (std::floor ((100.0 * (cnt+1) / num_seq_pairs) + 0.5));
	  if (CTAGGING(7,FSA))
	    CTAG(7,FSA) << "Processed sequence pair '" << seq_db.get_seq (i).name << "' and '" << seq_db.get_seq (j).name << "' (" << percent_done << "% complete)." << endl;
	}

	// remove the forward table
	delete pFW;
      }

      if (learn_gap) {

	// get the expected transition counts
	const double trMM = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trMM")];
	const double trMI1 = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trMI1")];
	const double trMD1 = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trMD1")];

	// these guys don't exist if only 1 set of indel states, so do this to avoid a segfault...
	const double trMI2 = params.is_indel2 ? (double) baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trMI2")] : 0.;
	const double trMD2 = params.is_indel2 ? (double) baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trMD2")] : 0.;

	const double trI1I1 = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trI1I1")];
	const double trI1M = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trI1M")];
	const double trD1D1 = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trD1D1")];
	const double trD1M = baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trD1M")];

	const double trI2I2 = params.is_indel2 ? (double) baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trI2I2")] : 0.;
	const double trI2M = params.is_indel2 ? (double) baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trI2M")] : 0.;
	const double trD2D2 = params.is_indel2 ? (double) baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trD2D2")] : 0.;
	const double trD2M = params.is_indel2 ? (double) baum_welch.transitionBaumWelchCount00[baum_welch.transitionIndex ("trD2M")] : 0.;

	// add DOUBLE_TINY to ensure sanity
	double match_to_indel1 = trMI1 + trMD1 + DOUBLE_TINY;
	double indel1_to_indel1 = trI1I1 + trD1D1 + DOUBLE_TINY;
	double match_to_indel2 = trMI2 + trMD2 + DOUBLE_TINY;
	double indel2_to_indel2 = trI2I2 + trD2D2 + DOUBLE_TINY;
	double from_match = trMM + match_to_indel1 + match_to_indel2 + DOUBLE_TINY;
	double from_indel1 = trI1I1 + trI1M + trD1D1 + trD1M + DOUBLE_TINY;
	double from_indel2 = trI2I2 + trI2M + trD2D2 + trD2M + DOUBLE_TINY;

	// log transition counts
	if (CTAGGING(-1,FSAEM)) {
	  CL << "Transition counts:" << endl
	     << "match_to_indel1 = " << match_to_indel1 << endl
	     << "match_to_indel2 = " << match_to_indel2 << endl
	     << "indel1_to_indel1 = " << indel1_to_indel1 << endl
	     << "indel2_to_indel2 = " << indel2_to_indel2 << endl
	     << "from_match = " << from_match << endl
	     << "from_indel1 = " << from_indel1 << endl
	     << "from_indel2 = " << from_indel2 << endl;
	}

	// regularization
	if (regularize) {
	  // add pseudocounts, scaled by the number of sequence pairs
	  if (CTAGGING(2,FSAEM))
	    CL << "Adding transition pseudocounts" << endl;
	  match_to_indel1 += pseudocounts.gap_open1 * num_seq_pairs;
	  match_to_indel2 += pseudocounts.gap_open2 * num_seq_pairs;
	  from_match += pseudocounts.gap_open1 * num_seq_pairs + pseudocounts.gap_open2 * num_seq_pairs;
	  indel1_to_indel1 += pseudocounts.gap_extend1 * num_seq_pairs;
	  indel2_to_indel2 += pseudocounts.gap_extend2 * num_seq_pairs;
	  from_indel1 += pseudocounts.gap_extend1 * num_seq_pairs;
	  from_indel2 += pseudocounts.gap_extend2 * num_seq_pairs;

	  // log counts post-regularization
	  if (CTAGGING(-1,FSAEM)) {
	    CL << "Transition counts after regularization:" << endl
	       << "match_to_indel1 = " << match_to_indel1 << endl
	       << "match_to_indel2 = " << match_to_indel2 << endl
	       << "indel1_to_indel1 = " << indel1_to_indel1 << endl
	       << "indel2_to_indel2 = " << indel2_to_indel2 << endl
	       << "from_match = " << from_match << endl
	       << "from_indel1 = " << from_indel1 << endl
	       << "from_indel2 = " << from_indel2 << endl;
	  }
	}

	// re-estimate indel parameters
	params.gap_open1 = match_to_indel1 / from_match;
	params.gap_extend1 = indel1_to_indel1 / from_indel1;
	params.gap_open2 = match_to_indel2 / from_match;
	params.gap_extend2 = indel2_to_indel2 / from_indel2;

	// force indel parameters to be sane:
	// Problems can arise b/c the Insert->Match transitions are calculated internally in the DP code as 
	//    1.0 - gap_extend - to_end
	// and Match->Match probs as
	//    1.0 - gap_open - to_end
	// If gap_extend becomes 1 (as is wont to happen), then this quantity becomes negative and Very Bad Things can happen.
	// Make sure that it doesn't:
	// Force all transition probabilities to be in the range [DOUBLE_TINY, 1.0 - DOUBLE_TINY].

	// gap-open probability
	params.gap_open1 = (params.gap_open1 > (1.0 - DOUBLE_TINY)) ? 1.0 - DOUBLE_TINY : params.gap_open1;
	params.gap_open1 = (params.gap_open1 < DOUBLE_TINY) ? DOUBLE_TINY : params.gap_open1;
	params.gap_open2 = (params.gap_open2 > (1.0 - DOUBLE_TINY)) ? 1.0 - DOUBLE_TINY : params.gap_open2;
	params.gap_open2 = (params.gap_open2 < DOUBLE_TINY) ? DOUBLE_TINY : params.gap_open2;

	// gap-extend probability
	params.gap_extend1 = (params.gap_extend1 > (1.0 - DOUBLE_TINY)) ? 1.0 - DOUBLE_TINY : params.gap_extend1;
	params.gap_extend1 = (params.gap_extend1 < DOUBLE_TINY) ? DOUBLE_TINY : params.gap_extend1;
	params.gap_extend2 = (params.gap_extend2 > (1.0 - DOUBLE_TINY)) ? 1.0 - DOUBLE_TINY : params.gap_extend2;
	params.gap_extend2 = (params.gap_extend2 < DOUBLE_TINY) ? DOUBLE_TINY : params.gap_extend2;

	// gap-close probability
	if (1.0 - params.gap_extend1 - DOUBLE_TINY < DOUBLE_TINY) // hardcode in parameter to_end = DOUBLE_TINY
	  params.gap_extend1 -= DOUBLE_TINY;
	if (1.0 - params.gap_extend2 - DOUBLE_TINY < DOUBLE_TINY) // hardcode in parameter to_end = DOUBLE_TINY
	  params.gap_extend2 -= DOUBLE_TINY;

	// match->match probability
	if (1.0 - params.gap_open1 - params.gap_open2 - DOUBLE_TINY < 0) { // hardcode in parameter to_end = DOUBLE_TINY
	  params.gap_open1 -= DOUBLE_TINY / 2;
	  params.gap_open2 -= DOUBLE_TINY / 2;
	}

	// assert sane
	assert (params.gap_open1 >= 0. && params.gap_open1 <= 1.);
	assert (params.gap_extend1 >= 0. && params.gap_extend1 <= 1.);
	assert (params.gap_open2 >= 0. && params.gap_open2 <= 1.);
	assert (params.gap_extend2 >= 0. && params.gap_extend2 <= 1.);
	
      }

      // re-estimate emission parameters
      if (learn_emit) {

	// accumulate total counts for normalization and log expected counts
	double singletotal = 0;
	double pairtotal = 0;
	if (CTAGGING(-1,FSAEM))
	  CL << "Emission counts:" << endl;
	for (size_t i = 0; i < alphabet_size; ++i) {
	  if (CTAGGING(-1,FSAEM)) {
	    CL << " " << emit1 << "[" << params.alphabet_string[i] << "]: " << baum_welch.emissionBaumWelchCount10[i][baum_welch.emissionIndex (emit1)] << endl
	       << " " << emit2 << "[" << params.alphabet_string[i] << "]: " << baum_welch.emissionBaumWelchCount01[i][baum_welch.emissionIndex (emit2)] << endl;
	  }
	  singletotal += baum_welch.emissionBaumWelchCount10[i][baum_welch.emissionIndex (emit1)] + baum_welch.emissionBaumWelchCount01[i][baum_welch.emissionIndex (emit2)];
	}
	for (size_t i = 0; i < alphabet_size; ++i) {
	  for (size_t j = 0; j < alphabet_size; ++j) {
	    if (CTAGGING(-1,FSAEM))
	      CL << " " << emit12 << "[" << params.alphabet_string[i] << params.alphabet_string[j] << "]: " << baum_welch.emissionBaumWelchCount11[i][j][baum_welch.emissionIndex (emit12)] << endl;
	    pairtotal += baum_welch.emissionBaumWelchCount11[i][j][baum_welch.emissionIndex (emit12)];
	  }
	}
      
	// regularization
	if (regularize) {
	  // add pseudocounts, scaled by the number of sequence pairs
	  if (CTAGGING(2,FSAEM))
	    CL << "Adding emission pseudocounts" << endl;
	  for (size_t i = 0; i < alphabet_size; ++i) {
	    baum_welch.emissionBaumWelchCount10[i][baum_welch.emissionIndex (emit1)] += pseudocounts.single_dist[i] * num_seq_pairs; // scale by number of sequence pairs b/c pseudocounts are for a single pairwise comparison
	    baum_welch.emissionBaumWelchCount01[i][baum_welch.emissionIndex (emit2)] += pseudocounts.single_dist[i] * num_seq_pairs;
	    singletotal += 2 * pseudocounts.single_dist[i] * num_seq_pairs;
	    for (size_t j = 0; j < alphabet_size; ++j) {
	      baum_welch.emissionBaumWelchCount11[i][j][baum_welch.emissionIndex (emit12)] += pseudocounts.pair_dist[i][j] * num_seq_pairs;
	      pairtotal += pseudocounts.pair_dist[i][j] * num_seq_pairs;
	    }
	  }

	  // log counts post-regularization
	  if (CTAGGING(-1,FSAEM)) {
	    CL << "Emission counts after regularization:" << endl;
	    for (size_t i = 0; i < alphabet_size; ++i) {
	      CL << " " << emit1 << "[" << params.alphabet_string[i] << "]: " << baum_welch.emissionBaumWelchCount10[i][baum_welch.emissionIndex (emit1)] << endl
		 << " " << emit2 << "[" << params.alphabet_string[i] << "]: " << baum_welch.emissionBaumWelchCount01[i][baum_welch.emissionIndex (emit2)] << endl;
	    }
	    for (size_t i = 0; i < alphabet_size; ++i)
	      for (size_t j = 0; j < alphabet_size; ++j)
		CL << " " << emit12 << "[" << params.alphabet_string[i] << params.alphabet_string[j] << "]: " << baum_welch.emissionBaumWelchCount11[i][j][baum_welch.emissionIndex (emit12)] << endl;
	  }

	}

	// re-estimate emission params
	for (size_t i = 0; i < alphabet_size; ++i) {
	  params.single_dist[i] = (baum_welch.emissionBaumWelchCount10[i][baum_welch.emissionIndex (emit1)] + baum_welch.emissionBaumWelchCount01[i][baum_welch.emissionIndex (emit2)]) / singletotal;
	  for (size_t j = 0; j < alphabet_size; ++j)
	    params.pair_dist[i][j] = baum_welch.emissionBaumWelchCount11[i][j][baum_welch.emissionIndex (emit12)] / pairtotal;
	}

      }

      // log all new params
      if (CTAGGING(1,FSAEM)) {
	CL << "Parameter estimates after EM round " << iter << "/" << iterations << endl;
	if (learn_gap)
	  params.show_transition_matrix (CL);
	if (learn_emit)
	  params.show_emission (CL);
      }

      params.assert_normalized();

      // check that fractional increase in log-likelihood is sufficient to continue
      // (conditional to cover case of first iteration)
      // in below function, we DON'T use std::log in order to get compiler to call Algebra<BFloatMethods>'s friend method log 
      const double inc = std::abs (static_cast<double> (log (curr_fw_total) / ((std::abs (log (prev_fw_total)) < DOUBLE_TINY) ? 1.0 : log (prev_fw_total)) - 1.0));
      if ((curr_fw_total <= prev_fw_total) || (inc < min_inc))	// fractional improvement too low, so stop
	{
	  // did our probabilities fail to increase?
	  if ((curr_fw_total + DOUBLE_TINY < prev_fw_total) && CTAGGING(7,FSAEM)) // ignore rounding errors
	    CL << "WARNING: per-site log-likelihood dropped from " << log (prev_fw_total) / sites
	       << " to " << log (curr_fw_total) / sites << " during EM." << endl;
	  if (CTAGGING(4,FSAEM))
	    CL << "Failed EM fractional improvement threshold of " << min_inc << "; stopping (inc = " << inc << ")." << endl;
	  break;
	}
      prev_fw_total = curr_fw_total;

    }

    return true;
  }


  /**
   * \brief Accessor classes for the standard DP algorithms (nucleotide).
   */
  class NucleotideWithoutBanding {

  public:
  
    typedef NucleotideAlignDPTable    DPTable;
    typedef NucleotideAlignBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Forward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Backward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelch (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

  };

  /**
   * \brief Accessor classes for the banded DP algorithms (nucleotide).
   */
  class NucleotideWithBanding {

  public:

    typedef NucleotideAlignWithBandingDPTable    DPTable;
    typedef NucleotideAlignWithBandingBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::ForwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelchBanding (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

  };

  /**
   * \brief Accessor classes for the standard DP algorithms (amino acid).
   */
  class AminoAcidWithoutBanding {

  public:
  
    typedef AminoAcidAlignDPTable    DPTable;
    typedef AminoAcidAlignBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Forward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Backward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelch (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

  };

  /**
   * \brief Accessor classes for the banded DP algorithms (amino acid).
   */
  class AminoAcidWithBanding {

  public:

    typedef AminoAcidAlignWithBandingDPTable    DPTable;
    typedef AminoAcidAlignWithBandingBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::ForwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelchBanding (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

  };

  /**
   * \brief Accessor classes for the standard DP algorithms (nucleotide with 2 sets of indel states).
   */
  class NucleotideIndel2WithoutBanding {

  public:
  
    typedef NucleotideIndel2AlignDPTable    DPTable;
    typedef NucleotideIndel2AlignBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Forward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Backward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelch (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

  };

  /**
   * \brief Accessor classes for the banded DP algorithms (nucleotide with 2 sets of indel states).
   */
  class NucleotideIndel2WithBanding {

  public:

    typedef NucleotideIndel2AlignWithBandingDPTable    DPTable;
    typedef NucleotideIndel2AlignWithBandingBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::ForwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelchBanding (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

  };

  /**
   * \brief Accessor classes for the standard DP algorithms (amino acid with 2 sets of indel states).
   */
  class AminoAcidIndel2WithoutBanding {

  public:
  
    typedef AminoAcidIndel2AlignDPTable    DPTable;
    typedef AminoAcidIndel2AlignBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Forward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::Backward (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelch (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix);
    }

  };

  /**
   * \brief Accessor classes for the banded DP algorithms (amino acid with 2 sets of indel states).
   */
  class AminoAcidIndel2WithBanding {

  public:

    typedef AminoAcidIndel2AlignWithBandingDPTable    DPTable;
    typedef AminoAcidIndel2AlignWithBandingBaumWelch  BaumWelchCounters;

    static bfloat Forward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::ForwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat Backward (DPTable** ppOutTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBanding (ppOutTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

    static bfloat BackwardBaumWelch (BaumWelchCounters& bw, DPTable* pInTable, Params& params, const std::string& iSequence1, const std::string& iSequence2) {
      return ::BackwardBaumWelchBanding (bw, pInTable, iSequence1, iSequence2, params.single_dist, params.pair_dist, params.transition_matrix, params.bandwidth);
    }

  };

}

#endif /* FSA_MODEL_INCLUDED */

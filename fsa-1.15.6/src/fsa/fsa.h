
/**
 * \file fsa.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 * Jaeyoung Do wrote the parallelization and database code.
 *
 * All nucleotide sequences are treated internally as DNA, meaning that 'U' gets converted to 'T'
 * (this is done when degenerate characters are randomized in build_multiple_alignment).
 * Messing this up will make bad things happen!
 * NB: The case of the input sequences is preserved (in order to not destroy potential softmasking);
 * however, the inference engines treat upper and lower-case characters identically.
 *
 * Dirichlet regularization scales:
 * The emission regularization scales correspond precisely to the total number of
 * pseudocount emissions because the seed distribution for pseudocount calculation
 * is normalized to 1.
 * The default values are (approximately) equal to the number of free parameters in a
 * symmetric pair emission matrix
 *   = 1/2 * n * (n - 1) + n
 *   => 1/2 * 4 * 3 + 4 = 10 for nucleotides
 *   => 1/2 * 20 * 19 + 20 = 210 for amino acids
 * The transition regularization scale is chosen completely arbitrarily.
 */

#ifndef FSA_FSA_INCLUDED
#define FSA_FSA_INCLUDED


#define DNA_ALPHABET_STRING "ACGT"
#define PROTEIN_ALPHABET_STRING "ACDEFGHIKLMNPQRSTVWY"

#define MODEL_DEFAULT 1
#define TIME_DEFAULT 0.4
#define ALPHA_R_DEFAULT 1.3
#define ALPHA_Y_DEFAULT 1.3
#define BETA_DEFAULT 1.

#define EM_MIN_INC_DEFAULT 0.1
#define EM_MAX_ITER_DEFAULT 3
#define REGULARIZATION_NUC_TRANSITION_SCALE_DEFAULT 3.
#define REGULARIZATION_NUC_EMISSION_SCALE_DEFAULT 10.
#define REGULARIZATION_AA_TRANSITION_SCALE_DEFAULT 3.
#define REGULARIZATION_AA_EMISSION_SCALE_DEFAULT 210.

#define MIN_NUC_EMISSION_TRAINING_DATA 60
#define MIN_NUC_TRANSITION_TRAINING_DATA 60
#define MIN_AA_EMISSION_TRAINING_DATA 1596
#define MIN_AA_TRANSITION_TRAINING_DATA 60

#define LONG_DNA_DEFAULT 500                ///< length to trigger anchoring

#define DEFAULT_POSTERIOR_PROB_CUTOFF 0.01  ///< minimum posterior probability value that is maintained in the sparse matrix representation
#define DEFAULT_DEGREE 5                    ///< default number of pairwise comparisons to use per sequence

#define MAX_NUM_PARALLELIZED_JOBS 100
#define MIN_NUM_PARALLELIZED_JOBS 1

#define MIN_HARDMASK_LENGTH 10

#include <iostream>
#include <fstream>
#include <algorithm>

#include "util/misc.h"
#include "util/opts_list.h"
#include "util/memcheck.h"
#include "seq/sequence.h"
#include "seq/alignment.h"

#include "fsa/anchors.h"
#include "fsa/model.h"

#include "annealing/alignment_DAG.h"
#include "annealing/dotplot.h"
#include "annealing/SparseMatrix.h"
#include "annealing/tree_weights.h"

#include "manager/manager.h"


namespace fsa {

  /**
   * \brief Represent a running FSA program.
   */
  struct FSA {

    // options
    Opts_list opts;          ///< options list

    // output
    bool write_stockholm;    ///< output Stockholm instead of MFA
    bool write_params;       ///< write bubbleplot-format files showing emission parameters during training
    bool write_divergences;  ///< write divergences to file
    bool write_dotplots;     ///< write dotplot-format files showing post. prob. matrices
    bool output_for_gui;     ///< log intermediate alignments formatted for GUI displayer
    std::string gui_prefix;  ///< prefix for gui output file

    // model options
    bool nucprot;            ///< align nucleotide sequence in protein space
    bool is_indel2;          ///< 1 or 2 sets of indel states
    double gap_open1;        ///< gap-open probability of set 1 of indel states
    double gap_extend1;      ///< gap-extend probability of set 1 of indel states
    double gap_open2;        ///< gap-open probability of set 2 of indel states
    double gap_extend2;      ///< gap-extend probability of set 2 of indel states
    bool ragged_ends;        ///< allow for easy insertions and deletions at ends of sequences
    int model;               ///< substitution model (also used as seed for learning)
    double time;             ///< Jukes-Cantor/Tamura-Nei evolutionary time parameter
    double alpha_R;          ///< Tamura-Nei alpha_R parameter
    double alpha_Y;          ///< Tamura-Nei alpha_Y parameter
    double beta;             ///< Tamura-Nei beta parameter
    sstring load_probs_file; ///< load pairwise posterior probabilities from file

    // parameter estimation options
    bool learn_gap;                           ///< learn indel parameters
    bool learn_emit_all;                      ///< learn emit parameters over all sequences
    bool learn_emit_bypair;                   ///< learn emit parameters for each pair
    bool nolearn;                             ///< no learning (redundant but convenient)
    bool regularize;                          ///< regularize learned parameters with Dirichlet distribution specified by model
    double regularization_emission_scale;     ///< scaling factor (per sequence pair) for the Dirichlet pseudocounts for emissions (nuc or aa)
    double regularization_transition_scale;   ///< scaling factor (per sequence pair) for the Dirichlet pseudocounts for transitions (nuc or aa)
    double em_min_inc;                        ///< minimal fractional increase in log-likelihood per round EM
    int em_max_iter;                          ///< max number of rounds of EM
    int min_emit_training_data;               ///< minimum amount of sequence data for training emission probs
    int min_gap_training_data;                ///< minimum amount of sequence data for training indel probs

    // sequence annealing options
    int num_refinement_steps;            ///< number of iterative refinement steps
    bool maxsn;                          ///< maximum sensitivity
    bool use_tgf;                        ///< use tgf or maxstep heuristic for weighting (no command-line control)
    double gap_factor;                   ///< the gap factor
    bool enable_dynamic_weights;         ///< dynamic edge-weight calculation
    double minprob;                      ///< minimum posterior probability to store in sparse matrix
    sstring tree_weights_file;           ///< weights for sequence pairs
    bool require_homology;               ///< require some (potential) homology between all sequences considered

    // alignment speedup options
    int bandwidth;                       ///< banding width
    bool fast_defaults;                  ///< only look at 10 * Erdos-Renyi threshold percent of sequence pairs for alignment
    int num_minimum_spanning_trees;      ///< number of minimum spanning trees to use for alignment
    int num_maximum_spanning_trees;      ///< number of maximum spanning trees to use for alignment
    int num_minimum_spanning_palm_trees; ///< number of minimum spanning palm trees to use for alignment
    int degree;                          ///< number of closest sequences to use for alignment
    int kmer_length;                     ///< length of k-mers to use when determining sequence similarity
    bool refalign;                       ///< minimal alignment (star)
    double fraction_alignment_pairs;     ///< fraction of all (n choose 2) pairs to consider during alignment inference
    int num_alignment_pairs;             ///< total number of all (n choose 2) pairs to consider during alignment inference
    int max_ram;                         ///< maximum RAM (in megabytes) to use for DP

    // anchoring options
    bool anchored;                       ///< use anchor annealing
    bool use_translated;                 ///< perform anchoring in protein space
    int anchor_minlen;                   ///< minimum anchor length
    int anchor_max_join_length;          ///< maximum separation of parallel, adjacent anchor to concatenate
    bool hardmasked;                     ///< sequence is hardmasked
    sstring mercator_constraint_file;    ///< input Mercator constraints
    bool exonerate;                      ///< use exonerate
    int exonerate_minscore;              ///< minimum score for exonerate anchors
    bool softmasked;                     ///< sequence is softmasked

    // parallelization options
    bool parallelizing;                  ///< collect post. prob. matrices and cand. edges simutaneously
    int num_parallelized_jobs;           ///< number of jobs to run simultaneously
    bool noannealing;                    ///< false when the annealing step is not needed
	
    // database options
    sstring db_hostname;                 ///< database server host name 
    sstring db_hostaddr;                 ///< database server host IP address
    sstring db_name;                     ///< database name
    int db_port;                         ///< database server port
    sstring db_user;                     ///< database user name
    sstring db_password;                 ///< database password
    bool read_posteriors_from_db;        ///< do not compute posteriors; instead read them from database
    int db_max_ram;                      ///< maximum RAM (in megabytes) to use when the database mode
    bool write_db;                       ///< write post. prob. matrices and candidate edges to database

    // sequence composition
    bool is_dna;                         ///< is it nucleotide sequence?
    bool is_protein;                     ///< is it amino acid sequence?

    // sequence alphabet
    Alphabet alphabet;                   ///< store alphabet for randomizing degenerate characters, etc.
    std::string alphabet_string;         ///< store alphabet in alphabetical order

    // sequence data
    Sequence_database seq_db;            ///< hold all input sequence data
    Sequence_database seq_db_internal;   ///< hold all processed (nondegenerate, etc.) input sequence data (FSA's internal format)
    Sequence_pairs alignment_seq_pairs;  ///< hold sequence pairs to be used for alignment
    Tree_weights tree_weights;           ///< weights for sequence pairs during anchor and sequence annealing
    size_t num_seqs;                     ///< total number of input sequences
    size_t num_seq_pairs;                ///< total number of possible input sequence pairs

    // parallelization
    int w_num_seq_pairs;                 ///< the number of seq pairs to be considered
    int w_prev_seq_pairs_sum;            ///< the number of seq pairs that have been considered by other workers
    bool w_worker;                       ///< is this run by a worker?
    std::pair<int, int> w_start_seq_pair;     ///< the first pair of sequences to be considered

    // parameters
    std::vector<Params> alignment_params;     ///< hold (trained) parameters for each alignment_seq_pair

    /**
     * \brief Constructor.
     */
    FSA (int argc, char** argv);

    /**
     * \brief Initialize command-line options.
     */
    void init_opts();

    /**
     * \brief Input sequence data.
     *
     * Detects whether nucleotide (DNA or RNA) or amino acid sequence.
     * Also adjusts the default regularization scales as appropriate.
     */
    void input_data();

    /**
     * \brief Construct preset options for DNA, RNA and proteins.
     *
     * Sets up default HMM structure (1 or 2 sets of indel states),
     * learning strategy and memory and speed usage before command-line
     * options (which may override these values) are parsed.
     */
    void set_up_defaults();

    /**
     * \brief Parse options.
     *
     * Provides opportunity for command-line options to override defaults.
     */
    void parse_opts();

    /**
     * \brief Strip hardmasked sequence, randomize degenerate characters, etc.
     */
    void assemble_sequence_data();

    /**
     * \brief Choose sequence pairs to consider for inference and alignment.
     * \see Sequence_pair_selector
     */
    void choose_seq_pairs();

    /**
     * \brief Build multiple alignment.
     */
    void build_multiple_alignment();

    /**
     * \brief Build anchored multiple alignment.
     */
    void build_anchored_multiple_alignment();

    /**
     * \brief Run FSA!
     */
    int run();

    /**
     * \brief Init indel parameters.
     *
     * This is a bit weird: It initializes the values to the defaults,
     * then checks to see if command-line values were passed.  If so, 
     * then it uses those instead.
     */
    void init_indel_params (Params& params);

    /** 
     * \brief Init substitution parameters (both nuc and amino acid).
     *
     * \param params parameters to initialize
     * \param char_freq empirical character frequencies
     * \param model substitution model to use
     * \param time evolutionary time parameter of substitution model
     */
    void init_subst_params (Params& params, const std::vector<double>& char_freq, const unsigned model = MODEL_DEFAULT, const double time = TIME_DEFAULT);

    /**
     * \brief Init pseudocounts, scaled versions of the passed seed distribution.
     *
     * Pseudocounts for a single pairwise alignment.
     * \param pseudo pseodocounts object
     * \param seed seed distribution (scaled to get pseudocounts)
     * \param transition_scale scaling factor for transition pseudocounts
     * \param emission_scale scaling factor for emission pseudocounts
     */
    void init_pseudocounts (Params& pseudo, const Params& seed, const double transition_scale, const double emit_scale);

    /**
     * \brief Interface to Model::train_params_engine.
     *
     * Decides whether to learn gap parameters based on FSA::learn_gap.
     * Doesn't learn if the estimated memory usage is greater than max_ram.
     * \param params Params object to train
     * \param seq_db_train sequence data
     * \param subset sequence pairs to consider
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     * \param ragged_ends 
     * \param pseudocounts pseudocounts
     * \param learn_emit learn emission parameters
     * \return whether training completed
     * \see Model::train_params_engine
     */
    bool train_params (Params& params, const Sequence_database& seq_db_train, const Sequence_pairs& subset,
		       const bool left_match, const bool right_match, const bool ragged_ends,
		       const Params& pseudocounts, bool learn_emit);

    /**
     * \brief Interface to Model::get_alignment_post_probs_engine.
     *
     * \param params parameters
     * \param xseq first sequence
     * \param yseq second sequence
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     * \param ragged_ends 
     * \return computed Post_probs object
     * \see Model::get_pairwise_post_probs_engine
     */
    Post_probs get_pairwise_post_probs (Params& params, const Sequence& xseq, const Sequence& yseq,
					const bool left_match, const bool right_match, const bool ragged_ends);

    /**
     * \brief Interface to Model::get_pairwise_dotplot_engine.
     *
     * \param params parameters
     * \param xseq first sequence
     * \param yseq second sequence
     * \param left_match force homology at start of sequences
     * \param right_match force homology at end of sequences
     * \param ragged_ends 
     * \return computed Dotplot object
     * \see Model::get_pairwise_dotplot_engine
     */
    Dotplot get_pairwise_dotplot (Params& params, const Sequence& xseq, const Sequence& yseq,
				  const bool left_match, const bool right_match, const bool ragged_ends);

    /**
     * \brief Perform pairwise inference.
     *
     * Train parameters and compute posterior alignment probabilities.
     * \see train_params
     * \see get_pairwise_post_probs
     */
    Post_probs perform_pairwise_inference (Params& params, const Sequence& xseq, const Sequence& yseq,
					   const bool left_match, const bool right_match, const bool ragged_ends,
					   const Params &pseudocounts);

    /**
     * \brief Perform anchored pairwise inference.
     *
     * Train parameters and compute posterior alignment probabilities for 
     * all anchored subsequences.
     * \see train_params
     * \see get_pairwise_post_probs
     */
    Post_probs perform_anchored_pairwise_inference (const Params& params_seed, const Sequence& xseq, const Sequence& yseq,
						    const Params& pseudocounts,
						    const Anchors& anchors);

    /**
     * \brief Load pairwise probabilities of alignment from a file.
     *
     * File must be formatted as the .probs files created by FSA.
     * \param sparse_matrices populate this object with SparseMatrix* read from file
     */
    void load_probabilities_from_file (const std::string& filename,
				       std::vector<std::vector<SparseMatrix*> >& sparse_matrices);

    /**
     * \brief Show pairwise distance estimates.
     *
     * Reports a pairwise divergence of -1 if the calculation wasn't done 
     * or the substitution matrix wasn't positive definite.
     */
    void show_divergences (std::ostream& o) const;

    /**
     * \brief Estimates the amount of RAM (in Mb) needed for DP.
     *
     * Size of state space is hard-coded in.
     * \param xlen length of first sequence X
     * \param ylen length of second second Y
     * \return Estimated lower-bound on size of DP matrix in megabytes.
     */
    size_t estimate_ram_needed (size_t xlen, size_t ylen) const;

    void init_for_mw_worker (const std::pair<int, int> start_seq_pair, const int prev_seq_pairs_sum, const int num_seq_pairs); 

    /**
     * \brief Do we store a sequence pair for alignment?
     *
     * Handles case of parallelization.
     */
    bool is_valid_worker_seq_pair (const int cnt, const int i, const int j) const; 

  };

}

#endif /* FSA_FSA_INCLUDED */

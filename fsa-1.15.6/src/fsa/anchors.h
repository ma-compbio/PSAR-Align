
/**
 * \file anchors.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef FSA_ANCHORS_INCLUDED
#define FSA_ANCHORS_INCLUDED

#include <iostream>
#include <fstream>

#include "seq/sequence.h"
#include "seq/alignment.h"
#include "annealing/dotplot.h"
#include "annealing/SparseMatrix.h"
#include "annealing/alignment_DAG.h"
#include "fsa/constraints.h"
#include "fsa/model.h"
#include "fsa/sequence_pair_selector.h"

#define MEGABYTE 1048576

#define ANCHOR_NUC_MINLEN_DEFAULT 10
#define ANCHOR_AA_MINLEN_DEFAULT 7

#define EXONERATE_MINSCORE_DEFAULT 100

namespace fsa {

  /**
   * \brief Represent a single ungapped anchor.
   *
   * Assumes a 0-based coordinate system and a fully-closed interval.
   */
  struct Anchor {

  public:

    Interval xcoords;           ///< interval of the anchor in sequence X
    Interval ycoords;           ///< interval of the anchor in sequence Y
    unsigned x;                 ///< center of the interval xcoords
    unsigned y;                 ///< center of the interval ycoords
    size_t length;              ///< length of the (ungapped) matched sequence

  private:

    double __score;             ///< score of the anchor (higher is better)
    bool __external_scoring;    ///< was the anchor scored according to some external measure?
    bool __immutable;           ///< an anchor which should never be pruned or re-scored

  public:

    /**
     * \brief Constructor.
     */
    Anchor() { }

    /**
     * \brief Constructor.
     */
    Anchor (Interval xcoords, Interval ycoords);

    /**
     * \brief Constructor.
     */
    Anchor (Interval xcoords, Interval ycoords, double score);

    /**
     * \brief Set score.
     *
     * Immutable anchors cannot be rescored.
     * \param score new score
     */
    inline void set_score (const double score) {
      if (!is_immutable())
	__score = score;
    }

    /**
     * \brief Get score.
     */
    inline double get_score() const {
      return __score;
    }

    /**
     * \brief Was the anchor scored according to some external measure?
     *
     * If so, then it shouldn't be updated according to length.
     */
    inline bool is_external_scoring() const {
      return __external_scoring;
    }

    /**
     * \brief Mark the score as external.
     */
    inline void set_external_scoring() {
      __external_scoring = true;
    }

    /**
     * \brief Mark anchor as immutable.
     *
     * Immutable anchors should never be pruned or rescored.
     * Useful primarily for constraints which we want to force
     * (we assume that they are consistent,
     * although nothing will break if they are not).
     * Sets anchor score to 1.0 - DOUBLE_VERY_TINY.
     */
    inline void set_immutable() {
      __score = 1.0 - DOUBLE_VERY_TINY;
      __immutable = true;
    }

    /**
     * \brief Is the anchor score immutable?
     */
    inline bool is_immutable() const {
      return __immutable;
    }

    /**
     * \brief Is it parallel to another anchor?
     */
    inline bool is_parallel (const Anchor& a) const {
      return (xcoords.start - a.xcoords.start == ycoords.start - a.ycoords.start);
    }

    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<< (std::ostream& o, const Anchor& anchor)  {
      o << "[" << anchor.xcoords.start << ", " << anchor.xcoords.end << "] ~ ["
	<< anchor.ycoords.start << ", " << anchor.ycoords.end << "]"
	<< " : " << anchor.x << " ~ " << anchor.y
	<< " => " << std::setprecision (PRECISION_DEFAULT) << anchor.get_score();
      if (anchor.is_external_scoring())
	o << " (external_scoring)";
      if (anchor.is_immutable())
	o << " (immutable)";
      return o;
    }

    /**
     * \brief Re-calculate x, y and length.
     */
    void update();

    /**
     * \brief Map score to a high probability.
     *
     * Linear map from a score in the interval [0,1] to a 
     * high probability in [1 - 0.01, 1].
     * The point of this is: When performing anchor annealing, we want an actual
     * score/probability distribution on candidate anchors.
     * However, post-anchor annealing, the probability distribution over anchored
     * characters is a series of step functions, ie all the mass is concentrated
     * in single aligned pairs.  We therefore no longer want the 'score'
     * but rather a probability near 1, where the mapping preserves
     * the ordering of the scores.
     * \return probability of anchor
     * \see FSA::build_anchored_multiple_alignment
     */
    inline double score_to_prob() const {
      return ((1.0 - DOUBLE_TINY) + DOUBLE_TINY * get_score());
    }

    /**
     * \brief Compute the p-value associated with two anchored subsequences.
     *
     * Computed p-value is that associated with the null model,
     * so a lower p-value is more significant.
     * \param xseq complete sequence X
     * \param yseq complete sequence Y
     * \return p-value of anchored subseqs under a null model
     */
    double p_value (const std::string& xseq, const std::string& yseq) const;

    /**
     * \brief Compute the p-value associated with two anchored subsequences.
     *
     * Computed p-value is that associated with the null model,
     * so a lower p-value is more significant.
     * \param xseq complete sequence X
     * \param yseq complete sequence Y
     * \return p-value of anchored subseqs under a null model
     */
    double p_value (const std::string* xseq, const std::string* yseq) const;

    /**
     * \brief Lexical order, preferring x coordinate, on the left-hand boundaries.
     *
     * Sorted first by the left-hand x coordinate,
     * then the left-hand y coordinate.
     */
    bool operator< (const Anchor& r) const {
      return lexical_comparison_x (*this, r);
    }

    /**
     * \brief Equality.
     */
    bool operator== (const Anchor& r) const {
      if ((xcoords.start == r.xcoords.start) && (xcoords.end == r.xcoords.end) && (ycoords.start == r.ycoords.start) && (ycoords.end == r.ycoords.end))
	return true;
      else
	return false;
    }

    /**
     * \brief Lexical order on x coordinate.
     */
    static bool lexical_comparison_x (const Anchor& l, const Anchor& r) {
      if (l.xcoords.start == r.xcoords.start)
	return (l.xcoords.end < r.xcoords.end);
      else
	return (l.xcoords.start < r.xcoords.start);
    }

    /**
     * \brief Lexical order on y coordinate.
     */
    static bool lexical_comparison_y (const Anchor& l, const Anchor& r) {
      if (l.ycoords.start == r.ycoords.start)
	return (l.ycoords.end < r.ycoords.end);
      else
	return (l.ycoords.start < r.ycoords.start);
    }

    /**
     * \brief Lexical order, preferring x, then y, on the left-hand boundaries.
     */
    static bool lexical_comparison_x_y (const Anchor& l, const Anchor& r) {
      if (l.xcoords.start == r.xcoords.start)
	return (l.ycoords.start < r.ycoords.start);
      else
	return (l.xcoords.start < r.xcoords.start);
    }

    /**
     * \brief Lexical order, preferring y coordinate, then x, on the left-hand boundaries.
     */
    static bool lexical_comparison_y_x (const Anchor& l, const Anchor& r) {
      if (l.ycoords.start == r.ycoords.start)
	return (l.xcoords.start < r.xcoords.start);
      else
	return (l.ycoords.start < r.ycoords.start);
    }

    /**
     * \brief Lexical order, preferring x, then y, on the centroids.
     */
    static bool lexical_comparison_x_y_centroid (const Anchor& l, const Anchor& r) {
      if (l.x == r.x)
	return (l.y < r.y);
      else
	return (l.x < r.x);
    }

  };


  /**
   * \brief A vector of anchors for two sequences.
   */
  struct Anchors {

  public:

    /**
     * \brief Constructor.
     */
    Anchors()
    : xseq (NULL), yseq (NULL) { }

    /**
     * \brief Constructor.
     * \param xseq sequence X
     * \param yseq sequence Y
     */
    Anchors (const std::string& xseq, const std::string& yseq)
    : xseq (&xseq), yseq (&yseq) { }

    /**
     * \brief Constructor.
     * \param xseq sequence X
     * \param yseq sequence Y
     */
    Anchors (const std::string* xseq, const std::string* yseq)
    : xseq (xseq), yseq (yseq) { }

    /**
     * \brief Convert constraints to anchors.
     *
     * For every constraint, creates an anchor for each end of the constraint.
     * If sequence is hardmasked, then constraints whose ends fall in 
     * hardmasked regions are pruned until they don't.
     * Sequences are for hardmasking.
     * \param xseq_orig Sequence for sequence X
     * \param yseq_orig Sequence for sequence Y
     * \param hardmasked sequence is hardmasked
     */
    static Anchors convert_constraints_to_anchors (const Constraints& constraints,
						   const Sequence& xseq_orig, const Sequence& yseq_orig,
						   const bool hardmasked);

    /**
     * Display.
     */
    void show (std::ostream& o) const;

    /**
     * \brief Output formatted for GUI.
     *
     * Format is 
     * (0 ~ 1) == [1,6] ~ [3,8] => 0.9
     * meaning that sequences 0 and 1 have an anchor which
     * spans [1,6] in sequence 0 (X) and [3,8] in sequence 1 (Y)
     * and has a score of 0.9.
     * \param xidx index of sequence X
     * \param yidx index of sequence Y
     */
    void write_gui_output (std::ostream& o, const size_t xidx, const size_t yidx) const;

    /**
     * \brief Convert these anchors to a Post_probs format.
     *
     * Enforces lexical ordering, preferring x coordinate, on centroids of anchors.
     * For subsequent conversion to a SparseMatrix object.
     * Recommended function; memory-efficient.
     */
    const Post_probs convert_to_post_probs() const;

    /**
     * \brief Creates a map from centroid coordinates x and y specified by an anchor to the anchor index in (*this).
     *
     * A set of anchors is nondegenerate if the centroids x and y
     * of anchors uniquely specify the containing anchor.
     * \return false if the set of anchors is degenerate (more than one spanning anchor in (*this)).
     */
    bool create_spanning_map();

    /**
     * \brief Returns the anchor index if an anchor which anchors coordinates x and y exists; returns this->size() otherwise.
     *
     * \param x coordinate anchored in first sequence X
     * \param y coordinate anchored in second sequence Y
     */
    size_t exists_spanning_anchor (const unsigned x, const unsigned y) const;

    /**
     * \brief Returns the anchor which anchors coordinates x and y.
     *
     * Assumes that exists_spanning_anchor (x, y) has been used
     * to check that such an anchor does exist.
     * Calls exists_spanning_anchor (x, y).
     * \param x coordinate anchored in first sequence X
     * \param y coordinate anchored in second sequence Y
     * \return anchor which anchors coordinates x and y
     */
    const Anchor get_spanning_anchor (const unsigned x, const unsigned y);

    /**
     * \brief Resolve parallel degenerate, overlapping and adjacent anchors.
     *
     * Drop degenerate anchors which are entirely spanned by a larger parallel anchor.
     * Merge overlapping parallel anchors.
     * Concatenate nearby (semi-adjacent) parallel anchors which are separated by <= max_join_distance (ungapped).
     * If degenerate anchors are immutable,
     * then the degenerate anchor is kept and the containing anchor discarded
     * (the presumption being that immutable anchors are "special").
     * Calls enforce_immutable_unique and rebuilds the spanning_map.
     * \param max_join_distance maximum separation between parallel anchors to concatenate
     * \return true if non-degenerate spanning map created successfully
     * \see enforce_immutable_unique
     */
    bool resolve_parallel (const size_t max_join_distance = 0, const bool concatenate_immutable = false);

    /**
     * \brief Remove anchors whose centroids overlap with the immutable anchors.
     *
     * This must be called before impose_normalized_probability_distribution
     * in order to make sure that the distribution can be normalized.
     * Does not rebuild the spanning_map.
     */
    void enforce_immutable_unique();

    /**
     * \brief Remove overlaps between anchors.
     *
     * Resolve overlapping anchors favor of the highest-scoring one.
     * Internally enforces lexical ordering of anchors.
     * This is absolutely crucial for efficient DP on unanchored regions.
     * Does not rebuild the spanning_map!
     * \param minlen minimum anchor length
     */
    void remove_overlaps (const size_t minlen = ANCHOR_NUC_MINLEN_DEFAULT);

    /**
     * \brief Remove overlaps in X coordinate between anchors.
     *
     * Assumes that the anchors are already properly sorted.
     * This is the responsibility of the calling function.
     * \param minlen minimum anchor length
     * \param which 0 for X, 1 for Y
     */
    void remove_overlaps (const size_t minlen, const unsigned which);

    /**
     * \brief Update score of an anchor and map to [0,1].
     *
     * Updates the score of the anchor (*this)[idx].
     * Scores are computed as (1 - p_value) of the anchor.
     * They therefore fall in the interval [0,1] and can be interpreted as probabilities, 
     * but do NOT generally form a probability distribution over anchored positions.
     * The scores for a particular sequence position aren't required to be normalized to sum to <= 1.
     * \param idx index of the anchor in (*this)
     * \see impose_normalized_probability_distribution()
     */
    void update_score (const size_t idx);

    /**
     * \brief Convert scores to a normalized probability distribution.
     *
     * Calls update_score() on all anchors.
     * Guarantees that the distribution over anchor centroids is normalized, ie
     *   sum_j P(x_i ~ y_j) <= 1 for every i
     *   sum_i P(x_i ~ y_j) <= 1 for every j
     * where x_i and y_j are centroids of anchored intervals.
     * This normalization is accomplished as follows:
     * If we have candidate anchors A ~ A1 and A ~ A2,
     *   X .....A.......
     *   Y ..A1.....A2..
     * then P(A ~ A1) and P(A ~ A2) as reported by update_score()
     * are scaled by their fractional respective weights.
     * More precisely, an anchor with score P(x_i ~ y_k) becomes
     *    max_j P(x_i ~ y_j) * ( P(x_i ~ y_k) / sum_j P(x_i ~ y_j) )
     * This transformation has the properties:
     * - if there is no overlap, then the score is unchanged
     * - if one anchor has a much higher score than overlapping ones,
     *    then its score is only reduced a little bit, whereas the scores
     *    of the others are drastically reduced
     * - it maps to [0, 1] (assuming that the scores are originally in [0, 1])
     * Because each can be individually interpreted as probabilities (<= 1),
     * this guarantees a normalized distribution.
     * This transformation penalizes non-unique (overlapping) anchors while
     * maintaining their relative ordering.
     * IMPORTANT: The scores of immutable anchors cannot be changed,
     * so this will NOT give a normalized distribution if immutable anchors
     * overlap other anchors.  enforce_immutable_unique must be called first.
     * \see ensure_immutable_unique
     */
    void impose_normalized_probability_distribution();

    /**
     * \brief Number of anchors.
     */
    size_t size() const { return __anchors.size(); }

    /**
     * \brief Store an anchor.
     */
    void store (const Anchor& a) {
      __anchors.push_back (a);
    }

    /**
     * \brief Store anchors.
     */
    void store (const Anchors& new_anchors) {
      __anchors.insert (__anchors.end(),
			new_anchors.begin(), new_anchors.end());
    }

    /**
     * \brief Data access operator.
     */
    const Anchor& operator[] (const size_t i) const {
      assert (i < __anchors.size());
      return __anchors[i];
    }

    /**
     * \brief Data access operator.
     */
    Anchor& operator[] (const size_t i) {
      assert (i < __anchors.size());
      return __anchors[i];
    }

    /**
     * \brief Get reference to last anchor.
     */
    Anchor& back() {
      return __anchors.back();
    }
    
    /**
     * \brief Get iterator to start of __anchors.
     */
    std::vector<Anchor>::iterator begin() {
      return __anchors.begin();
    }

    /**
     * \brief Get iterator to start of __anchors.
     */
    std::vector<Anchor>::const_iterator begin() const {
      return __anchors.begin();
    }

    /**
     * \brief Get iterator to end of __anchors.
     */
    std::vector<Anchor>::iterator end() {
      return __anchors.end();
    }

    /**
     * \brief Get iterator to end of __anchors.
     */
    std::vector<Anchor>::const_iterator end() const {
      return __anchors.end();
    }


  private:

    // pointers are ugly but necessary to have an assignment operator :(
    const std::string* xseq;        ///< first sequence X
    const std::string* yseq;        ///< second sequence Y

    std::vector<Anchor> __anchors;  ///< hold Anchor objects
    /**
     * \brief Map from the anchored coordinates x and y specified by an anchor to the anchor index in (*this).
     *
     * If coordinates x and y are anchored by an anchor A, then we
     * say that anchor A "spans" x and y.
     */
    std::map<std::pair<unsigned, unsigned>, size_t, Util::Duple_less<unsigned, unsigned> > __spanning_map;

  };


  /**
   * \brief Toy anchoring.
   *
   * Use with toy file anchors_cabbage.fasta.
   */
  struct Cabbage_adapter {

  public:

    /**
     * \brief Constructor.
     */
    Cabbage_adapter (const Sequence& xseq, const Sequence& yseq)
    : xseq (xseq), yseq (yseq) { }

    /**
     * \brief Get a set of anchors.
     */
    Anchors get_candidate_anchors() const;

  private:

    /**
     * \brief Named sequences.
     */
    const Sequence& xseq;
    const Sequence& yseq;

  };


  /**
   * \brief Anchoring with MUMmer.
   *
   * Called as 'mummer -mum'.
   */
  struct Mummer_adapter {

  public:

    /**
     * \brief Constructor.
     *
     * \param xseq sequence X
     * \param yseq sequence Y
     */
    Mummer_adapter (const Sequence& xseq, const Sequence& yseq,
		    bool use_translated = false)
    : xseq (xseq), yseq (yseq),
      __use_translated (use_translated) { }

    /**
     * \brief Get a set of anchors.
     *
     * Wrapper function for call_mummer();
     * needed to handle case of using translated anchors.
     * Anchors may be degenerate, overlapping, etc.
     * Does NOT call create_spanning_map() before returning list of anchors.
     * \param minlen minimum length of an exact match
     */
    Anchors get_candidate_anchors (const size_t minlen = ANCHOR_NUC_MINLEN_DEFAULT) const;

  private:

    /**
     * \brief Actually call MUMmer and parse output.
     *
     * Forks a MUMmer process and reads from STDOUT via a pipe.
     * Note that the parameters xseq and yseq are NOT necessarily the same
     * as the member variables xseq and yseq.
     * \param xseq sequence X to run on
     * \param yseq sequence Y to run on
     * \return Anchors found by MUMmer
     */
    Anchors call_mummer (const std::string& xseq, const std::string& yseq, const size_t minlen = ANCHOR_NUC_MINLEN_DEFAULT) const;

    const Sequence& xseq; ///< Sequence X
    const Sequence& yseq; ///< Sequence Y

    /**
     * \brief find anchors in protein space
     *
     * This should only be set if the sequences are nucleotide sequence!
     * Otherwise results will be meaningless.
     */
    const bool __use_translated;

  };

  /**
   * \brief Anchoring with exonerate.
   *
   * Called as '--querytype dna --targettype dna --gappedextension false
   * --saturatethreshold 5 --bigseq true --showvulger false --showalignment false
   * --ryo "%qab %qae %qS %tab %tae %tS %s\n"'.
   */
  struct Exonerate_adapter {

  public:

    /**
     * \brief Constructor.
     *
     * \param xseq sequence X
     * \param yseq sequence Y
     */
    Exonerate_adapter (const Sequence& xseq, const Sequence& yseq,
		       bool use_translated = false, bool softmasked = false)
    : xseq (xseq), yseq (yseq),
      __use_translated (use_translated), __softmasked (softmasked) { }

    /**
     * \brief Get a set of anchors.
     *
     * Wrapper function for call_exonerate.
     * Anchors may be degenerate, overlapping, etc.
     * Does NOT call create_spanning_map before returning list of anchors.
     * \param minscore minimum score of an alignment
     */
    Anchors get_candidate_anchors (const int minscore = EXONERATE_MINSCORE_DEFAULT) const;

  private:

    /**
     * \brief Actually call exonerate and parse output.
     *
     * Forks an exonerate process and reads from STDOUT via a pipe.
     * Note that the parameters xseq and yseq are NOT necessarily the same
     * as the member variables xseq and yseq.
     * \param xseq sequence X to run on
     * \param yseq sequence Y to run on
     * \return alignments found by exonerate
     */
    Anchors call_exonerate (const std::string& xseq, const std::string& yseq,
			    const int minscore = EXONERATE_MINSCORE_DEFAULT) const;

    const Sequence& xseq; ///< Sequence X
    const Sequence& yseq; ///< Sequence Y

    /**
     * \brief find anchors in protein space
     *
     * This should only be called if the sequences are nucleotide sequence!
     * Otherwise results will be meaningless.
     */
    const bool __use_translated;

    /**
     * \brief is the sequence softmasked?
     *
     * exonerate uses softmasking information if present.
     */
    const bool __softmasked;

  };


  /**
   * \brief Adapter class to perform anchor annealing.
   *
   * Note that the sequence information containers
   * hold sequences with the case preserved.
   * This is ok for MUMmer, since it is case-insensitive,
   * and is essential for preserving softmasking for exonerate.
   */
  struct Anchor_resolver {

  public:

    /**
     * \brief Constructor.
     *
     * \param seq_db sequence data
     * \param seq_db_internal sequence data (hardmasked, nondegenerate)
     * \param seq_pairs sequence pairs in seq_db and seq_db_internal to consider
     */
    Anchor_resolver (const Sequence_database& seq_db, const Sequence_database& seq_db_internal,
		     const Sequence_pairs& seq_pairs)
    : seq_db (seq_db), seq_db_internal (seq_db_internal),
      seq_pairs (seq_pairs),
      __constraints_set (Constraints_set (seq_db))
    { }

    /**
     * \brief Add Mercator constraints.
     *
     * \param filename filename of Mercator constraint file
     */
    void add_mercator_constraints (const std::string& filename);

    /**
     * \brief Perform anchor annealing using sequence pair weights to get a set of consistent anchors.
     *
     * Each entry in std::vector<Anchors> corresponds to a set of pairwise anchors;
     * all anchors in the corresponding complete "multiple anchoring" are consistent.
     * Anchor annealing is performed on the "reduced" sequences obtained by
     * storing only the positions corresponding to the centroids of candidate anchors
     * (this dramatically speeds up depth-first search on the sequences).
     * While the same algorithm is used, anchor annealing is qualitatively different
     * from sequence annealing.  The reason is this: Sequence annealing is designed
     * for the case where we have information P (x ~ y) that characters are aligned.
     * In anchor annealing, in contrast, we have p-values of anchors occuring
     * under a null model rather than an actual probability distribution over
     * anchored positions.  We therefore want to greedily call anchors with the smallest
     * p-values and ignore "gap" information entirely (which isn't well-defined
     * for anchors).  We can do exactly this by using the maxstep weighting scheme 
     * with a gap factor of 0 (yes, it's a bit of a hack).
     * Note that we always call Anchors::impose_normalized_probability_distribution
     * before performing anchor annealing in order to heuristically reduce
     * the p-values associated with non-unique anchors.
     * \param tree_weights weights for sequence pairs during anchor annealing
     * \param minlen minimum length of anchor
     * \param max_join_distance maximum separation between parallel anchors to concatenate
     * \param use_translated find anchors in protein space
     * \param use_exonerate call exonerate to find anchors
     * \param softmasked sequence is softmasked (used by exonerate)
     * \param hardmasked sequence is hardmasked (used for Mercator coord mapping)
     * \param num_refinement_steps number of iterative refinement steps
     * \param output_for_gui output anchoring information for MAD GUI
     * \param gui_prefix prefix for GUI output filename
     * \return set of resolved anchors
     * \see Anchors::impose_normalized_probability_distribution()
     */
    const std::vector<Anchors> get_resolved_anchors (const Tree_weights& tree_weights,
						     const size_t minlen = ANCHOR_NUC_MINLEN_DEFAULT, const size_t max_join_distance = 0, bool const use_translated = false,
						     const bool use_exonerate = false, const int minscore = EXONERATE_MINSCORE_DEFAULT, const bool softmasked = false,
						     const bool hardmasked = true,
						     const size_t num_refinement_steps = 0,
						     const bool output_for_gui = false, const std::string gui_prefix = "") const;
    /**
     * \brief Perform anchor annealing on an equidistant star phylogeny to get a set of consistent anchors.
     *
     * All sequence pairs are weighted equally.
     */
    const std::vector<Anchors> get_resolved_anchors (const size_t minlen = ANCHOR_NUC_MINLEN_DEFAULT, const size_t max_join_distance = 0, bool const use_translated = false,
						     const bool use_exonerate = false, const int minscore = EXONERATE_MINSCORE_DEFAULT, const bool softmasked = false,
						     const bool hardmasked = true,
						     const size_t num_refinement_steps = 0,
						     const bool output_for_gui = false, const std::string gui_prefix = "") const;


  private:

    /**
     * \brief Hold all input sequence data.
     */
    const Sequence_database& seq_db;

    /**
     * \brief Hold all input sequence data (hardmasked & nondegenerate).
     */
    const Sequence_database& seq_db_internal;

    /**
     * \brief Hold sequence pairs to find pairwise anchors for.
     * 
     * Indices refer to sequences in seq_db and seq_db_internal
     * (they are indexed identically).
     */
    const Sequence_pairs& seq_pairs;

    /**
     * \brief Constraint information.
     */
    Constraints_set __constraints_set;

  };

}

#endif /* FSA_ANCHORS_INCLUDED */

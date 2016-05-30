
/**
 * \file alignment_DAG.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Ariel Schwartz and Robert Bradley.
 * Sudeep Juvekar wrote the iterative refinement code do_lateral_refinement.
 * Jaeyoung Do wrote the parallelization and database code.
 */



// to do:

// - try converting dfs_f and dfs_b to a true postorder dfs
//    in order to return R_{forward,backward} in topological order
//    => avoid two sorts in reorder() ?



/* Notes */

// I have empirically observed that using the sort algorithm on a vector<Column*>
//  is *much* faster than using list<Column*>::sort.  I have no idea why,
//  but it means that R_forward and R_backward had better be vectors!

// If I decide that I want to keep track of ppv as well as ama (expected accuracy),
//  then add a member to Column and have update_weight_{tgf,maxstep} return a pair<p_match,p_gap>
//  so that I can update those on the fly.


#ifndef ALIGNMENT_DAG_INCLUDED
#define ALIGNMENT_DAG_INCLUDED

#define ACCURACY_ANNOT "Accuracy"
#define GUI_FILENAME_DAG_SUFFIX ".gui"
#define GUI_FILENAME_PROB_SUFFIX ".probs"
#define INVALID_EDGE -1e10 // to do: this could be fragile; replace with negative return value and test for that?

#include <stack>
#include <list>
#include <map>
#include <queue>
#include <iostream>
#include <algorithm>

#include "config.h"
#ifdef HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#endif

#ifdef HAVE_TR1_UNORDERED_SET
#include <tr1/unordered_set>
#endif

#include "util/hash_fcn.h"
#include "seq/sequence.h"
#include "seq/alignment.h"
#include "annealing/SparseMatrix.h"
#include "annealing/tree_weights.h"

#include "manager/manager.h"

namespace fsa {

  // Note:
  // All indexing of sequence data is 0-based here,
  // but SparseMatrix uses 1-based indexing.

  /**
   * \brief (sequence, position) duple
   */
  typedef std::pair<size_t, unsigned> Seq_pos;

  /**
   * \brief Map of sequences to positions.
   *
   * I have experimented with using a hash_map instead,
   * but observed better performance in practice with a standard map.
   */
  typedef std::map<size_t, unsigned> Seq_pos_map;


  /**
   * \brief Column of the alignment.
   *
   * These are the nodes of the DAG.
   */
  class Column {

  private:

    size_t index;               ///< the position of the column in the alignment (0-based)
    size_t orig_index;          ///< the original index (when DAG first initialized)
    Seq_pos_map seq_pos_map;    ///< map of (seq_id,pos) pairs in the column
    Column* merged_into;        ///< the Column which this column has merged into
    bool marked;                ///< indicates if the column has been visited by the dfs 
    bool dead;                  ///< mark a column as dead (dropped from the DAG)

  public:

    static Manager* manager;   ///< manager for database and MW

    /**
     * \brief Constructor.
     *
     * \param index the initial index of the column
     */
    Column (size_t index)
      : index (index), orig_index (index),
      merged_into (this), marked (false), dead (false)
      { }

    /**
     * \brief Compares two columns based on their current indices.
     *
     */
    bool operator< (Column const &r) const {
      return (index < r.index);
    }

    /**
     * \brief Compares two columns based on their current indices.
     *
     */
    bool operator== (Column const &r) const {
      return (index == r.index);
    }

    /**
     * \brief Number of sequences in column.
     *
     */
    inline size_t size() const {
      return seq_pos_map.size();
    }

    /**
     * \brief Has the column has been marked as visited?
     *
     */
    inline bool is_marked() const {
      return marked;
    }

    /**
     * \brief Mark the column as visited.
     *
     */
    inline void mark() {
      marked = true;
    }

    /**
     * \brief Unmark the column as visited.
     *
     */
    inline void unmark() {
      marked = false;
    }

    /**
     * \brief Is the column dead?
     *
     */
    inline bool is_dead() const {
      return dead;
    }

    /**
     * \brief Mark the column as dead.
     *
     * Columns are marked as dead when they have been merged into another column
     * (ie, the node which they represent is no longer present in the DAG).
     */
    inline void set_dead() {
      dead = true;
    }

    /**
     * \brief Returns the current index of the column.
     *
     * The index of the column is its position in the total order specified by the alignment.
     * 0-based.
     */
    inline size_t get_index() const {
      return index;
    }

    /**
     * \brief Set the index.
     *
     */
    inline void set_index (const size_t idx) {
      index = idx;
    }

    /**
     * \brief Returns the original index of the column.
     *
     */
    inline const size_t get_orig_index() const {
      return orig_index;
    }

    /**
     * \brief Get the Column* which this column has been merged into.
     *
     */
    inline Column* get_merged_into() const {
      return merged_into;
    }

    /**
     * \brief Set the Column* which this column has been merged into.
     *
     */
    inline void set_merged_into (Column* col) {
      merged_into = col;
    }

    /**
     * \brief Get the sequence positions represented by this column.
     *
     */
    inline const Seq_pos_map& get_seq_pos_map() const {
      return seq_pos_map;
    }

    /**
     * \brief Does this column contain a particular sequence?
     */
    inline bool contains_seq (const size_t i) const {
      return seq_pos_map.find (i) != seq_pos_map.end();
    }

    /**
     * \brief Add a sequence position.
     *
     */
    inline void add_seq_pos (const Seq_pos& seq_pos) {
      seq_pos_map.insert (seq_pos);
    }

    /**
     * \brief Remove the entry for a particular sequence.
     */
    inline void erase_seq_pos (const size_t seq) {
      if (seq_pos_map.erase (seq) != 1)
	THROWEXPR ("Invalid attempt to erase a (sequence, position) pair.");
    }

    /**
     * \brief Get the expected accuracy of the column.
     *
     * Normalization is by the sum-of-pairs number of characters in the column
     * (but only if we have posterior information available for a pair).
     * Note that the denominator has to be a float due to the weights.
     * \param tree_weights weights for sequence pairs during annealing
     * \return (numerator, denominator)
     */
    std::pair<float, float> get_accuracy_normalized (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, const size_t num_seqs) const;

    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<< (std::ostream& o, const Column& col) {

      o << "column " << col.index << ": ";

      const Seq_pos_map seqs = col.seq_pos_map;
      // catch case of empty (all gaps) column
      if (seqs.begin() == seqs.end()) {
	return o;
      }
      // else display aligned characters
      Seq_pos_map::const_iterator iter = seqs.begin();
      o << "(" << iter->first << ", " << iter->second << ")";
      for (++iter; iter != seqs.end(); ++iter)
	o << " ~ (" << iter->first << ", " << iter->second << ")";

      return o;
    }

    /**
     * \brief Get the expected sum-of-pairs score of this column.
     *
     * Calculates the value of the objective function
     * (accuracy modified by the gap factor) for this column.
     * \param tree_weights weights for sequence pairs during annealing
     * \param gap_factor gap factor
     * @see change_in_expected_score
     */
    float expected_score (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor, size_t num_seqs) const;

    /**
     * \brief Get the change in expected sum-of-pairs score from aligning a character here.
     *
     * Given a character (sequence and position seq_pos) which is not aligned in this column,
     * returns the /change/ in expected value of the objective function (accuracy modified by gap factor) 
     * associated with aligning the character as part of this column.
     * Note that we can use this function to calculate the expected change resulting from
     * removing seq_pos from this column, calculating the expected change, and re-adding,
     * even when seq_pos already is aligned in this column, by setting skip_seq_pos = true.
     * \param seq_pos a character (sequence, position) which isn't aligned in this column
     * \param tree_weights weights for sequence pairs during annealing
     * \param gap_factor gap factor
     * \param skip_seq_pos if false, then return INVALID_EDGE if the character is already present in the column; else ignore
     * \return INVALID_EDGE if there already is a character from the sequence in seq_pos aligned here
     */
    float change_in_expected_score (const Seq_pos& seq_pos, std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor, size_t num_seqs, bool skip_seq_pos = false) const;

  };

  /**
   * \brief Edge between two columns.
   *
   * Store information about a candidate edge (between columns, i.e., nodes of the DAG).
   * Calculates edge weights.
   */
  class Edge {

  public:

    Column* source;              ///< source column (outgoing edge)
    Column* dest;                ///< destination column (incoming edge)
    float weight;                ///< column weight

    unsigned short weight_size;  ///< size of source and dest when the weight was last calculated

    static Manager* manager;     ///< manager for database and MW

    /**
     * \brief Constructor.
     *
     * \param source source column
     * \param dest dest column
     * \param weight weight of the edge
     */
    Edge (Column* source, Column* dest,
	  float weight, unsigned short weight_size)
      : source (source), dest (dest),
      weight (weight), weight_size (weight_size)
    { }

    /**
     * \brief If requested, updates weight based on tgf.
     *
     * Updates change in expected sum-of-pairs score (accuracy modified by gap factor)
     * from adding the edge.
     * Note that while annealing with tgf weights is a true steepest ascent algorithm,
     * the updated weights can increase iff the edge is inconsistent
     * (in this case the steepest ascent condition doesn't apply, since it's only
     * a guarantee for accepted edges).
     * \param tree_weights weights for sequence pairs during annealing
     * \param update_weight update weight of the edge
     * \return true if inconsistent edge
     */
    bool update_weight_tgf (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, bool update_weight = true);

    /**
     * \brief If requested, updates weight based on maxstep.
     *
     * Updates change in expected sum-of-pairs score (accuracy modified by gap factor)
     * from adding the edge.
     * \param tree_weights weights for sequence pairs during annealing
     * \param update_weight update weight of the edge
     * \return true if inconsistent edge
     */
    bool update_weight_maxstep (std::vector<std::vector<SparseMatrix *> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor, bool update_weight = true);

    /**
     * \brief Compare two edges based on their weights.
     */
    bool operator< (Edge const e2) const {
      return weight < e2.weight;
    }

    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<<(std::ostream& o,const Edge& edge) {
      o << *edge.source << " -> " << *edge.dest << "; weight = " << edge.weight;
      return o;
    }

    /**
     * \brief Does the weight need to be updated?
     * \return true if the source or dest column has changed size since the weight was last calculated
     */
    inline bool weight_outdated() const {
      return (weight_size != (source->size() + dest->size()));
    }

    /**
     * \brief Output for GUI.
     *
     * Only for use after calling add_edge();
     * assumes that source and dest of this edge have been reset with get_merged_into().
     */
    inline void write_gui_output (std::ostream& o) {
      assert (source->get_merged_into() == source);
      assert (dest->get_merged_into() == source);
      o << source->get_orig_index() << " -> " << dest->get_orig_index() << endl;
    }

  };


  /**
   * \brief Function object for defining a binary comparison operator for column pointers.
   *
   * Comparison is based on the corresponding columns' indices.
   * @see Column::operator<().
   */
  class greater_index : std::binary_function<Column*, Column*, bool> {
  public:
    bool operator() (Column* x, Column* y) const { return (*y  < *x); }
  };

  /**
   * \brief Function object for defining a binary comparison operator for column pointers.
   *
   * Comparison is based on the corresponding columns' indices.
   * @see Column::operator<().
   */
  class smaller_index : std::binary_function<Column*, Column*, bool> {
  public:
    bool operator() (Column* x, Column* y) const { return (*x  < *y); }
  };

  /**
   * \brief Function object for defining a binary comparison operator for column pointers.
   *
   * Comparison is based on the corresponding edges' weights.
   */
  class smaller_weight : std::binary_function<Edge*, Edge*, bool> {
  public:
    bool operator() (Edge* x, Edge* y) const { return (x->weight < y->weight); }
  };

  /**
   * \brief Function object for hashing pairs of columns.
   *
   * Cast to 64-bit representation (even though may actually be a 32-bit representation
   * on many systems) for safety.
   */
  class column_pair_hash : std::unary_function<std::pair<Column*, Column*>, bit64_t> {
  public:
    bit64_t operator() (std::pair<Column*, Column*> col_pair) const {
      return Hash_functions::bit64_t_pair_hash (reinterpret_cast<bit64_t> (col_pair.first), reinterpret_cast<bit64_t> (col_pair.second));
    }
  };

  /**
   * \brief Map a position in a sequence to the containing Column*.
   */
  typedef std::vector<Column*> Seq_pos_col_map;

  /**
   * \brief Map a sequence position to the containing Column*.
   *
   * Indexed by [sequence][position].
   */
  typedef std::vector<Seq_pos_col_map > Seq_pos_col_maps;

  /**
   * \brief Store a multiple alignment as a DAG and perform sequence annealing.
   *
   * Imposes a total order on the DAG at all times with 
   * an explicit vector representation of the columns.
   */
  class Alignment_DAG {

  private:

    std::vector<Column*> columns;                  ///< the current columns of the alignment
    Seq_pos_col_maps seq_pos_col_maps;             ///< mapping from sequence positions pair<seq, pos> to the containing Column*

    const Sequence_database& seq_db;               ///< sequence data
    size_t num_seqs;                               ///< number of sequences in the alignment
    size_t num_columns;                            ///< number of columns in the alignment

    float expected_accuracy;                       ///< expected accuracy of the alignment

  public:

    /**
     * \brief Constructor.
     *
     * Initialize alignment DAG to the null alignment.
     * Initializes to the "maximal chain decomposition" null alignment,
     * which imposes the total order corresponding to a chain decomposition
     * where the number of chains is maximized.
     * This also corresponds to maximizing the possible number of gap-opens.
     */
    Alignment_DAG (const Sequence_database& seq_db);

    /**
     * \brief Constructor.
     *
     * Initialize alignment DAG to the passed Stockholm alignment.
     */
    Alignment_DAG (const Sequence_database& seq_db, const Stockholm& stock);

    /**
     * \brief Destructor.
     */
    ~Alignment_DAG();

    /**
     * \brief Get the columns of the alignment.
     */
    const std::vector<Column*>& get_columns() const {
      return (*this).columns;
    }

    /**
     * \brief Drop all-gaps columns of the alignment.
     */
    void drop_all_gaps_columns();

    /**
     * \brief Depth-first search for toplogical sorting a la Tarjan.
     *
     * Computes the Dilworth chain decomposition for the graph.
     * Uses a non-recursive (external stack-based) implementation of a postorder traversal of the graph.
     * Pretty tricky.
     */
    void dfs_topological_sort();

    /**
     * \brief Perform sequence annealing with an equidistant star phylogeny.
     *
     * All sequence pairs are weighted equally.
     * \param sparse_matrices pairwise posterior probabilities
     * \param manager database manager
     * \param use_tgf use tgf instead of maxstep weighting
     * \param gap_factor gap factor
     * \param enable_dynamic_weights re-weight edges after merges 
     * \param edge_weight_threshold the threshold for accepting edges
     * \param output_for_gui write GUI-formatted output to disk
     * \param gui_prefix prefix for GUI filename
     */
    void anneal (std::vector<std::vector<SparseMatrix*> >& sparse_matrices,
		 Manager& manager,
		 const bool use_tgf,
		 const float gap_factor, const bool enable_dynamic_weights, const float edge_weight_threshold,
		 const size_t num_refinement_steps,
		 const bool output_for_gui, const std::string gui_prefix);


    /**
     * \brief Perform sequence annealing using sequence pair weights.
     *
     * \param sparse_matrices pairwise posterior probabilities
     * \param tree_weights weights for sequence pairs during annealing
     * \param manager database manager
     * \param use_tgf use tgf instead of maxstep weighting
     * \param gap_factor gap factor
     * \param enable_dynamic_weights re-weight edges after merges 
     * \param edge_weight_threshold the threshold for accepting edges
     * \param output_for_gui write GUI-formatted output to disk
     * \param gui_prefix prefix for GUI filename
     */
    void anneal (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights,
		 Manager& manager,
		 const bool use_tgf,
		 const float gap_factor, const bool enable_dynamic_weights, const float edge_weight_threshold,
		 const size_t num_refinement_steps,
		 const bool output_for_gui, const std::string gui_prefix);

    /**
     * \brief Display the DAG in human-readable format.
     *
     * Show the current total order.
     * If aligned_only, then only show columns with aligned characters.
     */
    void show (std::ostream& o, bool aligned_only = false) const;

    /**
     * \brief Get the Stockholm object corresponding to the current alignment.
     *
     * Mark up with unnormalized column accuracy scores.
     * \param dilworth_resort impose a total order corresponding to a minimal chain decompostion of the DAG
     */
    Stockholm get_stockholm() const;

    /**
     * \brief Get the Stockholm object corresponding to the current alignment.
     *
     * Mark up with normalized column accuracy scores.
     * \param tree_weights weights for sequence pairs during annealing
     */
    Stockholm get_stockholm (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights) const;


  private:

    /**
     * \brief Get the Column* specified by seq_pos.
     */
    inline Column* get_seq_pos_col (const Seq_pos& seq_pos) const {
      assert ((seq_pos.first >= 0) && (seq_pos.first < seq_pos_col_maps.size()));
      assert ((seq_pos.second >= 0) && (seq_pos.second < seq_pos_col_maps[seq_pos.first].size()));
      assert (seq_pos_col_maps[seq_pos.first][seq_pos.second] != 0);

      return seq_pos_col_maps[seq_pos.first][seq_pos.second];
    }

    /**
     * \brief Set the Column* specified by seq_pos.
     */
    inline void set_seq_pos_col (const Seq_pos& seq_pos, Column* col) {
      assert ((seq_pos.first >= 0) && (seq_pos.first < seq_pos_col_maps.size()));
      assert ((seq_pos.second >= 0) && (seq_pos.second < seq_pos_col_maps[seq_pos.first].size()));

      seq_pos_col_maps[seq_pos.first][seq_pos.second] = col;
    }

    /**
     * \brief Output for GUI.
     *
     * Should only be called before annealing begins.
     * \param show_edges show edges of the DAG
     */
    void write_gui_output (std::ostream& o, const bool show_edges = false) const;

    /**
     * \brief Forward depth-first search of the Pearce-Kelly algorithm.
     *
     * Search forwards, ie in the direction from l_bound to upper bound.
     * Look for cycles in the graph, ie a path l_bound ~> u_bound.
     * Store visited nodes (columns) in R_forward.
     * Non-recursive (external stack-based) implementation of DFS.
     * Note that the traversal is neither preorder nor postorder, but rather something arbitrary
     * (because the columns reached by the iteration over sequence positions
     * won't be sorted w.r.t. the total order of the DAG).
     * \param l_bound lower bound of the DFS
     * \param u_bound upper bound of the DFS
     * \param R_forward visited nodes
     */
    bool dfs_f (Column* l_bound, Column* u_bound, std::vector<Column*>& R_forward);

    /**
     * \brief Backward depth-first search of the Pearce-Kelly algorithm.
     *
     * Searches backwards, ie in the direction from u_bound to u_bound.
     * Store visited nodes (columns) in R_backward.
     * Non-recursive (external stack-based) implementation of DFS.
     * Note that the traversal is neither preorder nor postorder, but rather something arbitrary
     * (because the columns reached by the iteration over sequence positions
     * won't be sorted w.r.t. the total order of the DAG).
     * \param l_bound lower bound of the (backward) DFS
     * \param u_bound upper bound of the (backward) DFS
     * \param R_backward visited nodes
     */
    void dfs_b (Column* u_bound, Column* l_bound, std::vector<Column*>& R_backward);

    /**
     * \brief Reorder procedure of the Pearce-Kelly algorithm.
     *
     * \param R_forward nodes reached by a forward DFS
     * \param R_backward nodes reached by a backward DFS
     */
    void reorder (std::vector<Column*>& R_forward, std::vector<Column*>& R_backward);

    /**
     * \brief Unmark all columns in the vector.
     *
     * \param cols columns to unmark
     */
    void unmark (std::vector<Column*>& cols);

    /**
     * \brief Check whether a candidate edge induces a cycle.
     *
     * Performs cycle-checking with a forward DFS from the destination node.
     * Stores nodes reached by the forward DFS in R_forward for
     * later use by add_edge.
     * \param edge edge to add
     * \param l_bound lower bound on DFS
     * \param u_bound upper bound on DFS
     * \param R_forward nodes reached by forward DFS
     * \see dfs_f
     */
    bool induces_cycle (Edge* edge, Column* l_bound, Column* u_bound, std::vector<Column*>& R_forward);

    /**
     * \brief Merge the source and dest columns of the edge into the source column.
     *
     * \param edge edge to merge
     */
    void merge (Edge* edge);

    /**
     * \brief Add a new edge to the alignment DAG.
     *
     * \param edge edge to add
     * \param l_bound lower bound on DFS
     * \param u_bound upper bound on DFS
     * \param R_forward nodes reached by forward DFS
     */
    void add_edge (Edge* edge, Column* l_bound, Column* u_bound, std::vector<Column*>& R_forward);

    /**
     * \brief Perform one iteration of lateral refinement.
     *
     * Iterate through all columns of the alignment, re-aligning characters 
     * according to the per-column change in expected accuracy.
     * \param tree_weights weights for sequence pairs during annealing
     */
    void do_lateral_refinement (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, const float gap_factor);

    /**
     * \brief Get the overall and per-column normalized expected accuracy (formatted as a string for #=GC annotation).
     *
     * The sparse matrices are necessary for calculating the denominator for normalization.
     * \param tree_weights weights for sequence pairs during annealing
     */
    std::pair<float, std::string> get_accuracy_annotation_normalized (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights) const;

    /**
     * \brief Get the expected sum-of-pairs score of the alignment.
     *
     * Calculates the value of the objective function (accuracy modified by the gap factor)
     * over all columns.
     * \param tree_weights weights for sequence pairs during annealing
     */
    float expected_score (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, const float gap_factor) const;

    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<< (std::ostream& o,const Alignment_DAG& dag) {
      for (std::vector<Column*>::const_iterator iter = dag.columns.begin(); iter != dag.columns.end(); ++iter)
	o << **iter << endl;
      return o;
    }

  };

}

#endif /* ALIGNMENT_DAG_INCLUDED */

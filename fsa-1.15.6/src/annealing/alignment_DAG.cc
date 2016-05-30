
/**
 * \file alignment_DAG.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Ariel Schwartz and Robert Bradley.
 * Sudeep Juvekar wrote the iterative refinement code do_lateral_refinement.
 * Jaeyoung Do wrote the parallelization and database code.
 */

#include <limits>

#include "annealing/alignment_DAG.h"
#include "manager/manager.h"

using namespace fsa;


// use a hash map if available; otherwise use a standard map
#ifdef HAVE_TR1_UNORDERED_MAP
typedef std::tr1::unordered_map<std::pair<Column*, Column*>, Edge*, column_pair_hash> Edge_table;
#else
typedef std::map<std::pair<Column*, Column*>, Edge*> Edge_table;
#endif

Manager*   Edge::manager = NULL;
Manager* Column::manager = NULL;

std::pair<float, float> Column::get_accuracy_normalized (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, const size_t num_seqs) const {

  float p_match = 0.;
  float p_gap = 0.;
  float denom = 0.;    // denominator of expected accuracy calculation: number of (sum-of-pairs) characters
  for (Seq_pos_map::const_iterator seq_pos = seq_pos_map.begin(); seq_pos != seq_pos_map.end(); seq_pos++) {
    const size_t i = seq_pos->first;
    const unsigned ii = seq_pos->second;

    for (size_t j = 0; j < num_seqs; j++) {

      if (i == j)
	continue;

      // do we have the (i,j) posterior probability matrix available?
      // if not, go to the next matrix
      if (sparse_matrices[i][j] == 0 && !manager->get_sparse_matrix (sparse_matrices, i, j)) 
	continue;
      const SparseMatrix* ijMatrix = sparse_matrices[i][j];

      // does this column contains sequence j?
      const Seq_pos_map::const_iterator jj = seq_pos_map.find (j);

      // if seq j isn't present in this column, then increment p_gap accordingly
      if (jj == seq_pos_map.end()) {
	p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (0, ii + 1);
	denom += tree_weights (i, j);
      }

      // if it is, then increment p_match
      else if (i < j) {
	p_match += tree_weights (i, j) * ijMatrix->get_match_prob (ii + 1, jj->second + 1);
	denom += tree_weights (i, j) * 2;
      }

    }
  }

  return std::make_pair (2 * p_match + p_gap, denom);
}

float Column::expected_score (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor, size_t num_seqs) const {

  float p_match = 0.;
  float p_gap = 0.;
  for (Seq_pos_map::const_iterator seq_pos = seq_pos_map.begin(); seq_pos != seq_pos_map.end(); seq_pos++) {
    const size_t i = seq_pos->first;
    const unsigned ii = seq_pos->second;

    for (size_t j = 0; j < num_seqs; j++) {

      if (i == j) 
	continue;

      // do we have the (i,j) posterior probability matrix available?
      // if not, go to the next matrix

      if (sparse_matrices[i][j] == 0 && !manager->get_sparse_matrix (sparse_matrices, i, j))
	continue;
      const SparseMatrix* ijMatrix = sparse_matrices[i][j];

      // does this column contains sequence j?
      const Seq_pos_map::const_iterator jj = seq_pos_map.find (j);

      // if seq j isn't present in this column, then increment p_gap accordingly
      if (jj == seq_pos_map.end())
	p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (0, ii + 1);

      // if it is, then increment p_match
      else if (i < j) {
	p_match += tree_weights (i, j) * ijMatrix->get_match_prob (ii + 1, jj->second + 1);
      }
    }
  }

  return (2 * p_match + gap_factor * p_gap);
}

float Column::change_in_expected_score (const Seq_pos& seq_pos, std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor, size_t num_seqs, bool skip_seq_pos /* = false */) const {

  // sequence and position of the character to be aligned
  const size_t i = seq_pos.first;
  const unsigned ii = seq_pos.second;

  // if so requested, make sure that sequence i doesn't already have a character in this column
  if (!skip_seq_pos && contains_seq (i))
    return INVALID_EDGE;

  // if not, then calculate the change in expected accuracy associated with aligning it to this column
  float p_match = 0.;
  float p_gap = 0.;
  for (size_t j = 0; j < num_seqs; j++) {

    if (i == j)
      continue;

    // do we have the (i,j) posterior probability matrix available?
    // if not, go to the next matrix
    if (sparse_matrices[i][j] == 0  && !manager->get_sparse_matrix (sparse_matrices, i, j))
      continue;
    const SparseMatrix* ijMatrix = sparse_matrices[i][j];

    // does this column contains sequence j?
    const Seq_pos_map::const_iterator jj = seq_pos_map.find (j);

    // 2 cases for sequence j:
    // if gapped, then increment p_gap
    if (jj == seq_pos_map.end())
      p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (0, ii + 1);

    // if aligned, then increment p_match and decrement p_gap
    else {
      p_match += tree_weights (i, j) * ijMatrix->get_match_prob (ii + 1, jj->second + 1);
      p_gap -= tree_weights (i, j) * ijMatrix->get_gap_prob (1, jj->second + 1);
    }

  }

  return (2 * p_match + gap_factor * p_gap);
}

bool Edge::update_weight_tgf (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, bool update_weight /* = true */) {

  // NB: It's crucial that at this point the Edge be between
  // the current (after merges) source and dest columns;
  // this must be assured by the calling function.

  float p_match = 0.;
  float p_gap = 0.;

  // for all pairs of positions (i, ii) ~ (j, jj) in each column which will be aligned after the merge
  const Seq_pos_map& source_seq_pos_map = source->get_seq_pos_map();
  const Seq_pos_map& dest_seq_pos_map = dest->get_seq_pos_map();

  // loop over seqs i,j and associated seq positions ii,jj in columns joined by edge
  for (Seq_pos_map::const_iterator source_seq_pos = source_seq_pos_map.begin(); source_seq_pos != source_seq_pos_map.end(); ++source_seq_pos) {
    const size_t i = source_seq_pos->first;       // sequence i
    const unsigned ii = source_seq_pos->second;   // position ii

    for (Seq_pos_map::const_iterator dest_seq_pos = dest_seq_pos_map.begin(); dest_seq_pos != dest_seq_pos_map.end(); ++dest_seq_pos) {
      const size_t j = dest_seq_pos->first;       // sequence j
      const unsigned jj = dest_seq_pos->second;   // position jj

      // simplest consistency check
      if (i == j)
	return true;

      // do we have the (i,j) posterior probability matrix available?
      // if not, go to the next matrix
      if (sparse_matrices[i][j] == 0 && !manager->get_sparse_matrix (sparse_matrices, i, j))
	continue;
      const SparseMatrix* ijMatrix = sparse_matrices[i][j];

      p_match += tree_weights (i, j) * ijMatrix->get_match_prob (ii + 1, jj + 1);
      p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (0, ii + 1); // so p_gap gets incremented as:
      p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (1, jj + 1); // for each (seq,position) pair (i,ii) and (j,jj), 
      //  add P(position ii in seq i is gapped) and P(position jj in seq j is gapped)

      // Here's a long-winded explanation of p_match and p_gap:
      // Consider merging column 1 with n non-gap characters and column 2 with m non-gap characters.
      // P_match should have (n + m choose 2) terms, 
      // one for each pair of aligned characters in the newly-created column.  
      // P_gap should have 2 * (n + m choose 2) terms, 
      // two for each pair of aligned characters in the newly-created column.
      // Because P_gap is actually incremented at each merging operation rather than being calculated anew, 
      // the incremental increase in P_gap after a merge operation has (n * m + m * n) = (2 * n * m terms).

    }
  }

  // update weight if requested
  if (update_weight) {
    if (CTAGGING(3,ANNEALING_VERBOSE)) {
      CL << " updated weight: " << weight << " (p_match = " << p_match << "; p_gap = " << p_gap << ") -> " << 2 * p_match / p_gap << endl;
    }

    // update weight of edge:
    // calculated like eq (4) in the AMAP paper,
    // but gap factor scaled s.t. value of 1 (instead of 1/2) gives accuracy
    weight = 2 * p_match / p_gap;

  }

  // keep track of the size of the nodes
  // (for tracking whether the re-weighting needs to be redone
  // after edge re-insertion onto heap)
  weight_size = source->size() + dest->size();

  return false;
}

bool Edge::update_weight_maxstep (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor, bool update_weight /* = true */) {

  // NB: It's crucial that at this point the Edge be between
  // the current (after merges) source and dest columns;
  // this must be assured by the calling function.

  float p_match = 0.;
  float p_gap = 0.;

  // for all pairs of positions (i, ii) ~ (j, jj) in each column which will be aligned after the merge
  const Seq_pos_map& source_seq_pos_map = source->get_seq_pos_map();
  const Seq_pos_map& dest_seq_pos_map = dest->get_seq_pos_map();

  // for all pairs of positions (i, ii) ~ (j, jj) in each column which will be aligned after the merge
  for (Seq_pos_map::const_iterator source_seq_pos = source_seq_pos_map.begin(); source_seq_pos != source_seq_pos_map.end(); ++source_seq_pos) {
    const size_t i = source_seq_pos->first;
    const unsigned ii = source_seq_pos->second;

    for (Seq_pos_map::const_iterator dest_seq_pos = dest_seq_pos_map.begin(); dest_seq_pos != dest_seq_pos_map.end(); ++dest_seq_pos) {
      const size_t j = dest_seq_pos->first;
      const unsigned jj = dest_seq_pos->second;

      // simplest consistency check
      if (i == j)
	return true;

      // do we have the (i,j) posterior probability matrix available?
      // if not, go to the next matrix
      if (sparse_matrices[i][j] == 0 && !manager->get_sparse_matrix (sparse_matrices, i, j))
	continue;
      const SparseMatrix* ijMatrix = sparse_matrices[i][j];

      p_match += tree_weights (i, j) * ijMatrix->get_match_prob (ii + 1, jj + 1);
      p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (0, ii + 1);
      p_gap += tree_weights (i, j) * ijMatrix->get_gap_prob (1, jj + 1);

    }
  }

  // update weight if requested
  if (update_weight) {
    if (CTAGGING(3,ANNEALING_VERBOSE)) {
      CL << " updated weight: " << weight << " (p_match = " << p_match << "; p_gap = " << p_gap << ") -> " << (2 * p_match - gap_factor * p_gap) << endl;
    }

    // update weight of edge
    weight = (2 * p_match - gap_factor * p_gap) / (source_seq_pos_map.size() * dest_seq_pos_map.size());

  }

  // keep track of the size of the nodes
  // (for tracking whether the re-weighting needs to be redone
  // after edge re-insertion onto heap)
  weight_size = source->size() + dest->size();

  return false;

}

Alignment_DAG::Alignment_DAG (const Sequence_database& seq_db)
  : seq_db (seq_db),
    num_seqs (seq_db.size()), num_columns (0),
    expected_accuracy (0) {

  // get maximum sequence length
  size_t maxlen = 0;
  for (size_t seq = 0; seq < seq_db.size(); ++seq) {
    const Sequence& sequence = seq_db.get_seq (seq);
    maxlen = std::max (maxlen, sequence.length());
  }

  // initialize seq_pos_col_maps for sequences
  // use 0-based indexing
  seq_pos_col_maps.resize (num_seqs);
  for (size_t r = 0; r < num_seqs; ++r)
    seq_pos_col_maps[r].resize (seq_db.get_seq (r).length(), NULL);

  // store sequence data
  num_columns = 0;
  for (unsigned pos = 0; pos < maxlen; ++pos) {        // position within the sequences
    for (size_t seq = 0; seq < num_seqs; ++seq) {

      // are we done with this sequence?
      const Sequence& sequence = seq_db.get_seq (seq);
      if (pos >= sequence.length())
	continue;

      // store character
      columns.push_back (new Column (num_columns++));  // 0-based indexing for Column::index
      const Seq_pos seq_pos (seq, pos);
      (*(--columns.end()))->add_seq_pos (seq_pos);
      set_seq_pos_col (seq_pos, *(--columns.end()));
    }
  }

}

Alignment_DAG::Alignment_DAG (const Sequence_database& seq_db, const Stockholm& stock)
  : seq_db (seq_db),
    num_seqs (seq_db.size()), num_columns (stock.columns()),
    expected_accuracy (0) {

  // initialize seq_pos_col_maps for sequences
  // use 0-based indexing
  seq_pos_col_maps.resize (num_seqs);
  for (size_t seq = 0; seq < num_seqs; ++seq)
    seq_pos_col_maps[seq].resize (seq_db.get_seq (seq).length(), NULL);

  // create columns
  for (size_t c = 0; c < num_columns; ++c)
    columns.push_back (new Column (c));  // 0-based indexing for Column::index

  // store sequence data
  for (size_t seq = 0; seq < num_seqs; ++seq) {

    // reset position within current sequence
    unsigned pos = 0;

    // iterate over columns, storing (seq, pos) information as we go
    for (std::vector<Column*>::iterator col = columns.begin(); col != columns.end(); ++col) {

      // if seq has a character in column col
      if (!stock.is_gapped (seq, (*col)->get_index())) {
	const Seq_pos seq_pos (seq, pos++);
	(*col)->add_seq_pos (seq_pos);
	set_seq_pos_col (seq_pos, *col);
      }

      // else seq must be gapped in column col, so store nothing

    }

  }

  // drop all-gaps columns (if they were present in the input Stockholm alignment)
  drop_all_gaps_columns();

}

Alignment_DAG::~Alignment_DAG() {

  // columns
  for (std::vector<Column*>::iterator col = columns.begin(); col != columns.end(); ++col)
    (*col)->set_merged_into (NULL);
  for (std::vector<Column*>::iterator col = columns.begin(); col != columns.end(); ++col) {
    Column* colPtr = *col;
    *col = NULL;
    delete colPtr;
  }
  columns.clear();

  Edge::manager   = NULL;
  Column::manager = NULL;

}

void Alignment_DAG::drop_all_gaps_columns() {

  // drop all-gaps columns
  std::vector<Column*> sans_allgaps;
  size_t idx = 0;
  for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {
    if ((*col)->size()) {
      (*col)->set_index (idx++);
      sans_allgaps.push_back (*col);
    }
  }

  columns.assign (sans_allgaps.begin(), sans_allgaps.end());
  num_columns = columns.size();

}

void Alignment_DAG::dfs_topological_sort() {

  CTAG(7,ANNEALING ANNEALING_VERBOSE DAG) << "Finding the most-parsimonious indel structure." << endl;

  // unmark all nodes
  unmark (columns);

  // initialize stack of nodes by pushing on columns holding the first characters of all sequences
  std::stack<Column*> nodes;
  std::vector<bool> on_stack (columns.size(), false);
  std::vector<int> finishing_times (columns.size(), -1);
  std::vector<Column*> starting_cols;
  for (size_t s = 0; s < num_seqs; ++s) {
    if (seq_db.get_seq (s).length() == 0)            // handle case of 0-length sequences
      continue;
    Column* node = get_seq_pos_col (Seq_pos (s, 0)); // remember 0-based indexing
    assert (!node->is_dead());
    if (!on_stack[node->get_index()]) {
      starting_cols.push_back (node);
      on_stack[node->get_index()] = true;
    }
  }
  // now sort these starting columns (MUST preserve the partial order!)
  std::sort (starting_cols.begin(), starting_cols.end(), smaller_index());
  for (std::vector<Column*>::iterator col = starting_cols.begin(); col != starting_cols.end(); ++col)
    nodes.push (*col);

  // perform a postorder traversal to get finishing times for columns
  unsigned time = 0;
  while (!nodes.empty()) {

    // get next node on stack
    Column* node = nodes.top();
    nodes.pop();
    assert (!node->is_dead());

    // if the node which we've just popped is marked, then we're finished with it
    if (node->is_marked()) {
      finishing_times[node->get_index()] = time++;
      continue;
    }

    // else mark the node and push it on the stack
    node->mark();
    nodes.push (node);

    // continue DFS: explore node's neighbors
    std::vector<Column*> neighbors;
    const Seq_pos_map& seq_pos_map = node->get_seq_pos_map();
    for (Seq_pos_map::const_iterator seq_pos_iter = seq_pos_map.begin(); seq_pos_iter != seq_pos_map.end(); ++seq_pos_iter) {

      // get the next position in this sequence
      const Seq_pos next_seq_pos = Seq_pos (seq_pos_iter->first, seq_pos_iter->second + 1);

      // have we reached the end of this sequence?
      if (next_seq_pos.second >= seq_db.get_seq (next_seq_pos.first).length())
	continue;

      // if not, then continue DFS on neighbors
      Column* w = get_seq_pos_col (next_seq_pos);

      // store the neighbors (if they're not already on the stack)
      if (!on_stack[w->get_index()]) {
	neighbors.push_back (w);
	on_stack[w->get_index()] = true;
      }

    }

    // now a *critical* step: use the current total order
    // to impose the correct partial order on the neighboring nodes
    // (otherwise we won't be doing a postorder traversal => things will break mysteriously)
    std::sort (neighbors.begin(), neighbors.end(), smaller_index()); // not doing this will MESS EVERYTHING UP!!

    // then store
    for (std::vector<Column*>::iterator column = neighbors.begin(); column != neighbors.end(); ++column)
      nodes.push (*column);

  } // end DFS loop

  // construct list of live columns for sorting
  std::vector<Column*> columns_tmp;
  for (std::vector<Column*>::iterator column = columns.begin(); column != columns.end(); ++column) {
    if (!(*column)->is_dead()) {
      columns_tmp.push_back (*column);
      assert (finishing_times[(*column)->get_index()] >= 0);        // assert node finished in postorder traversal
      (*column)->set_index (finishing_times[(*column)->get_index()]); // store finishing times as indices
    }
  }
  assert (num_columns == columns_tmp.size());

  // now (reverse) sort according to finishing times
  std::sort (columns_tmp.begin(), columns_tmp.end(), greater_index());

  // store and make column indices sane as we go
  columns.clear();
  for (size_t idx = 0; idx < columns_tmp.size(); ++idx) {
    columns.push_back (columns_tmp[idx]);
    columns[idx]->set_index (idx);
  }

}

Stockholm Alignment_DAG::get_stockholm() const {

  // store alignment as matrix of booleans,
  // with each row representing the alignment for a particular sequence
  // (see Alignment_row::Row_path)
  std::vector<Alignment_row::Row_path*> align_matrix (num_seqs, NULL);
  for (std::vector<Alignment_row::Row_path*>::iterator row_path = align_matrix.begin(); row_path != align_matrix.end(); ++row_path)
    *row_path = new Alignment_row::Row_path (num_columns, false);

  for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {

    // skip dead columns
    if ((*col)->is_dead())
      continue;

    for (size_t seq = 0; seq < num_seqs; seq++) {
      if ((*col)->contains_seq (seq))
	(*align_matrix[seq])[(*col)->get_index()] = true;
    }

  }

  // now add the sequence data
  Stockholm stock (const_cast<Sequence_database&> (seq_db)); // hacky cast away const
  for (size_t r = 0; r < num_seqs; ++r)
    stock.set_row (seq_db.get_seq (r).name, align_matrix[r]);

  stock.assert_flush();

  return stock;
}

Stockholm Alignment_DAG::get_stockholm (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights) const {

  // get the base alignment
  Stockholm stock = get_stockholm();

  // clear out the unnormalized accuracy annotations and 
  // store the properly normalized ones
  stock.clear_annot();
  std::pair<float, std::string> accuracy_pair = get_accuracy_annotation_normalized (sparse_matrices, tree_weights);
  std::string acc = Util::to_string (accuracy_pair.first);
  stock.add_gf_annot (ACCURACY_ANNOT, acc);
  stock.set_gc_annot (ACCURACY_ANNOT, accuracy_pair.second);

  stock.assert_flush();

  return stock;
}

void Alignment_DAG::anneal (std::vector<std::vector<SparseMatrix*> >& sparse_matrices,
			    Manager& manager,
			    const bool use_tgf,
			    const float gap_factor, const bool enable_dynamic_weights, const float edge_weight_threshold,
			    const size_t num_refinement_steps,
			    const bool output_for_gui, const std::string gui_prefix) {

  // initialize dummy Tree_weights (returns all weights as 1.0)
  Tree_weights tree_weights;
  anneal (sparse_matrices, tree_weights,
	  manager,
	  use_tgf,
	  gap_factor, enable_dynamic_weights, edge_weight_threshold,
	  num_refinement_steps,
	  output_for_gui, gui_prefix);

}

void Alignment_DAG::anneal (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights,
			    Manager& manager,
			    const bool use_tgf,
			    const float gap_factor, const bool enable_dynamic_weights, const float edge_weight_threshold,
			    const size_t num_refinement_steps,
			    const bool output_for_gui, const std::string gui_prefix) {

  // set db manager to Column class and Edge class
  Edge::manager   = &manager;
  Column::manager = &manager;

  // log
  CTAG(6,ANNEALING DAG ANNEALING_VERBOSE) << "Creating candidate edge list." << endl;

  // pre-count the number of sequence pairs which we have SparseMatrix entries for
  // as well as the maximum possible number of edges
  size_t seq_pairs = 0;
  size_t num_max_edges = 0;
  for (size_t i = 0; i < num_seqs; ++i)
    for (size_t j = i + 1; j < num_seqs; ++j) {
      if ((sparse_matrices[i][j] != 0) || (sparse_matrices[j][i] != 0 || manager.is_sparse_matrix_available (i, j))) {
	++seq_pairs;
	num_max_edges += std::min (seq_db.get_seq (i).length(), seq_db.get_seq (j).length());
      }
    }

  // gui output
  const std::string gui_filename = gui_prefix + GUI_FILENAME_DAG_SUFFIX;
  std::ofstream gui_file;
  // show initial DAG
  if (output_for_gui) {
    gui_file.open (gui_filename.c_str());
    if (!gui_file.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << gui_filename << "'.");
    write_gui_output (gui_file);
    gui_file << "; Merges (added edges)" << endl
	     << "; Format is:" << endl
	     << ";   source_column -> dest_column => accuracy_change" << endl
	     << "; (new accuracy is source_column_accuracy + dest_column_accuracy + accuracy_change)" << endl
	     << "; (note that all accuracies are unnormalized)" << endl
	     << endl;
  }

  // gui output: post prob matrices
  if (output_for_gui) {
    const std::string prob_filename = gui_prefix + GUI_FILENAME_PROB_SUFFIX;
    std::ofstream prob_file;
    prob_file.open (prob_filename.c_str());
    if (!prob_file.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << prob_filename << "'.");
    for (size_t i = 0; i < num_seqs; ++i)
      for (size_t j = i + 1; j < num_seqs; ++j) {
	if (sparse_matrices[i][j] != 0)
	  sparse_matrices[i][j]->write_gui_output (prob_file);
	else if (sparse_matrices[j][i] != 0)
	  sparse_matrices[j][i]->write_gui_output (prob_file);
	else if (manager.get_sparse_matrix (sparse_matrices, i, j))
	  sparse_matrices[i][j]->write_gui_output (prob_file);

      }
    prob_file.close();
    CTAG(9,ANNEALING ANNEALING_VERBOSE) << "Created GUI output file '" << prob_filename << "'." << endl;
  }

  // create priority queue of candidate edges
  std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> edges;

  // if relevant, get edges from database
  if (manager.is_edges_available()) {

#if defined(HAVE_CONDOR) || defined(HAVE_POSTGRES)

    // log
    CTAG(8,ANNEALING DAG ANNEALING_VERBOSE) << "Getting candidate edges from database." << endl;
    manager.get_edges (edges, seq_pos_col_maps);	
    // log
    CTAG(8,ANNEALING DAG ANNEALING_VERBOSE) << "Got all candidate edges." << endl;

#endif

  }

  else {
    // loop through upper-triangle of sequences:
    //  calculate weights (maxstep or tgf) for each potential edge;
    //  store edges whose weights meet edge_weight_threshold and
    //  have a weight >= the gap_factor (if using tgf)
    CTAG(9,ANNEALING DAG ANNEALING_VERBOSEs) << "Assembling list of candidate edges." << endl;
    size_t seq_pairs_processed = 0;

    for (size_t i = 0; i < num_seqs; ++i) {
      const size_t xlen = seq_db.get_seq (i).length();

      for (size_t j = i + 1; j < num_seqs; ++j) {

	const SparseMatrix* ijMatrix = sparse_matrices[i][j];
	// do we have the (i,j) posterior probability matrix available?
	// if not, go to the next matrix
	if (ijMatrix == 0)
	  continue;

	// for all entries in the SparseMatrix, calculate the corresponding weight
	// and store the edge
	for (unsigned ii = 0; ii < xlen; ii++) { // note 0-based indexing for sequences
	  const float p_gap_ii = ijMatrix->get_gap_prob (0, ii + 1);

	  for (std::vector<Matrix_entry>::iterator rowPtr = ijMatrix->GetRowPtr (ii + 1),
		 rowEnd = rowPtr + ijMatrix->GetRowSize (ii + 1); rowPtr != rowEnd; rowPtr++) {
	    const unsigned jj = rowPtr->first - 1; // convert from SparseMatrix's 1-based coords to the 0-based coords which we use here

	    const float p_match = tree_weights (i, j) * rowPtr->second;
	    if (!p_match)
	      continue;
	    const float p_gap = tree_weights (i, j) * (p_gap_ii + ijMatrix->get_gap_prob (1, jj + 1));
	    const float weight = use_tgf ? (2 * p_match / p_gap) : (2 * p_match - gap_factor * p_gap);

	    // if the edge already scores too low to ever be accepted, then skip it
	    if ((weight < edge_weight_threshold) || (use_tgf && weight < gap_factor))
	      continue;

	    // else create and store it
	    // magic number 2 = 1 + 1 = size of columns

	    Edge* edge = new Edge (get_seq_pos_col (Seq_pos (i, ii)), get_seq_pos_col (Seq_pos (j, jj)),
				   weight, 2);
	    edges.push (edge);

	    if (CTAGGING(-1,ANNEALING_VERBOSE)) {
	      CL << "Stored candidate edge: " << *edge << "." << endl;
	      CL << "  p_match = " << p_match << "; p_gap = " << p_gap << endl;
	    }
	  }
	}

	// log progress
	++seq_pairs_processed;
	const unsigned percent_done = static_cast<unsigned> (std::floor ((100.0 * seq_pairs_processed / seq_pairs) + 0.5));
	if (CTAGGING(7,ANNEALING DAG ANNEALING_VERBOSE) && (seq_pairs_processed % 100 == 0)) {
	  CTAG(7,ANNEALING DAG ANNEALING_VERBOSE) << "Stored candidate edges from " << seq_pairs_processed << "/" << seq_pairs << " sequence pairs (" << percent_done << "% done)." << endl;
	}

      }
    }

  }

  // assert a greedy algorithm
  //#ifndef NDEBUG  
  //  float prev_weight = std::numeric_limits<float>::max();
  //#endif
  // see comment near endif directive

  // keep track of candidate edges:
  // - edges which were previously re-inserted in the heap
  //   (in order to prevent unnecessary re-weighting of duplicate edges)
  // - edges which were found to be inconsistent
  //   (for these, store a NULL value)
  Edge_table edge_table;

  // anneal:
  // loop through the priority queue of edges,
  // greedily adding edges based on their (updated) weights
  CTAG(9,ANNEALING DAG ANNEALING_VERBOSE) << "Annealing by adding edges to the DAG." << endl;

  size_t edges_annealed = 0; // for logging progress
  const size_t num_orig_edges = manager.get_edges_size (edges);
  size_t sequential_skips = 0; // number of sequential edges skipped because they created a cycle
  while (!edges.empty()) {

    // get candidate edge
    Edge *edge = manager.get_next_top_edge (edges, seq_pos_col_maps);

    // log
    if (CTAGGING(3,ANNEALING_VERBOSE)) {
      CL << "Popped edge: " << *edge << endl;
    }

    // lookup current columns for edge:
    // note that this can be linear in the number of merges,
    // although this shouldn't be too bad for common alignments
    while (edge->source != edge->source->get_merged_into())
      edge->source = edge->source->get_merged_into();
    while (edge->dest != edge->dest->get_merged_into())
      edge->dest = edge->dest->get_merged_into();

    // fast check for redundant edges
    // (edges specifying columns which have already been constructed)
    // (wrapped in braces to prevent compiler complaints about
    // 'jump to label skipping initialization of variables' blah blah blah...)
    {
      if (edge->source == edge->dest) {
	if (CTAGGING(3,ANNEALING_VERBOSE)) {
	  CL << "Inconsistent (redundant) edge: " << *edge << endl;
	}
	goto DELETE_EDGE;
      }
    }

    // now DFS for full-blown consistency-checking
    // this is also where the actual Pearce-Kelly takes place
    // (assuming the edge is consistent)
    {

      // check whether this is a duplicate edge:
      // do we already have an edge corresponding to merging these columns
      // in our table of candidate edges edge_table?
      Edge_table::iterator edge_in_table_iter = edge_table.find (std::make_pair (edge->source, edge->dest));

      // use symmetry: edge can be in either direction
      if (edge_in_table_iter == edge_table.end())
	edge_in_table_iter = edge_table.find (std::make_pair (edge->dest, edge->source));

      // have we seen the current edge previously?
      bool in_table = false;
      Edge* edge_in_table = NULL;
      if (edge_in_table_iter != edge_table.end()) {
	in_table = true;
	edge_in_table = edge_in_table_iter->second;
      }

      // is this a previously-seen inconsistent edge?
      if (in_table && (edge_in_table == 0)) { // test for NULL
	if (CTAGGING(3,ANNEALING_VERBOSE)) {
	  CL << "Inconsistent (previously-seen) edge: " << *edge << endl;
	}
	goto DELETE_EDGE;
      }

      // is this a duplicate edge?
      // if we already have a different edge corresponding to this merge
      // (with an updated, "current" weight),
      // then this is a duplicate edge, so get rid of it
      // (it conveys no new information)
      if (in_table && (edge_in_table != edge)) {
	assert (edge_in_table_iter != edge_table.end());
	assert (edge_in_table != 0);

	// if the weight is current, meaning that the columns
	// joined by the edges have not been affected by other merges
	// since the weight was last calculated, then the weight
	// for edge_in_table is correct and we can just delete
	// the current edge
	if (!edge_in_table->weight_outdated()) {
	  if (CTAGGING(3,ANNEALING_VERBOSE)) {
	    CL << "Ignoring edge: " << *edge << " (" << edge << " != " << edge_in_table << ")" << endl;
	  }
	  goto DELETE_EDGE;
	}

	// if the weight isn't current, meaning that the columns
	// joined by the edges have been affected by other merges
	// since the weight was last calculated, then we need to re-weight it
	// in order to ensure that we maintain greediness
	// (otherwise it won't be a steepest ascent algorithm)
	// use a quasi-hack to accomplish this:
	// replace the stored edge with the current one;
	else {
	  edge_table.erase (edge_in_table_iter);
	  edge_in_table_iter = (edge_table.insert (std::make_pair (std::make_pair (edge->source, edge->dest), edge))).first; // get new iterator
	  if (CTAGGING(3,ANNEALING_VERBOSE)) {
	    CL << "Erased and re-inserted edge to force re-weighting." << endl;
	  }
	}

      }

      // recalculate edge weight if using dynamic edge weights
      // as a byproduct, calculate the change in expected accuracy from adding the edge
      if (edge->weight_outdated() && enable_dynamic_weights) {

	// the second argument enable_dynamic_weights is now obsolete;
	// originally the update_weight_{tgf,maxstep} functions were used
	// to calculate changes to the expected accuracy after each merge,
	// but this functionality is retired for speed when 
	// dynamic weighting is disabled
	bool conflict = false;
	if (use_tgf) {
	  conflict = (edge->update_weight_tgf) (sparse_matrices, tree_weights, enable_dynamic_weights);
	}
	else {
	  conflict = (edge->update_weight_maxstep) (sparse_matrices, tree_weights, gap_factor, enable_dynamic_weights);
	}

	// if an inconsistency detected, record it and delete edge
	if (conflict) {

	  // log
	  if (CTAGGING(3,ANNEALING_VERBOSE)) {
	    CL << "Inconsistent edge: " << *edge << endl;
	  }

	  // if present, mark record as inconsistent
	  if (in_table)
	    edge_in_table_iter->second = (Edge*) NULL;
	  // else record it as inconsistent
	  else
	    edge_table.insert (std::make_pair (std::make_pair (edge->source, edge->dest), (Edge*) NULL));

	  goto DELETE_EDGE;
	}

	// if priority queue sorting is wrong, then push and go to the next
	if (edge->weight < edges.top()->weight) {

	  // if this is the first time we've seen this edge, record it
	  if (!in_table) {
	    if (CTAGGING(3,ANNEALING_VERBOSE)) {
	      CL << "Recorded new edge: " << *edge << endl;
	    }
	    edge_table.insert (std::make_pair (std::make_pair (edge->source, edge->dest), edge));
	  }

	  // log
	  if (CTAGGING(3,ANNEALING_VERBOSE)) {
	    CL << "Re-inserted edge: " << *edge << endl;
	  }

	  // push back on heap
	  manager.push_edge (edges, edge);
	  
	  continue;
	}

      } // 	if (edge->weight_outdated() && enable_dynamic_weights)

      // if best edge has weight < edge_weight_threshold or has
      // a weight < gap_factor (if using tgf), then we're done with sequence annealing
      // only relevant for dynamic weights;
      // edges whose initial weights don't satisfy these criteria are never
      // stored on the heap as candidate edges
      if (enable_dynamic_weights
	  && ((edge->weight < edge_weight_threshold) || (use_tgf && (edge->weight < gap_factor)))) {
	while (!edges.empty()) {
	  edge = edges.top();
	  edges.pop();
	  delete edge;
	}
	break;
      }

      // now actually perform the consistency check and Pearce-Kelly

      // NB: Do NOT NOT NOT forget to unmark nodes after visiting!
      // Crucial and easy to forget.
      std::vector<Column*> R_forward;

      // get the lower and upper bounds of the depth-first searches
      // (we have already imposed a total order on columns)
      Column *l_bound, *u_bound;
      if (*(edge->source) < *(edge->dest)) {
	l_bound = edge->source;
	u_bound = edge->dest;
      } else {
	l_bound = edge->dest;
	u_bound = edge->source;
      }

      // if edge is inconsistent (induces cycle)
      if (induces_cycle (edge, l_bound, u_bound, R_forward)) {

	// log
	if (CTAGGING(3,ANNEALING_VERBOSE)) {
	  CL << "Inconsistent edge: " << *edge << endl;
	}

	// erase record if present
	if (in_table)
	  edge_in_table_iter->second = (Edge*) NULL;
	// else record it as inconsistent
	else
	  edge_table.insert (std::make_pair (std::make_pair (edge->source, edge->dest), (Edge*) NULL));

	goto DELETE_EDGE;
      }

      // if edge is consistent (doesn't induce cycle)
      else {

	// add edge with the Pearce-Kelly algorithm
	add_edge (edge, l_bound, u_bound, R_forward);

	// log
	if (CTAGGING(4,ANNEALING DAG ANNEALING_VERBOSE)) {
	  CL << "Added edge: " << *edge << endl;
	}

	// we've added the edge, so remove the record
	if (in_table) {
	  edge_table.erase (edge_in_table_iter);
	  if (CTAGGING(3,ANNEALING_VERBOSE)) {
	    CL << "Erased added edge: " << *edge_in_table << endl;
	  }
	}

#ifndef NDEBUG
	//	if (use_tgf)
	//	  assert (edge->weight - prev_weight < DOUBLE_TINY);
	//	prev_weight = edge->weight;
	// comment in to check that sequence annealing is truly steepest ascent
	// note that this check can fail due to floating-point error if the weights are large (ie the max 10000)
	//  -- RKB 10/21/08
#endif

	++edges_annealed;     // increment counter
	sequential_skips = 0; // reset counter

	//log
	if (CTAGGING(-1,DAG)) {
	  CL << "DAG after annealing " << edges_annealed << "/" << num_orig_edges << " candidate edges." << endl;
	  this->show (CL, true);
	}

	// write gui output if requested
	if (output_for_gui) {
	  edge->write_gui_output (gui_file);
	}

      } // end if edge doesn't induce cycle

    } // end block for DFS



    // clean up: delete edge and log progress
    // goto is ugly but appropriate here...
  DELETE_EDGE:

    delete edge;

    // check whether we've skipped more than 10 * (maximum possible number of remaining edges to anneal) edges;
    // if so, then stop annealing
    if (sequential_skips >= 10 * (num_max_edges - edges_annealed)) {
      CTAG(7,ANNEALING ANNEALING_VERBOSE) << "Skipped more than 10 * (maximum possible number of remaining edges to anneal) sequential edges; breaking." << endl;
      while (!edges.empty()) {
	edge = edges.top();
	edges.pop();
	delete edge;
      }
      break;
    }

    const size_t edges_size = manager.get_edges_size (edges);

    // log progress
    if (CTAGGING(7,ANNEALING ANNEALING_VERBOSE) && (((edges_annealed % 10000) == 0) || (((num_orig_edges - edges_size) % 10000) == 0))) {
      const unsigned percent_done = static_cast<unsigned> (std::floor ((100.0 * edges_annealed / num_orig_edges) + 0.5));
      const unsigned percent_remaining = static_cast<unsigned> (std::floor ((100.0 * edges_size / num_orig_edges) + 0.5));
      CTAG(7,ANNEALING ANNEALING_VERBOSE) << "Annealed " << edges_annealed << "/" << num_orig_edges << " (" << percent_done << "%) candidate edges;"
					  << " " << edges_size << "/" << num_orig_edges << " (" << percent_remaining << "%) remaining candidate edges to process." << endl;
    }

  } // while (!edges.empty())

  // we're done with annealing!
  CTAG(8,ANNEALING DAG ANNEALING_VERBOSE) << "Finished adding edges to the DAG." << endl;

  // perform iterative refinement if requested
  // use "lateral shift"
  if (num_refinement_steps > 0) {

    CTAG(9,ANNEALING REFINEMENT) << "Beginning iterative refinement." << endl;

    // do hill-climbing refinement
    float curr_score = expected_score (sparse_matrices, tree_weights, gap_factor);
    size_t step = 0;
    for (step = 0; step < num_refinement_steps; ++step) {

      // refine
      do_lateral_refinement (sparse_matrices, tree_weights, gap_factor);

      // log
      if (CTAGGING(8,ANNEALING REFINEMENT) && ((step % 5) == 0)) {
	CTAG(8,ANNEALING REFINEMENT) << "Completed " << step + 1 << " iterative refinement steps." << endl;
      }

      // did we improve?
      // if not, then stop iterative refinement
      const float score = expected_score (sparse_matrices, tree_weights, gap_factor);
      if (score > curr_score)
	curr_score = score;
      else {
	if (score + DOUBLE_TINY < curr_score) // warn if score has decreased (shouldn't be able to!)
	  CTAG(8,ANNEALING REFINEMENT) << "WARNING: Expected score dropped from " << curr_score << " to " << score << " during iterative refinement!" << endl;
	break;
      }

    }

    // log (if didn't just do so)
    if (CTAGGING(8,ANNEALING REFINEMENT) && ((step % 5) != 0)) {
      CTAG(8,ANNEALING REFINEMENT) << "Finished iterative refinement (" << step + 1 << " steps)." << endl;
    }

  }

  // close gui file
  if (output_for_gui) {
    gui_file.close();
    CTAG(9,ANNEALING) << "Created GUI output file '" << gui_filename << "'." << endl;
  }

}

void Alignment_DAG::show (std::ostream& o, const bool aligned_only /* = false */) const {

  for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {

    if ((*col)->is_dead())
      continue;
    if (aligned_only && ((*col)->get_seq_pos_map().size() < 2))
      continue;

    // show the column and its forward and back edges
    o << **col << endl;

    // incoming cols
    const Seq_pos_map& seq_pos_map = (*col)->get_seq_pos_map();
    for (Seq_pos_map::const_iterator seq_pos_iter = seq_pos_map.begin(); seq_pos_iter != seq_pos_map.end(); ++seq_pos_iter) {

      // have we reached the beginning of this sequence?
      // (remember that sequence coordinates are 0-based)
      if (seq_pos_iter->second == 0)
	continue;

      // get the previous position in this sequence
      const Seq_pos prev_seq_pos = Seq_pos (seq_pos_iter->first, seq_pos_iter->second - 1);

      // show this column
      const Column* w = get_seq_pos_col (prev_seq_pos);
      o << " <- " << *w << endl;

    }

    // check outgoing cols
    for (Seq_pos_map::const_iterator seq_pos_iter = seq_pos_map.begin(); seq_pos_iter != seq_pos_map.end(); ++seq_pos_iter) {

      // get the next position in this sequence
      const Seq_pos next_seq_pos = Seq_pos (seq_pos_iter->first, seq_pos_iter->second + 1);

      // have we reached the end of this sequence?
      if (next_seq_pos.second >= seq_db.get_seq (next_seq_pos.first).length())
	continue;

      // show this column
      const Column* w = get_seq_pos_col (next_seq_pos);
      o << " -> " << *w << endl;

    }

  }

}

void Alignment_DAG::write_gui_output (std::ostream& o, bool show_edges /* = false */) const {

  o << "; Initial DAG" << endl
    << "; Format is:" << endl
    << ";   column: (sequence, position) => initial_accuracy" << endl
    << "; sequence is 0-based and position is 0-based" << endl
    << endl;

  // first show column indices and the corresponding sequence positions, as well as initial accuracies
  for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {

    assert (!(*col)->is_dead());

    // show the column and its forward and back edges
    const Seq_pos_map& seq_pos_map = (*col)->get_seq_pos_map();
    assert (seq_pos_map.size() == 1);

    o << (*col)->get_orig_index() << ": (" << seq_pos_map.begin()->first << ", " << seq_pos_map.begin()->second << ")" << endl;
  }

  o << endl;

  // show edges if requested
  if (show_edges) {

    o << "; Initial edges" << endl;
    o << "; Format is:" << endl;
    o << ";   source_column -> dest_column" << endl;
    o << endl;

    // now show forward edges
    for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {

      const Seq_pos_map& seq_pos_map = (*col)->get_seq_pos_map();

      // get the next position in this sequence
      const Seq_pos next_seq_pos = Seq_pos (seq_pos_map.begin()->first, seq_pos_map.begin()->second + 1);

      // have we reached the end of this sequence?
      if (next_seq_pos.second >= seq_db.get_seq (next_seq_pos.first).length())
	continue;

      // else show this edge
      const Column* w = get_seq_pos_col (next_seq_pos);
      o << (*col)->get_orig_index() << " -> " << w->get_orig_index() << endl;

    }

    o << endl;
  }

}

bool Alignment_DAG::dfs_f (Column* l_bound, Column* u_bound, std::vector<Column*>& R_forward) {

  // init stack of nodes with the start node l_bound
  std::stack<Column*> nodes;
  const size_t offset = l_bound->get_index();
  std::vector<bool> on_stack ((u_bound->get_index() - offset + 1), false);

  nodes.push (l_bound);
  on_stack[l_bound->get_index() - offset] = true;

  while (!nodes.empty()) {

    // get next node on stack
    Column* node = nodes.top();
    nodes.pop();

    // if we've already seen node, then continue with the next
    if (node->is_marked())
      continue;

    // else mark as visited and store it
    node->mark();
    R_forward.push_back (node);

    // continue DFS: explore node's neighbors
    const Seq_pos_map& seq_pos_map = node->get_seq_pos_map();
    for (Seq_pos_map::const_iterator seq_pos_iter = seq_pos_map.begin(); seq_pos_iter != seq_pos_map.end(); ++seq_pos_iter) {

      // get the next position in this sequence
      const Seq_pos next_seq_pos = Seq_pos (seq_pos_iter->first, seq_pos_iter->second + 1);

      // have we reached the end of this sequence?
      if (next_seq_pos.second >= seq_db.get_seq (next_seq_pos.first).length())
	continue;

      // if not, then continue DFS on neighbors
      Column* w = get_seq_pos_col (next_seq_pos);

      // have we found a cycle?
      if (*w == *u_bound)
	return true;

      // if not, then continue DFS
      // it's critical that the conditional is ordered like this to prevent overflowing on_stack
      // if *w > *u_bound!
      else if (!(w->is_marked()) && (*w < *u_bound) && !on_stack[w->get_index() - offset]) {
	nodes.push (w);
	on_stack[w->get_index() - offset] = true;
      }

    }

  }

  // DFS terminated successfully, so no cycles found!
  return false;
}

void Alignment_DAG::dfs_b (Column* u_bound, Column* l_bound, std::vector<Column*>& R_backward) {

  // init stack of nodes with the start node u_bound
  std::stack<Column*> nodes;
  const size_t offset = l_bound->get_index();
  std::vector<bool> on_stack ((u_bound->get_index() - offset + 1), false);

  nodes.push (u_bound);
  on_stack[l_bound->get_index() - offset] = true;

  while (!nodes.empty()) {

    // get next node on stack
    Column* node = nodes.top();
    nodes.pop();

    // if we've already seen node, then continue with the next
    if (node->is_marked())
      continue;

    // else mark as visited and store it
    node->mark();
    R_backward.push_back (node);

    // continue DFS: explore node's neighbors
    const Seq_pos_map& seq_pos_map = node->get_seq_pos_map();
    for (Seq_pos_map::const_iterator seq_pos_iter = seq_pos_map.begin(); seq_pos_iter != seq_pos_map.end(); ++seq_pos_iter) {

      // have we reached the beginning of this sequence?
      // (remember that sequence coordinates are 0-based)
      if (seq_pos_iter->second == 0)
	continue;

      // get the previous position in this sequence (remember that we're searching backwards)
      Seq_pos prev_seq_pos = Seq_pos (seq_pos_iter->first, seq_pos_iter->second - 1);

      // if not, then continue DFS on neighbors
      Column* w = get_seq_pos_col (prev_seq_pos);

      // note that there's no need to check for cycles here,
      // since we've already done so in the dfs_f() call

      // continue DFS
      if (!(w->is_marked()) && (*l_bound < *w) && !on_stack[w->get_index() - offset]) {
	nodes.push (w);
	on_stack[w->get_index() - offset] = true;
      }
    }

  }

}

void Alignment_DAG::reorder (std::vector<Column*>& R_forward, std::vector<Column*>& R_backward) {

  // sort nodes in deltaB and deltaF
  std::sort (R_backward.begin(), R_backward.end(), smaller_index());
  std::sort (R_forward.begin(), R_forward.end(), smaller_index());

  // assemble list of indices from nodes in deltaB and deltaF
  // R_indices is R, the pool of available indices, in the Pearce and Kelly paper
  std::vector<size_t> R_indices;
  std::vector<Column*> L_columns;
  // deltaB
  for (std::vector<Column*>::const_iterator column = R_backward.begin(); column != R_backward.end(); ++column) {
    L_columns.push_back (*column);
    R_indices.push_back ((*column)->get_index());
  }
  // deltaF
  for (std::vector<Column*>::const_iterator column = R_forward.begin(); column  != R_forward.end(); ++column) {
    L_columns.push_back (*column);
    R_indices.push_back ((*column)->get_index());
  }
  // sort indices
  std::sort (R_indices.begin(), R_indices.end());

  // allocate indices to nodes in deltaB and deltaF
  for (size_t i = 0; i < R_indices.size(); ++i) {
    Column* col = L_columns[i];
    const size_t index = R_indices[i];
    assert (index < this->columns.size()); // remember 0-based indexing for Column::index
    col->set_index (index);
    this->columns[index] = col;
  }

}

void Alignment_DAG::unmark (std::vector<Column*>& cols) {
  for (std::vector<Column*>::iterator column = cols.begin(); column != cols.end(); ++column)
    (*column)->unmark();
}

void Alignment_DAG::merge (Edge* edge) {
  Column* source = edge->source;
  Column* dest = edge->dest;

  // add dest sequence positions to source
  const Seq_pos_map& dest_seq_pos_map = dest->get_seq_pos_map();
  for (Seq_pos_map::const_iterator seq_pos = dest_seq_pos_map.begin(); seq_pos != dest_seq_pos_map.end(); ++seq_pos) {
    source->add_seq_pos (*seq_pos);
    set_seq_pos_col (*seq_pos, source);
  }

  // merge dest into source
  dest->set_merged_into (source);

}

bool Alignment_DAG::induces_cycle (Edge* edge, Column* l_bound, Column* u_bound, std::vector<Column*>& R_forward) {

  // does this edge induce a cycle in the graph?
  if (dfs_f (l_bound, u_bound, R_forward)) {

    // sanity check: must be at least one node (l_bound)
    assert (R_forward.size());

    // unmark marked nodes
    unmark (R_forward);

    return true;
  }

  // sanity check: must be at least one node (l_bound)
  assert (R_forward.size());

  return false;

}

void Alignment_DAG::add_edge (Edge* edge, Column* l_bound, Column* u_bound, std::vector<Column*>& R_forward) {

  Column* source = edge->source;
  Column* dest = edge->dest;

  // NB: It's crucial that at this point the Edge be between
  // the current (after merges) source and dest columns;
  // this must be assured by the calling function.

  // we assume that induces_cycle has been called,
  // thereby ensuring that the edge doesn't induce a cycle
  // and populating R_forward; now we just need to perform
  // the backward depth-first search, storing nodes as the search proceeds
  std::vector<Column*> R_backward;
  dfs_b (u_bound, l_bound, R_backward);

  // sanity check: must be at least one node (l_bound and u_bound respectively)
  assert (R_forward.size());
  assert (R_backward.size());

  // unmark visited nodes
  unmark (R_forward);
  unmark (R_backward);

  // impose a new total order via the Pearce-Kelly Reorder() procedure
  if (R_forward.size() == 1) {    // catch edge cases
    source = u_bound;
    dest = l_bound;
  } else if (R_backward.size() == 1) {
    source = l_bound;
    dest = u_bound;
  } else {
    reorder (R_forward, R_backward);
  }

  // modify the DAG: merge the (current) source and dest columns
  edge->source = source;
  edge->dest = dest;
  merge (edge);

  // mark the dest column as dead
  // (note that this is /much/ faster than erasing it!)
  dest->set_dead();

  // decrement the number of columns in the alignment accordingly
  --num_columns;

}

// (sequence position to shift, (original column, new column))
typedef std::pair<Seq_pos, std::pair<size_t, size_t> > Shift;


void Alignment_DAG::do_lateral_refinement (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, float gap_factor) {

  // create priority queue of sequence positions to shift
  std::priority_queue<std::pair<float, Shift> > shifts;

  // create vector of (non-dead) columns
  std::vector<Column*> tmp_columns;
  for (std::vector<Column*>::iterator col = columns.begin(); col != columns.end(); ++col) {
    if (!(*col)->is_dead())
      tmp_columns.push_back (*col);
  }

  // for each column
  //   for each (seq, pos) in the column, consider shifting it
  //   left or right
  for (size_t i = 0; i < tmp_columns.size(); i++) {

    Column* curr_col = tmp_columns[i];

    // consider shifting each character
    const Seq_pos_map& seq_pos_map = curr_col->get_seq_pos_map();
    for (Seq_pos_map::const_iterator seq_pos = seq_pos_map.begin(); seq_pos != seq_pos_map.end(); ++seq_pos) {

      // column with the greatest change in expected accuracy from adding *seq_pos
      // initialize to the original column holding this Seq_pos
      size_t best_j = i;
      // best change in expected accuracy from adding *seq_pos to another column
      float best_score_change;

      // initialize best_score_change to change in score corresponding to current column
      // (corresponding to removing *seq_pos from curr_col, calculating expected change and then re-adding)
      best_score_change = curr_col->change_in_expected_score (*seq_pos, sparse_matrices, tree_weights, gap_factor, num_seqs, true); // skip_seq_pos = true

      // now search over other columns: search right
      for (size_t j = i + 1; j < tmp_columns.size(); ++j) {

	// get the column
	Column* col = tmp_columns[j];
	assert (!col->is_dead());
	while (col != col->get_merged_into())
	  col = col->get_merged_into();

	// if col contains a character from the sequence in *seq_pos, then break
	// (can't slide it over)
	const float change = col->change_in_expected_score (*seq_pos, sparse_matrices, tree_weights, gap_factor, num_seqs);
	if (change == INVALID_EDGE)
	  break;

	// if an improvement, then record it
	if (change > best_score_change) {
	  best_score_change = change;
	  best_j = j;
	}

      }

      // search left
      for (int j = i - 1; j >= 0; --j) {

	// get the column
	Column* col = tmp_columns[j];
	assert (!col->is_dead());
	while (col != col->get_merged_into())
	  col = col->get_merged_into();

	// if col contains a character from the sequence in *seq_pos, then break
	// (can't slide it over)
	const float change = col->change_in_expected_score (*seq_pos, sparse_matrices, tree_weights, gap_factor, num_seqs);
	if (change == INVALID_EDGE)
	  break;

	if (change > best_score_change) {
	  best_score_change = change;
	  best_j = j;
	}

      }

      // if we've found a better position, then record it
      if (best_j != i) {
	const Shift shift (*seq_pos, std::pair<size_t, size_t> (i, best_j));
	shifts.push (std::pair<float, Shift> (best_score_change, shift));
      }

    }
  }

  // now shift the (seq, pos) according to the corresponding
  // change in expected score
  // for each candidate shift, check that it's still a valid shift; 
  // also check whether there's a new, better shift
  std::pair<float, Shift> score_shift;
  while (!shifts.empty()) {

    // get the top candidate shift
    score_shift = shifts.top();
    shifts.pop();
    const Shift& shift = score_shift.second;
    const Seq_pos& seq_pos = shift.first;

    // get the original column holding seq_pos
    const size_t orig_j = (shift.second).first;

    // initialize best column to original column
    size_t best_j = orig_j;
    // best change in expected accuracy from adding *seq_pos to another column
    float best_score_change;

    // initialize best_score_change to change in score corresponding to original column
    // (corresponding to removing *seq_pos from curr_col, calculating expected change and then re-adding)
    // Note that this is different from Sudeep's original method,
    // which initialized the scores to 0 rather than the score associated with
    // the current column (and so, I believe, was biased towards shifting)
    best_score_change = tmp_columns[orig_j]->change_in_expected_score (seq_pos, sparse_matrices, tree_weights, gap_factor, num_seqs, true); // skip_seq_pos = true

    // search right
    for (size_t j = orig_j + 1; j < tmp_columns.size(); ++j) {

      // get the column
      Column* col = tmp_columns[j];
      assert (!col->is_dead());
      while (col != col->get_merged_into())
	col = col->get_merged_into();

      // if col contains a character from the sequence in *seq_pos, then break
      // (can't slide it over)
      const float change = col->change_in_expected_score (seq_pos, sparse_matrices, tree_weights, gap_factor, num_seqs);
      if (change == INVALID_EDGE)
	break;

      // if an improvement, then record it
      if (change > best_score_change) {
	best_score_change = change;
	best_j = j;
      }

    }

    // search left
    for (int j = orig_j - 1; j >= 0; --j) {

      // get the column
      Column* col = tmp_columns[j];
      assert (!col->is_dead());
      while (col != col->get_merged_into())
	col = col->get_merged_into();

      // if col contains a character from the sequence in *seq_pos, then break
      // (can't slide it over)
      const float change = col->change_in_expected_score (seq_pos, sparse_matrices, tree_weights, gap_factor, num_seqs);
      if (change == INVALID_EDGE)
	break;

      // if an improvement, then record it
      if (change > best_score_change) {
	best_score_change = change;
	best_j = j;
      }

    }

    // if wrong order on priority queue (if lower expected change in score than next shift on queue)
    if (best_score_change < (shifts.top()).first) {
      shifts.push (std::pair<float, Shift> (best_score_change, shift));
      continue;
    }

    // if we've found an improvement, then make the change and move seq_pos!
    if (best_j != orig_j) {
      if (CTAGGING(-1,REFINEMENT)) {
	CL << "Converting " << endl
	   << "  " << *tmp_columns[orig_j] << endl
	   << "  " << *tmp_columns[best_j] << endl;
      }
      tmp_columns[orig_j]->erase_seq_pos (seq_pos.first);
      tmp_columns[best_j]->add_seq_pos (seq_pos);
      set_seq_pos_col (seq_pos, tmp_columns[best_j]);
      if (CTAGGING(-1,REFINEMENT)) {
	CL << "to " << endl
	   << "  " << *tmp_columns[orig_j] << endl
	   << "  " << *tmp_columns[best_j] << endl;
      }
    }

  }

  // we're done, so create new columns object
  columns.clear();
  size_t idx = 0;
  for (size_t i = 0; i < tmp_columns.size(); ++i) {
    Column* col = tmp_columns[i];
    // if empty, then skip
    if (col->get_seq_pos_map().empty())
      continue;
    columns.push_back (col);     // store
    col->set_index (idx++);      // update index
  }
  num_columns = columns.size();

}

std::pair<float, std::string> Alignment_DAG::get_accuracy_annotation_normalized (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights) const {

  float overall_acc = 0.;
  float overall_denom = 0.;
  std::string annot;
  for (std::vector<Column*>::const_iterator column = columns.begin(); column != columns.end(); ++column) {

    // skip dead columns
    if ((*column)->is_dead())
      continue;

    // If normalizing, then convert to per-pair accuracy as:
    // Divide by 2 * (n choose 2 = 1/2*n*(n-1)) + ((num_seqs - n) * n),
    // where the first part is the (sum-of-pairs) number of matched characters
    // and the second part the number of characters aligned to gaps.
    // The factor of 10.0 here is because accuracy as calculated during annealing is <= 1,
    // and we want to map it to the interval [0-10].
    const std::pair<float, float> accuracy_pair = (*column)->get_accuracy_normalized (sparse_matrices, tree_weights, num_seqs);

    // cover case of no information for the column (can occur when a sequence is empty in pairwise inference)
    if (accuracy_pair.second == 0.) {
      annot += '0';
      continue;
    }

    // else perform calculation as usual
    const float acc_norm = accuracy_pair.first / accuracy_pair.second;
    assert ((acc_norm >= 0.) && (acc_norm <= 1.00001)); // prevent rounding error for the case of 1.0 (yes, it happens!)
    annot += static_cast<char> ((acc_norm > 0) ? "0123456789"[ std::min (static_cast<int> (std::floor (acc_norm * 10.0)), 9) ] : '0'); // normalize to 0-9 (0 if expected accuracy < 0)

    // increment for overall calculation
    overall_acc += accuracy_pair.first;
    overall_denom += accuracy_pair.second;

  }

  // cover case of no information for the alignment
  if (overall_denom == 0.)
    return std::make_pair (0., annot);

  // normalize by number of (non-dead) columns
  overall_acc /= overall_denom;
  assert ((overall_acc >= 0.) && (overall_acc <= 1.00001));

  return std::make_pair (overall_acc, annot);

}

float Alignment_DAG::expected_score (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const Tree_weights& tree_weights, const float gap_factor) const {

  float score = 0;
  for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {
    if ((*col)->is_dead())
      continue;
    score += (*col)->expected_score (sparse_matrices, tree_weights, gap_factor, num_seqs);
  }

  return score;
}


/**
 * \file sequence_pair_selector.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <queue>

#include "util/logfile.h"
#include "fsa/sequence_pair_selector.h"

using namespace fsa;

Sequence_pair_selector::Sequence_pair_selector (const Sequence_database& seq_db,
						const Alphabet& alphabet,
						const size_t k)
: __seq_db (seq_db),
  __alphabet (alphabet),
  __k (k),
  __similarity_matrix (Sequence_similarity_matrix (seq_db,
						   alphabet,
						   k)),
  __added (std::vector<std::vector<bool> > (seq_db.size(), std::vector<bool> (seq_db.size(), false))),
  __num_selected (0)
{

  // log
  CTAG(9,FSA) << "Computed k-mer similarities (k = " << k << ") between input sequences." << endl;

  // log similarities if requested
  if (CTAGGING(-1,KMER)) {
    CL << "Sequence similarity matrix (k = " << k << "):" << endl;
    __similarity_matrix.write (CL);
  }
  

}

void Sequence_pair_selector::choose_minimum_spanning_tree (const size_t num_mst,
							   const bool compute_maximum_spanning_tree /* = false */) {

  // iteratively build disjoint minimum spanning trees
  size_t mst_added = 0;
  while (mst_added++ < num_mst) {

    // assemble the edge set
    MST_Kruskal::Edge_set edges;
    edges.reserve (__seq_db.size() * __seq_db.size());
    for (size_t i = 0; i < __seq_db.size(); ++i)
      for (size_t j = i + 1; j < __seq_db.size(); ++j) {

	// don't store previously-added edges
	if (__added[i][j])
	  continue;

	// convert the k-mer similarity to a weight
	// (k-mer similarities necessarily fall in the interval [0, 1])
	double weight;

	// if we're computing a maximum spanning tree,
	// then store the weight as a similarity,
	// so that my minimum spanning tree code will
	// choose edges with the lowest similarity and thereby
	// build a maximum spanning tree
	if (compute_maximum_spanning_tree)
	  weight = __similarity_matrix.get_similarity (i, j);

	// else we're computing a minimum spanning tree,
	// so convert the k-mer similarity to a k-mer distance
	else
	  weight = 1.0 - __similarity_matrix.get_similarity (i, j);

	// store
	edges.push_back (MST_Kruskal::Edge (i, j, weight));

      }

    // compute the MST
    MST_Kruskal::Edge_set mst = MST_Kruskal::compute_mst (edges, __seq_db.size(),
							  false); // die_on_failure = false

    // if no MST exists, then break (we're done)
    if (!mst.size())
      break;

    // log MST if requested
    if (CTAGGING(-1,KMER)) {

      // compute the total k-mer similarity, summed over branches of the tree
      double total_similarity = 0;
      for (MST_Kruskal::Edge_set::const_iterator edge = mst.begin(); edge != mst.end(); ++edge)
	total_similarity += __similarity_matrix.get_similarity (edge->u, edge->v);

      // display
      if (compute_maximum_spanning_tree)
	CL << "Maximum spanning tree on sequences (total k-mer similarity = " << total_similarity << "):" << endl;
      else
	CL << "Minimum spanning tree on sequences (total k-mer similarity = " << total_similarity << "):" << endl;
      for (MST_Kruskal::Edge_set::const_iterator edge = mst.begin(); edge != mst.end(); ++edge)
	CL << __seq_db.get_seq (edge->u).name << " -- " << __seq_db.get_seq (edge->v).name << " => " << __similarity_matrix.get_similarity (edge->u, edge->v) << endl;

    }

    // store the MST
    for (MST_Kruskal::Edge_set::const_iterator edge = mst.begin(); edge != mst.end(); ++edge) {

      const size_t u = edge->u;
      const size_t v = edge->v;

      // sanity check
      assert (u < v);
      __added[u][v] = true;
      ++__num_selected;

    }

  }

}

void Sequence_pair_selector::choose_maximum_spanning_tree (const size_t num_mst) {

  choose_minimum_spanning_tree (num_mst,
				true);  // compute_maximum_spanning_tree = true

}

void Sequence_pair_selector::choose_minimum_spanning_palm_tree (const size_t num_mst) {

  // check sane
  if (num_mst > __seq_db.size())
    THROWEXPR ("Cannot compute more minimum spanning palm trees than there are sequences.");

  // keep track of which sequences we've used as centers of palm trees
  // (so that we don't repeat ourselves)
  std::vector<bool> chosen_centers (__seq_db.size(), false);

  // iteratively build minimum spanning palm trees
  size_t mst_added = 0;
  while (mst_added++ < num_mst) {

    // find the sequence at the center of the minimum spanning palm tree
    size_t best_center = 0;
    double best_similarity = -1;
    for (size_t i = 0; i < __seq_db.size(); ++i) {

      // avoid repeats
      if (chosen_centers[i])
	continue;

      // calculate the total similarity (summed over all branches)
      // associated with a palm tree whose center is sequence i
      double tree_similarity = 0;
      for (size_t j = 0; j < __seq_db.size(); ++j) {

	if (i == j)
	  continue;

	tree_similarity += __similarity_matrix.get_similarity (i, j);

      }

      // is this a new best center for the palm tree?
      if (best_similarity < tree_similarity) {
	best_similarity = tree_similarity;
	best_center = i;
      }

    }

    // confirm that we did indeed find a valid MST
    if (best_similarity < 0)
      THROWEXPR ("Unable to compute a valid minimum spanning palm tree.");

    // log if requested
    if (CTAGGING(-1,KMER))
      CL << "Minimum spanning palm tree on sequences is centered on sequence '" << __seq_db.get_seq (best_center).name << "' (total k-mer similarity = " << best_similarity << ")." << endl;

    // store the MST
    for (size_t j = 0; j < __seq_db.size(); ++j) {

      if (best_center == j)
	continue;

      // keep track of how many new sequence pairs we're selecting here
      // (avoid over-counting previously-selected ones)
      if (!__added[best_center][j] && !__added[j][best_center])
	++__num_selected;

      // maintain i < j for consistency
      if (best_center < j)
	__added[best_center][j] = true;
      else
	__added[j][best_center] = true;

    }

    // record the center so that we don't repeat ourselves
    chosen_centers[best_center] = true;

  }

}

void Sequence_pair_selector::choose_palm_tree (const std::string& centerseq) {

  // check sane
  if (!__seq_db.exists_seq (centerseq))
    THROWEXPR ("Requested center sequence of palm tree does not exist: '" << centerseq << "'.");

  const size_t center = __seq_db.get_seq_index (centerseq);

  // store the MST
  for (size_t j = 0; j < __seq_db.size(); ++j) {

    if (center == j)
      continue;

    // keep track of how many new sequence pairs we're selecting here
    // (avoid over-counting previously-selected ones)
    if (!__added[center][j] && !__added[j][center])
      ++__num_selected;

    // maintain i < j for consistency
    if (center < j)
      __added[center][j] = true;
    else
      __added[j][center] = true;

  }

}

void Sequence_pair_selector::choose_kmer_similarity (const size_t degree) {

  // keep track of how many sequence pairs we've stored
  size_t num_added = 0;

  // now add sequence pairs based on their k-mer similarity until
  // we've hit the requested number degree per sequence
  for (size_t i = 0; i < __seq_db.size(); ++i) {

    // assemble a heap of edges, ordered by descending k-mer similarity
    // (STL priority_queue order such that the first element is the greatest,
    // according to the designated strict weak ordering)
    std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t> >, Util::Duple_less<double, size_t> > weighted_edges;
    for (size_t j = i + 1; j < __seq_db.size(); ++j) {
      if (!__added[i][j])
	weighted_edges.push (std::make_pair (__similarity_matrix.get_similarity (i, j), j));
    }

    // what's the current degree of sequence i?
    size_t degree_i = 0;
    for (size_t j = 0; j < __seq_db.size(); ++j) {
      if (__added[i][j] || __added[j][i])
	++degree_i;
    }

    // now add edges one by one until we reach the requested degree
    // (or run out of edges)
    while (!weighted_edges.empty()) {

      // have we reached the requested degree for sequence i?
      if (degree_i >= degree)
	break;

      // get the lowest-weight edge      
      const std::pair<double, size_t> weighted_edge = weighted_edges.top();
      weighted_edges.pop();
      const size_t j = weighted_edge.second;

      // store it
      assert (!__added[i][j]);
      __added[i][j] = true;
      ++degree_i;
      ++__num_selected;
      ++num_added;
      
    }

  }

  // log
  if (CTAGGING(-1,KMER) && degree > 1)
    CL << "Stored " << num_added << " sequence pairs to reach the requested degree = " << degree << "." << endl;

}

void Sequence_pair_selector::choose_random (size_t num) {

  // make sure that it's possible to store the desired number of sequence pairs without duplicates
  if ((__num_selected + num) > (__seq_db.size() * (__seq_db.size() - 1) / 2))
    num = (__seq_db.size() * (__seq_db.size() - 1) / 2) - __num_selected;

  // then randomly add more until we hit num_alignment_pairs
  size_t num_added = 0;
  while (num_added < num) {
      
    // randomly choose a sequence pair
    size_t i = static_cast<size_t> (std::floor (Util::rand (__seq_db.size() - 1)));
    size_t j = static_cast<size_t> (std::floor (Util::rand (__seq_db.size() - 1)));

    // check sane
    assert (i < __seq_db.size());
    assert (j < __seq_db.size());

    // skip identities
    if (i == j)
      continue;

    // require that i < j for consistency
    if (i > j)
      std::swap (i, j);
    assert (i < j);

    // skip duplicates
    if (__added[i][j])
      continue;

    // store
    __added[i][j] = true;
    ++num_added;
    ++__num_selected;

  }
  assert (num_added == num);

  // log
  CTAG(8,FSA) << "Stored " << num_added << " randomly-chosen sequence pairs." << endl;

}

void Sequence_pair_selector::choose_all() {

  for (size_t i = 0; i < __added.size(); ++i) {
    for (size_t j = i + 1; j < __added.size(); ++j)
      __added[i][j] = true;
  }

}

Sequence_pairs Sequence_pair_selector::get_chosen_sequence_pairs() const {

  Sequence_pairs sequence_pairs;
  for (size_t i = 0; i < __added.size(); ++i) {
    for (size_t j = i + 1; j < __added.size(); ++j)
      if (__added[i][j])
	sequence_pairs.push_back (Sequence_pair (i, j));
  }

  return sequence_pairs;

}

double Sequence_pair_selector::erdos_renyi_p_cutoff (const size_t nodes) {

  const double epsilon = 0.1;
  return ((1 + epsilon) * std::log (static_cast<double> (nodes)) / static_cast<double> (nodes));

}


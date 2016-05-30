
/**
 * \file mst.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include "math/mst.h"

using namespace fsa;

Disjoint_set::Disjoint_set (const size_t n) {
  
  __elements.reserve (n);
  for (size_t e = 0; e < n; ++e)
    __elements.push_back (Element (e, 0)); // parent = e, rank = 0

}

Disjoint_set::Index Disjoint_set::find_dj (const Disjoint_set::Index e) {

  Element& element = get_element (e);
  if (element.parent == e)
    return e;
  else {
    element.parent = find_dj (element.parent);
    return element.parent;
  }    

}

void Disjoint_set::union_dj (const Disjoint_set::Index e1, const Disjoint_set::Index e2) {

  Disjoint_set::Index root1index = find_dj (e1);
  Element& root1 = get_element (root1index);
  Disjoint_set::Index root2index = find_dj (e2);
  Element& root2 = get_element (root2index);

  if (root1.rank > root2.rank)
    root2.parent = root1index;
  else if (root1.rank < root2.rank)
    root1.parent = root2index;
  // unless e1 and e2 are already in the same set, merge them
  else if (root1 != root2) {
    root2.parent = root1index;
    root1.rank = root1.rank + 1;
  }

}

MST_Kruskal::Edge_set MST_Kruskal::compute_mst (MST_Kruskal::Edge_set& edges, const size_t N,
						const bool die_on_failure /* = true */) {

  MST_Kruskal::Edge_set mst;
  mst.reserve (N);

  // initialize forest:
  // each vertex is a tree
  Disjoint_set dj (N);

  // sort edges, ordered by increasing weight
  std::sort (edges.begin(), edges.end());

  // loop over edges, lowest weight first
  for (MST_Kruskal::Edge_set::const_iterator edge = edges.begin(); edge != edges.end(); ++edge) {

    // get vertices connected by this edge
    const size_t u = edge->u;
    const size_t v = edge->v;

    // check that the vertices are labeled as 0, ..., (N - 1)
    assert ((u < N) && (v < N));

    // if connecting u and v does not form a cycle
    // (if the trees containing them are distinct),
    // then do so
    if (dj.find_dj (u) != dj.find_dj (v)) {
      mst.push_back (*edge);
      dj.union_dj (u, v);
    }

    // have we assembled a MST?
    if (mst.size() == N - 1)
      break;

  }

  // confirm that we've created a spanning tree
  // (i.e., that the graph was complete to begin with)
  if (mst.size() != N - 1) {
    if (die_on_failure) {
      cerr << "ERROR: Unable to compute a minimal spanning tree.  Was the graph complete to begin with?" << endl;
      exit (1);
    } else
      return MST_Kruskal::Edge_set();
  }

  return mst;

}


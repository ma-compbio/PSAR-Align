
/**
 * \file mst.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef MATH_MST_INCLUDED
#define MATH_MST_INCLUDED

#include "util/misc.h"
#include "math/mathematics.h"

namespace fsa {

  /**
   * \brief Represent a disjoint set with UNION-FIND capabilities.
   *
   * Representation is a disjoint set forest, where each set is represented
   * by a tree and each node in the tree (element in the set) holds a reference
   * to the root of the tree, the representative element of the set.
   * This is essentially an implementation of the wikipedia description:
   *  http://en.wikipedia.org/wiki/Disjoint-set_data_structure
   */
  struct Disjoint_set {

  public:

    typedef size_t Index;

    /**
     * \brief Constructor.
     *
     * Create a disjoint set with the specified number of elements;
     * each element is in its own set.
     * Implicitly implements the make_set operation.
     */
    Disjoint_set (const size_t n);

    /**
     * \brief Find the set containing the element e.
     *
     * Can also be used to determine whether two elements are in the same set.
     * Implements the path compression heuristic.
     * \param e index of element
     * \return index of the representative of the set
     */
    Index find_dj (const Index e);

    /**
     * \brief Merge the sets containing two elements into a single set.
     *
     * Implements the union by rank heuristic.
     * Awkward name is because "union" is reserved by C++.
     * \param e1 index of element in set 1 to be merged
     * \param e2 index of element in set 2 to be merged
     */
    void union_dj (const Index e1, const Index e2);

  private:

    /**
     * \brief Represent a single element in the disjoint set.
     */
    struct Element {

      Element (const Index parent, const size_t rank)
      : parent (parent), rank (rank)
      { }

      Index parent;
      size_t rank;

      /**
       * \brief Equality iff parents are identical.
       * \see find_dj
       */
      bool operator== (const Element& e) const {
	return parent == e.parent;
      }

      bool operator!= (const Element& e) const {
	return !(*this == e);
      }

    };

    /**
     * \brief Get an element by its index.
     */
    Element& get_element (const Index e) {
      assert (e < __elements.size());
      return __elements[e];
    }

    std::vector<Element> __elements;  ///< elements in the disjoint set

  };

  /**
   * \brief Perform's Kruskal's algorithm for finding the minimum spanning tree of a graph.
   */
  struct MST_Kruskal {

  public:

    /**
     * \brief Represent a weighted, undirected edge.
     */
    struct Edge {

      /**
       * \brief Constructor.
       */
      Edge (const unsigned short u, const unsigned short v,
	    const float weight)
      : u (u), v (v), weight (weight)
      { }

      unsigned short u;
      unsigned short v;
      float weight;
      
      /**
       * \brief Comparison operator.
       */
      bool operator< (const Edge& e) const {
	return (weight < e.weight);
      }

      /**
       * \brief Comparison operator.
       */
      bool operator> (const Edge& e) const {
	return (weight > e.weight);
      }

      /**
       * \brief Output operator.
       */
      friend std::ostream& operator<< (std::ostream& o, const Edge& edge) {
	o << edge.u << " -- " << edge.v << " => " << edge.weight;
	return o;
      }

    };

    typedef std::vector<Edge> Edge_set;  ///< a set of edges

    /**
     * \brief Find the minimum spanning tree on the given set of edges.
     *
     * The vertices in edges MUST be numbered as 0,...,(N - 1)!
     * \param edges list of edges of graph
     * \param N number of nodes in the graph
     * \param die_on_failure die if the graph is not complete (a MST does not exist)
     * \return the edges composing the MST (empty if no MST found)
     */
    static Edge_set compute_mst (Edge_set& edges, const size_t N,
				 const bool die_on_failure = true);

  };

}

#endif /* MATH_MST_INCLUDED */

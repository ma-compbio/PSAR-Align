
/**
 * \file sequence_pair_selector.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef FSA_SEQUENCE_PAIR_SELECTOR_INCLUDED
#define FSA_SEQUENCE_PAIR_SELECTOR_INCLUDED

#include "math/mst.h"
#include "seq/sequence.h"
#include "seq/similarity_matrix.h"

namespace fsa {

  typedef std::pair<size_t, size_t> Sequence_pair;     ///< a sequence pair (i, j)
  typedef std::vector<Sequence_pair> Sequence_pairs;   ///< a vector of sequence pairs

  /**
   * \brief Select sequence pairs for pairwise comparisons.
   */
  struct Sequence_pair_selector {

  public:

    /**
     * \brief Constructor.
     * \param seq_db sequence data
     * \param alphabet alphabet which sequences are defined over
     * \param k word length to use when computing k-mer similarities
     */
    Sequence_pair_selector (const Sequence_database& seq_db,
			    const Alphabet& alphabet,
			    const size_t k);

    /**
     * \brief Choose sequence pairs based on their inclusion in minimum spanning trees.
     *
     * Iteratively constructs disjoint minimum spanning trees.
     * \param num_mst number of disjoint MSTs to construct
     * \param compute_maximum_spanning_tree compute the maximum, rather than minimum, spanning tree
     */
    void choose_minimum_spanning_tree (const size_t num_mst,
				       const bool compute_maximum_spanning_tree = false);

    /**
     * \brief Choose sequence pairs based on their inclusion in maximum spanning trees.
     *
     * Iteratively constructs disjoint maximum spanning trees.
     * \param num_mst number of disjoint MSTs to construct
     * \see choose_minimum_spanning_tree
     */
    void choose_maximum_spanning_tree (const size_t num_mst);

    /**
     * \brief Choose sequence pairs based on their inclusion in the minimum palm spanning tree.
     *
     * Iteratively construct minimum palm spanning trees with distinct centers.
     * \param num_palm_mst number of palm MSTs to construct
     */
    void choose_minimum_spanning_palm_tree (const size_t num_palm_mst);

    /**
     * \brief Choose sequence pairs which form a palm tree with the specified sequence at the center.
     *
     * \param seq sequence which is the center of the palm tree
     * 
     */
    void choose_palm_tree (const std::string& centerseq);

    /**
     * \brief Choose sequence pairs based on their k-mer similarity.
     *
     * \param degree every sequence must be included in degree sequence pairs
     */
    void choose_kmer_similarity (const size_t degree);

    /**
     * \brief Choose sequence pairs randomly.
     *
     * \param num number of sequence pairs to add randomly
     */
    void choose_random (size_t num);

    /**
     * \brief Choose all sequence pairs.
     */
    void choose_all();

    /**
     * \brief Count the number of selected sequence pairs.
     */
    size_t num_selected() const { return __num_selected; }

    /**
     * \brief Get a list of sequence pairs which we have selected.
     */
    Sequence_pairs get_chosen_sequence_pairs() const;

    /**
     * \brief Erdos-Renyi threshold probability.
     *
     * Probability above which an Erdos-Renyi random graph will almost surely be connected.
     * \param nodes number of nodes in graph
     */
    static double erdos_renyi_p_cutoff (const size_t nodes);

  private:

    const Sequence_database& __seq_db;  ///< sequence data
    const Alphabet& __alphabet;         ///< alphabet which sequences are defined over
    const size_t __k;                   ///< k-mer word length

    Sequence_similarity_matrix __similarity_matrix;  ///< k-mer similarity matrix
    std::vector<std::vector<bool> >__added;          ///< __added[i][j] == true iff sequence pair (i, j) or (j, i) is chosen
    size_t __num_selected;                           ///< number of sequence pairs selected

  };

}

#endif /* FSA_SEQUENCE_PAIR_SELECTOR_INCLUDED */

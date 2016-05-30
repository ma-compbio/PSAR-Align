
/**
 * \file tree_weights.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

// to do:
// - 

#ifndef ANNEALING_TREE_WEIGHTS_INCLUDED
#define ANNEALING_TREE_WEIGHTS_INCLUDED

#include <iostream>
#include <map>

#include "util/hash_fcn.h"
#include "seq/sequence.h"

namespace fsa {

  /**
   * \brief Represent a set of (symmetric) weights for sequence pairs.
   * 
   * Generally these weights will be derived using a tree relating
   * the sequences.
   * If used as an "empty class," without calling the method from_file
   * to populate the data structures, then the weights will all be set to 1.
   * This allows the method get_weight to be used transparently both when
   * weights are given and when they are not.
   */
  struct Tree_weights {

  public:

    /**
     * \brief Constructor.
     *
     * Returns all weights as 1 until populated with from_file.
     */
    Tree_weights()
    : __populated (false)
    { }

    /**
     * \brief Read weights from a file.
     *
     * Unless particular weights are present in the file,
     * sets them to 0.  Assumes that weights are symmetric,
     * i.e., that weight (i, j) == weight (j, i).
     */
    void from_file (const Sequence_database& seq_db, const std::string& filename);

    /**
     * \brief Normalize the weights such that they sum to (N choose 2).
     *
     * (N choose 2) is the number of iterations in a sum-of-pairs
     * operation over N sequences.
     */
    void normalize();

    /**
     * \brief Get the weight for a particular sequence pair.
     *
     * This function can be called on an empty instance of this class,
     * i.e., one for which the data structures have not been populated.
     * If they have not been, then all weights returned will be 1.
     * \param i index of first sequence
     * \param j index of second sequence
     * \return weight of pair, or 0 if undefined (or always 1 if used as an empty class)
     */
    inline float operator() (const size_t i, const size_t j) const;

    static const double min_sequence_weight;           ///< minimum weight for a sequence (summed over pairs)


  private:

    bool __populated;                            ///< have we read input data? (all weights = 1 unless true)
    std::vector<std::vector<float> > __weights;  ///< weights for each sequence pair
    
  };

  inline float Tree_weights::operator() (const size_t i, const size_t j) const {
    // return 1 if we haven't populated the weights
    if (!__populated)
      return 1.0;
    // else look up the appropriate weighta
    assert (i < __weights.size());
    assert (j < __weights[i].size());
    return __weights[i][j];
  }

}

#endif /* ANNEALING_TREE_WEIGHTS_INCLUDED */

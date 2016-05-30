
/**
 * \file constraints.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef FSA_CONSTRAINTS_INCLUDED
#define FSA_CONSTRAINTS_INCLUDED

#include <iostream>
#include <fstream>

#include <map>

#include "util/misc.h"
#include "seq/sequence.h"

namespace fsa {

  /**
   * \brief Represents a single constraint.
   *
   * Assumes a 0-based coordinate system and a fully-closed interval.
   * Constraints are like Anchors, but they are assumed to already be consistent.
   */
  struct Constraint {

    Interval xcoords;  ///< interval of the constraint in sequence X
    Interval ycoords;  ///< interval of the constraint in sequence Y

    /**
     * \brief Constructor.
     */
    Constraint (const Interval& xcoords, const Interval& ycoords)
    : xcoords (xcoords), ycoords (ycoords) { }
  
    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<< (std::ostream& o, const Constraint& constraint)  {
      o << "[" << constraint.xcoords.start << ", " << constraint.xcoords.end << "] ~ ["
	<< constraint.ycoords.start << ", " << constraint.ycoords.end << "]";
      return o;
    }

    /**
     * \brief Write constraint in Mercator format.
     *
     * Note that Mercator constraints are 0-based half-open [start, end),
     * whereas we represent them as [start, end].
     */
    inline void write_mercator (const std::string& xname, const std::string& yname, std::ostream& o) const {
      o << xname << ' ' << xcoords.start << ' ' << xcoords.end + 1 << ' '
	<< yname << ' ' << ycoords.start << ' ' << ycoords.end + 1 << endl;
    }

  };

  /**
   * \brief Hold a vector of constraints for two sequences.
   */
  struct Constraints {

    /**
     * \brief Constructor.
     */
    Constraints() { }

    /**
     * Constructor.
     *
     * \param xname name of sequence X
     * \param yname name of sequence X
     */
    Constraints (const std::string& xname, const std::string& yname)
      : __xname (xname), __yname (yname) { }

    /**
     * \brief Write constraint file in Mercator format.
     */
    void write_mercator (std::ostream& o) const;

    /**
     * \brief Number of constraints.
     */
    size_t size() const { return __constraints.size(); }

    /**
     * \brief Store a constraint.
     */
    void store (const Constraint& c) {
      __constraints.push_back (c);
    }

    /**
     * \brief Get iterator to start of __constraints.
     */
    std::vector<Constraint>::iterator begin() {
      return __constraints.begin();
    }

    /**
     * \brief Get iterator to start of __constraints.
     */
    std::vector<Constraint>::const_iterator begin() const {
      return __constraints.begin();
    }

    /**
     * \brief Get iterator to end of __constraints.
     */
    std::vector<Constraint>::iterator end() {
      return __constraints.end();
    }

    /**
     * \brief Get iterator to end of __constraints.
     */
    std::vector<Constraint>::const_iterator end() const {
      return __constraints.end();
    }


  private:
  
    const std::string __xname;               ///< name of sequence X
    const std::string __yname;               ///< name of sequence Y

    std::vector<Constraint> __constraints;   ///< Contraint objects

  };

  /**
   * \brief Represent a set of contraints for many sequences.
   */
  struct Constraints_set {

    /**
     * \brief Constructor.
     *
     * \param seq_db input sequence data
     */
    Constraints_set (const Sequence_database& seq_db);

    /**
     * \brief Read constraints from a Mercator file.
     *
     * Note that Mercator constraints are 0-based half-open [start, end).
     */
    void read_mercator (const std::string& filename);

    /**
     * \brief Write constraints to a Mercator-format file.
     *
     * Note that Mercator constraints are 0-based half-open [start, end).
     */
    void write_mercator (std::ostream& o) const;

    /**
     * \brief Get constraints for a particular sequence std::pair.
     *
     * \param i index for first sequence
     * \param j index for second sequence
     */
    inline Constraints get_constraints (unsigned i, unsigned j) const {
      if (constraints_list.find (std::make_pair (i, j)) == constraints_list.end())
	return Constraints();
      return (*constraints_list.find (std::make_pair (i, j))).second;
    }

    /**
     * \brief Number of sequence std::pairs which we have constraints for.
     */
    inline size_t size() const {
      return constraints_list.size();
    }

  private:

    /**
     * \brief Hold all input sequence data.
     */
    const Sequence_database& seq_db;

    /**
     * \brief constraints for sequence std::pairs
     */
    std::map<std::pair<size_t, size_t>, Constraints> constraints_list;

  };

}

#endif /* FSA_CONSTRAINTS_INCLUDED */

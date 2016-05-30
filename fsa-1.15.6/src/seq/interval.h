
/**
 * \file interval.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef SEQ_INTERVAL_INCLUDED
#define SEQ_INTERVAL_INCLUDED

#include "config.h"
#include "util/misc.h"
#include "seq/gff.h"

namespace fsa {

  /**
   * \brief Represent a 0-based, fully-closed interval.
   */
  struct Interval {

  public:

    Interval()
    { }

    /**
     * \brief Constructor.
     */
    Interval (const unsigned start, const unsigned end)
    : start (start), end (end)
    { }

    unsigned start;           ///< start coordinate of interval
    unsigned end;             ///< end coordinate of interval

    /**
     * \brief Length of interval.
     */
    size_t length() const;

  };

  /**
   * \brief Represent a single genomic interval.
   *
   * As with all FSA code unless otherwise noted,
   * this object is intended to use 0-based, fully-closed coordinates.
   */
  struct Genomic_interval : public Interval {

  public:

    Genomic_interval() { }

    /**
     * \brief Constructor.
     */
    Genomic_interval (const std::string& genome, const std::string& chromosome,
		      const unsigned start, const unsigned end, const char strand = GFF::unknown_strand)
      : Interval (start, end),
      genome (genome), chromosome (chromosome),
      strand (strand)
    { }

    std::string genome;       ///< genome
    std::string chromosome;   ///< chromosome name
    char strand;              ///< strand

    /**
     * \brief Write in a Mercator-like format.
     */
    std::string to_string() const {
      std::string str;
      return str;
    }

    /**
     * \brief Convert Genomic_intervals to GFFs.
     *
     * Assumes that the coordinates in intervals are 0-based and converts them to 1-based coordinates.
     */
    static GFF_database convert_to_gff_db (const std::vector<Genomic_interval>& intervals);

    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<< (std::ostream& o, const Genomic_interval& interval) {
      o << ((interval.genome == "") ? "." : interval.genome) << '\t' << interval.chromosome + '\t'
	<< Util::to_string (interval.start) << '\t' << Util::to_string (interval.end) << '\t'
	<< interval.strand;
      return o;
    }

    /**
     * \brief Compare two Genomic_interval objects based on their starting coordinates.
     */
    bool operator< (const Genomic_interval& r) const {
      if (genome == r.genome) {
	if (chromosome == r.chromosome)
	  return (start < r.start);
	return chromosome < r.chromosome;
      }
      return genome < r.genome;
    }

    /**
     * \brief Function object for binary comparison of Genomic_interval objects.
     *
     * \see Genomic_interval::operator<
     */
    struct Genomic_interval_less : std::binary_function<Genomic_interval, Genomic_interval, bool> {
    public:
      bool operator() (const Genomic_interval& l, const Genomic_interval& r) const {
	return l < r;
      }
    };

  };

  inline size_t Interval::length() const {
    // catch case of an empty interval
    if (start > end)
      return 0;
    return end - start + 1;
  }

}

#endif /* SEQ_INTERVAL_INCLUDED */


/**
 * \file dotplot.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Ian Holmes, Lars Barquist and Robert Bradley.
 */

#ifndef DOTPLOT_INCLUDED
#define DOTPLOT_INCLUDED

#define PROB_DISPLAY_CUTOFF 0.01

#include "util/misc.h"
#include "util/array2d.h"
#include "seq/sequence.h"

namespace fsa {

  /**
   * \brief Base class for probability matrices.
   *
   * Both alignment and fold (for e.g. Pair SCFGs) dotplots inherit from this class.
   */
  struct Dotplot : array2d<double> {

    std::string xseq; ///< x axis label
    std::string yseq; ///< y axis label

    Dotplot (size_t x, size_t y)
      : array2d<double> (x, y, 0.) { }

    /**
     * \brief Constructor.
     */
    Dotplot (const Sequence& x, const Sequence& y);

    /**
     * \brief Output method.
     *
     * For e.g. alignment dotplots, output is formatted as:
     * <cruft>
     * .   x1   x2   x3
     * y1  .2   .2   .2
     * y2  .2   .2   .2
     */
    void write_dotplot (const std::string& filename, const double cutoff = PROB_DISPLAY_CUTOFF) const;

    /**
     * \brief Output method for fold dotplots.
     */
    void write_dotplot (const std::string& prefix, const std::string& seqname, const double cutoff = PROB_DISPLAY_CUTOFF) const;

    /**
     * \brief Output method for alignment dotplots.
     */
    void write_dotplot (const std::string& prefix, const std::string& xseqname, const std::string& yseqname, const double cutoff = PROB_DISPLAY_CUTOFF) const;

  };

}

#endif /*DOTPLOT_INCLUDED*/	

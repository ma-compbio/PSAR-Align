
/**
 * \file SparseMatrix.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Chuong Do.
 * Robert Bradley wrote the Post_probs struct, corresponding
 * SparseMatrix constructor and write_gui_output.
 * Jaeyoung Do wrote a SparseMatrix constructor, 
 * add_probability and finishing_adding_probability
 * to construct a SparseMatrix with probabilities transfered
 * from Database or MW workers  on-the-fly.
 */

#ifndef SPARSEMATRIX_INCLUDED
#define SPARSEMATRIX_INCLUDED

// Note that the SparseMatrix class uses 1-based indexing for positions in sequences.
// This is IN CONTRAST to the rest of FSA, which uses 0-based indexing.

#define DOUBLE_TINY 0.001
#define DOUBLE_VERY_TINY 0.000001
#define PRECISION_DEFAULT 5

#include <iostream>
#include <iomanip>

#include "util/logfile.h"
#include "manager/db_misc.h"

namespace fsa {

  /**
   * \brief Sparse matrix entry type
   *
   * (column, value)
   */
  typedef std::pair<unsigned, float> Matrix_entry;

  /**
   * \brief Representation of a single posterior probability.
   *
   * Assumes 0-based coordinates.
   */
  struct Post_prob {

    unsigned x;     ///< coordinate in first sequence
    unsigned y;     ///< coordinate in second sequence
    float prob;     ///< posterior probability entry

    /**
     * \brief Constructor.
     */
    Post_prob (unsigned x, unsigned y, float prob)
    : x (x), y (y), prob (prob)
    { }

    /**
     * \brief Lexical order on the coordinates x and y.
     */
    bool operator< (const Post_prob& r) const {
      if (x == r.x)
	return y < r.y;
      else
	return x < r.x;
    }

    /**
     * \brief Output operator.
     *
     * Output in 0-based coordinates (which Post_probs uses).
     */
    friend std::ostream& operator<< (std::ostream& o, const Post_prob& post_prob)  {
      o << "[" << post_prob.x << ", " << post_prob.x << "] ~ [" << post_prob.y << ", " << post_prob.y << "]" << " => "
	<< std::setprecision (PRECISION_DEFAULT) << post_prob.prob;
      return o;
    }

  };

  /**
   * \brief Sparse representation of a set of posterior probabilities.
   */
  typedef std::vector<Post_prob> Post_probs;


  /**
   * \brief Sparse matrix.
   *
   * Uses 1-based coordinates internally.
   */
  class SparseMatrix {

    size_t xidx;                                                ///< numerical index of the first sequence
    size_t yidx;                                                ///< numerical index of the second sequence
    size_t seq1Length, seq2Length;                              ///< dimensions of matrix
    std::vector<size_t> rowSize;                                ///< rowSize[i] = # of cells in row i
    std::vector<Matrix_entry> data;                             ///< data values
    std::vector<std::vector<Matrix_entry>::iterator> rowPtrs;   ///< pointers to the beginning of each row

    std::vector<float> gapPosteriors;                           ///< gap posteriors (sum of pairs) for first and second sequences

    /**
     * \brief Default constructor.
     */
    SparseMatrix()
      { }

  public:

#ifdef HAVE_POSTGRES
    /**
     * \brief Constructor.
     *
     * Builds a sparse matrix with posterior probabilities from Database.
     */
    SparseMatrix (const size_t xidx, const size_t yidx,
		  const size_t xlen, const size_t ylen,
		  const size_t num_cells)
      : xidx (xidx), yidx (yidx),
      seq1Length (xlen), seq2Length (ylen) {
	
      // ensure sane
      assert (xlen >= 0);
      assert (ylen >= 0);

      // allocate memory
      data.resize (num_cells);
      rowSize.assign (xlen + 1, 0); rowSize[0] = 0;       // catch case of no entry!
      rowPtrs.resize (xlen + 1); rowPtrs[0] = data.end();
      gapPosteriors.resize (xlen + ylen + 2, 1.0);   // initialize the gap posteriors

    }

    void add_probability (const int pos1, const int pos2, const float prob, const int offset) {

      std::vector<Matrix_entry>::iterator dataPtr = data.begin();
      dataPtr+=offset;

      if (pos1 == 0) 
	gapPosteriors[seq1Length + pos2 + 1] = prob; 
      else if (pos2 == 0) 
	gapPosteriors[pos1] = prob;
      else {
	if (rowSize[pos1] == 0)
	  rowPtrs[pos1] = dataPtr;

	rowSize[pos1]++;

	dataPtr->first  = pos2;
	dataPtr->second = prob;
      }
    }
#endif

#ifdef HAVE_CONDOR
	
    /**
     * \brief Constructor.
     *
     * Builds a sparse matrix with posterior probabilities from workers.
     */
    SparseMatrix (const size_t xidx, const size_t yidx,
		  const size_t xlen, const size_t ylen,
		  int pre_len, int len, Sparse_Matrix_Buffer *sm_buffer)
      : xidx (xidx), yidx (yidx),
      seq1Length (xlen), seq2Length (ylen) {

      // ensure sane
      assert (xlen >= 0);
      assert (ylen >= 0);

      // allocate memory
      data.resize (len - xlen - ylen);
      rowSize.assign (xlen + 1, 0); rowSize[0] = 0;       // catch case of no entry!
      rowPtrs.resize (xlen + 1); rowPtrs[0] = data.end();
      gapPosteriors.resize (xlen + ylen + 2, 1.0);   // initialize the gap posteriors

      std::vector<Matrix_entry>::iterator dataPtr = data.begin();

      int row=0;	

      for (int cnt = pre_len; cnt < (pre_len + len); cnt++) {
	if (sm_buffer[cnt].pos1 == 0)  
	  gapPosteriors[seq1Length + sm_buffer[cnt].pos2 + 1] = sm_buffer[cnt].prob; 
	else if (sm_buffer[cnt].pos2 == 0)  
	  gapPosteriors[sm_buffer[cnt].pos1] = sm_buffer[cnt].prob;
	else {

	  if (row != sm_buffer[cnt].pos1) {
	    rowSize[row] = dataPtr - rowPtrs[row];
	    rowPtrs[sm_buffer[cnt].pos1] = dataPtr;
	    row = sm_buffer[cnt].pos1;
	  }
	  dataPtr->first = sm_buffer[cnt].pos2;
	  dataPtr->second = sm_buffer[cnt].prob;
	  dataPtr++;
	}
      } //end for cnt 
      rowSize[row] = dataPtr - rowPtrs[row]; //for the last tuple
    }
#endif

    /**
     * \brief Constructor.
     *
     *  Builds a sparse matrix from a Post_probs representation
     * of a posterior matrix.  The coordinates stored in the Post_prob entries
     * are assumed to be 0-based; note that this is in contrast to the 1-based coordinates
     * used internally by the SparseMatrix class.
     * The entries Post_prob must be in lexical order, ie
     * sorted by first coordinate (x), then second (y).
     * DOES NOT enforce a cutoff on posterior probabilities;
     * because the input is sparse, assumes that all entries should be kept.
     */
    SparseMatrix (const size_t xidx, const size_t yidx,
		  const size_t xlen, const size_t ylen,
		  const Post_probs& post_probs)
      : xidx (xidx), yidx (yidx),
      seq1Length (xlen), seq2Length (ylen) {

      // ensure sane
      assert (xlen >= 0);
      assert (ylen >= 0);

      // calculate memory required; count the number of cells in the posterior matrix
      const size_t num_cells = post_probs.size();

      // allocate memory
      data.resize (num_cells);
      rowSize.assign (xlen + 1, 0); rowSize[0] = 0;       // catch case of no entry!
      rowPtrs.resize (xlen + 1); rowPtrs[0] = data.end();
      gapPosteriors.resize (xlen + ylen + 2, 1.0);   // initialize the gap posteriors

      if (CTAGGING(-1, SPARSEMATRIX)) {
	CTAG(-1, SPARSEMATRIX) << "Initializing SparseMatrix:" << endl;
      }

      // read the data
      int xprev = -1; // the previous x-coordinate
      int xcurr = -1; // the current x-coordinate
      int yprev = -1; // the previous y-coordinate
      std::vector<Matrix_entry>::iterator dataPtr = data.begin(); // iterator through SparseMatrix::data
      // Note: When initializing the SparseMatrix like this, it's crucial that we iterate through the x-coordinate 
      // in increasing order => Post_probs must be sorted properly.  If we don't do it in order
      // then the pointers will get messed up.  Which is Very Bad.
      for (Post_probs::const_iterator iter = post_probs.begin(); iter != post_probs.end(); ++iter) {

	if (CTAGGING(-1, SPARSEMATRIX)) {
	  CL << " " << *iter << endl;
	}

	// pull out coordinates and convert to 1-based coordinates
	const int x = iter->x + 1;
	const int y = iter->y + 1;
	float prob = iter->prob;

	// check sane (now we're in 1-based coordinates)
	assert ((x >= 1) && (static_cast<size_t> (x) <= xlen));
	assert ((y >= 1) && (static_cast<size_t> (y) <= ylen));

	// don't allow prob to overflow (results in negative gap posteriors => negative weights => very bad)
	prob = (prob > 1.) ? 1. : prob;

	// if we just changed x coordinate
	if (x != xprev) {

	  // store the row size for the previous x coordinate
	  if (xprev != -1) // catch case of first row
	    rowSize[xprev] = dataPtr - rowPtrs[xprev];

	  // then continue on to the new x coordinate
	  xcurr = x;
	  rowPtrs[x] = dataPtr;
	  yprev = -1; // reset

	}

	// ensure that entries are ordered properly (we require lexical ordering)
	assert (xprev <= xcurr);
	// ensure no duplicate entries
	assert (yprev < y);

	// store this entry
	dataPtr->first = y;
	dataPtr->second = prob;
	++dataPtr;

	// decrement gap posteriors accordingly
	gapPosteriors[x] -= prob;
	gapPosteriors[xlen + y + 1] -= prob;

	// sanity checks: gap posteriors can't be negative beyond roundoff error
	// (ensure that we're actually storing a probability distribution)
	if (gapPosteriors[x] < -DOUBLE_TINY) {
	  (*this).Print (CL, true);
	  THROWEXPR ("ERROR: Negative gap posterior: gapPosteriors[x = " << x << "] = " << gapPosteriors[x]);
	}
	if (gapPosteriors[xlen + y + 1] < -DOUBLE_TINY) {
	  (*this).Print (CL, true);
	  THROWEXPR ("ERROR: Negative gap posterior: gapPosteriors[y = " << y << "] = " << gapPosteriors[xlen + y + 1]);
	}
	// if acceptable roundoff error, then round up to zero
	if (gapPosteriors[x] < 0.)
	  gapPosteriors[x] = 0.;
	if (gapPosteriors[xlen + y + 1] < 0.)
	  gapPosteriors[xlen + y + 1] = 0.;

	// prevent overflow
	if (gapPosteriors[x] < 1e-4)
	  gapPosteriors[x] = 1e-4;
	if (gapPosteriors[xlen + y + 1] < 1e-4)
	  gapPosteriors[xlen + y + 1] = 1e-4;

	// increment our indices
	xprev = xcurr;
	yprev = y;
      }

      // catch edge case of the last row
      if (xcurr != -1)
	rowSize[xcurr] = dataPtr - rowPtrs[xcurr];

    }

    /**
     * \brief Get row pointer.
     *
     * \return pointer to a particular row in the sparse matrix
     */
    std::vector<Matrix_entry>::iterator GetRowPtr (const unsigned row) const {
      assert (row >= 1 && row <= seq1Length);
      return rowPtrs[row];
    }

    /**
     * \brief Get match probability in matrix.
     *
     * \param row (0-based) row position
     * \param col (0-based) column position
     * \return value at a particular row, column
     */
    float get_match_prob (const unsigned row, const unsigned col) const {
      assert (row >= 1 && row <= seq1Length);
      assert (col >= 1 && col <= seq2Length);
      for (size_t i = 0; i < rowSize[row]; i++){
	if (rowPtrs[row][i].first == col) 
	  return rowPtrs[row][i].second;
      }
      return 0;
    }

    /**
     * \brief Get row size.
     *
     * \return number of entries in a particular row.
     */
    size_t GetRowSize (const unsigned row) const {
      assert (row >= 1 && row <= seq1Length);
      return rowSize[row];
    }

    /**
     * \brief Get length of X (sequence 1).
     *
     * \return first dimension of the matrix
     */
    size_t GetSeq1Length() const {
      return seq1Length;
    }

    /**
     * \brief Get length of Y (sequence 2).
     *
     * \return second dimension of the matrix
     */
    size_t GetSeq2Length() const {
      return seq2Length;
    }

    /// Get number of entries in matrix.
    /*
     * \return number of entries
     */
    size_t size() const {
      return data.size();
    }

    /**
     * \brief Write contents of object.
     *
     * Position indices are 0-based in output.
     * The SparseMatrix class uses 1-based sequence indexing internally;
     * this is converted to 0-based indexing when formatting output.
     */
    void Print (std::ostream &outfile, bool show_gaps = false) const {
      outfile << "Sparse match posteriors (0-based coords):" << endl;
      for (unsigned i = 1; i <= seq1Length; i++) {
	// skip zero entries
	if (rowSize[i] == 0)
	  continue;
	outfile << "  " << i - 1 << ":"; // convert to 0-based indexing
	for (unsigned j = 0; j < rowSize[i]; j++) {
	  outfile << " (" << rowPtrs[i][j].first - 1 << "," << rowPtrs[i][j].second << ")"; // convert to 0-based indexing
	}
	outfile << endl;
      }
      if (show_gaps) {
	outfile << "Gap posteriors 0 (0-based coords): ";
	for (unsigned i = 1; i <= seq1Length; i++){
	  outfile << " (" << i - 1 << "," << gapPosteriors[i] << ")";                  // convert to 0-based indexing
	}
	outfile << endl << "Gap posteriors 1 (0-based coords): ";
	for (unsigned i = 1; i <= seq2Length; i++){
	  outfile << " (" << i - 1 << "," << gapPosteriors[i + seq1Length + 1] << ")"; // convert to 0-based indexing
	}
	outfile << endl;
      }
    }

    /**
     * \brief Write output formatted for GUI.
     *
     * Sequence and position indices are both 0-based in output.
     * The SparseMatrix class uses 1-based sequence indexing internally;
     * this is converted to 0-based indexing when formatting output.
     */
    void write_gui_output (std::ostream& o) const {

      o << "; Sparse posterior probability matrix for sequences " << xidx << " and " << yidx << endl
	<< "; Format is:" << endl
	<< ";   (" << xidx << ", position_1) ~ (" << yidx << ", position_2) => prob" << endl
	<< "; which means that (" << xidx << ", position_1) is aligned to (" << yidx << ", position_2) with probability prob." << endl
	<< ";   (" << xidx << ", position_1) ~ (" << yidx << ", -1) => prob" << endl
	<< "; means that (" << xidx << ", position_1) is aligned to a gap in " << yidx << " with probability prob." << endl
	<< "; sequence is 0-based and position is 0-based" << endl
	<< endl;

      o << "; match posteriors" << endl;
      for (unsigned i = 1; i <= seq1Length; i++) {
	// skip zero entries
	if (rowSize[i] == 0)
	  continue;
	for (unsigned j = 0; j < rowSize[i]; j++) {
	  o << "(" << xidx << ", " << i - 1 << ") ~ (" << yidx << ", " << rowPtrs[i][j].first - 1 << ") => " << rowPtrs[i][j].second << endl; // convert to 0-based indexing
	}
      }
      o << endl;

      o << "; gap posteriors" << endl;
      for (unsigned i = 1; i <= seq1Length; i++)
	o << "(" << xidx << ", " << i - 1 << ") ~ (" << yidx << ", " << "-1" << ") => " << gapPosteriors[i] << endl; // convert to 0-based indexing
      o << endl;
      for (unsigned j = 1; j <= seq2Length; j++)
	o << "(" << xidx << ", " << "-1" << ") ~ (" << yidx << ", " << j - 1 << ") => " << gapPosteriors[seq1Length + j + 1] << endl; // convert to 0-based indexing

      o << endl;
    }

    /**
     * \brief Compute transpose of matrix.
     *
     * \return SparseMatrix which is the tranpose of the current matrix
     */
    SparseMatrix* ComputeTranspose() const {

      // create a new sparse matrix
      SparseMatrix* transpose = new SparseMatrix();
      const size_t numCells = data.size();

      transpose->xidx = yidx;
      transpose->yidx = xidx;
      transpose->seq1Length = seq2Length;
      transpose->seq2Length = seq1Length;

      // allocate memory
      transpose->data.resize (numCells);
      transpose->rowSize.resize (seq2Length + 1); transpose->rowSize[0] = 0;
      transpose->rowPtrs.resize (seq2Length + 1); transpose->rowPtrs[0] = transpose->data.end();
      transpose->gapPosteriors.resize(seq1Length + seq2Length + 2);

      // compute row sizes
      for (unsigned i = 1; i <= seq2Length; i++)
	transpose->rowSize[i] = 0;
      for (unsigned i = 0; i < numCells; i++)
	transpose->rowSize[data[i].first]++;

      // compute row ptrs
      for (unsigned i = 1; i <= seq2Length; i++)
	transpose->rowPtrs[i] = (i == 1) ? transpose->data.begin() : transpose->rowPtrs[i-1] + transpose->rowSize[i-1];

      // now fill in data
      std::vector<std::vector<Matrix_entry>::iterator> currPtrs = transpose->rowPtrs;

      for (unsigned i = 1; i <= seq1Length; i++){
	std::vector<Matrix_entry>::iterator row = rowPtrs[i];
	for (unsigned j = 0; j < rowSize[i]; j++){
	  currPtrs[row[j].first]->first = i;
	  currPtrs[row[j].first]->second = row[j].second;
	  currPtrs[row[j].first]++;
	}
      }

      for (unsigned i = 0; i <= seq1Length; i++)
	transpose->gapPosteriors[i + seq2Length + 1] = gapPosteriors[i];
      for (unsigned i = 0; i <= seq2Length; i++)
	transpose->gapPosteriors[i] = gapPosteriors[i + seq1Length + 1];

      return transpose;
    }

    /**
     * \brief Get gap probability in matrix.
     *
     * \param which 0 for sequence 0, 1 for sequence 1
     * \param pos (1-based) position in sequence
     * \return gap probability for pos
     */
    float get_gap_prob (const unsigned which, const unsigned pos) const {

      // check no overflow of sequence lengths
      assert ((which == 0 && pos <= seq1Length) || (which == 1 && pos <= seq2Length));

      return gapPosteriors[pos + which * (seq1Length + 1)];
    }

  };

}

#endif /* SPARSEMATRIX_INCLUDED */


/**
 * \file alignment.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 * Parsing of Stockholm-format alignments is based on Ian Holmes's
 * Stockholm class.
 */

#ifndef SEQ_ALIGNMENT_INCLUDED
#define SEQ_ALIGNMENT_INCLUDED

#include <set>
#include <algorithm>

#include "seq/sequence.h"

namespace fsa {

  /**
   * \brief Represent a single row of an alignment.
   *
   * All coordinates are 0-based.
   * Intervals are always fully-closed, [start, end].
   */
  struct Alignment_row {

  public:

    /**
     * \brief Represention of row of alignment.
     *
     * True indicates that a character is present at an alignment position;
     * false otherwise.
     */
    typedef std::vector<bool> Row_path;

    /**
     * \brief Constructor.
     */
    Alignment_row (Row_path* path);
    Alignment_row (const Alignment_row& parent);
    Alignment_row& operator= (const Alignment_row& parent);

    /**
     * \brief Destructor.
     */
    ~Alignment_row();

    /**
     * \brief Build mapping between sequence and alignment coordinates.
     */
    void build_coordinate_indices();

    /**
     * \brief Length of alignment row.
     */
    inline size_t length() const { return __row_path->size(); }

    /**
     * \brief Length of ungapped sequence.
     */
    inline size_t seqlength() const { return __seqlength; }

    /**
     * \brief Get gapped sequence corresponding to this row.
     * \param ungapped ungapped sequence
     */
    Sequence get_gapped_seq (const Sequence& ungapped, const char gap_char) const;

    /**
     * \brief Get (possibly-gapped) character for a column of this row.
     * \param col column of Alignment_row
     */
    inline char get_char (const size_t col,
			  const Sequence& ungapped, const char gap_char) const;

    /**
     * \brief Extract subalignment.
     * 
     * \param start start coordinate of alignment row
     * \param end end coordinate of alignment row
     */
    Alignment_row* subalignment (const unsigned start, const unsigned end) const;

    /**
     * \brief Extract Row_path subalignment.
     * 
     * \param start start coordinate of alignment row
     * \param end end coordinate of alignment row
     */
    Row_path* subalignment_row_path (const unsigned start, const unsigned end) const;

    /**
     * \brief Map a sequence position to the corresponding alignment coordinate.
     */
    inline unsigned map_seq_to_align (const unsigned pos) const;

    /**
     * \brief Map subalignment coordinates to the corresponding sequence interval.
     *
     * \param start start coordinate of alignment row
     * \param end end coordinate of alignment row
     * \return sequence interval, or (1, 0) if the interval is empty
     */
    inline Interval map_align_to_seq (const unsigned start, const unsigned end) const;

    /**
     * \brief Reverse alignment.
     *
     * Calls build_coordinate_indices.
     */
    void reverse();

    /**
     * \brief Is an alignment coordinate (column) gapped in this sequence?
     */
    inline bool is_gapped (const unsigned col) const {
      assert (col < __row_path->size());
      return !(*__row_path)[col];
    }

    /**
     * \brief Output operator (write Row_path).
     */
    friend std::ostream& operator<< (std::ostream& o, const Alignment_row& row) {
      o << Util::join (*row.__row_path, "") << endl;
      return o;
    }

  private:

    std::vector<size_t> __seq_to_align_coords_map;  ///< map sequence to alignment coordinates
    std::vector<size_t> __align_to_seq_coords_map;  ///< map alignment to sequence coordinates (next ungapped position)
    Row_path* __row_path;                           ///< alignment path
    size_t __seqlength;                             ///< length of ungapped sequence

  };

  /**
   * \brief Represent a matrix-form global alignment.
   *
   * All coordinates are 0-based.
   * An Alignment holds a reference to the Sequence_database
   * which contains the actual sequence information.
   * The reference is non-const because an Alignment
   * can modify the sequence information (e.g., when it reads 
   * in an alignment).
   */
  struct Alignment {

  public:

    static const char gap_char;   ///< gap character

    /**
     * \brief Constructor.
     * \param seq_db Sequence_database to hold actual sequence information
     */
    Alignment (Sequence_database& seq_db);
    Alignment (const Alignment& parent);
    Alignment& operator= (const Alignment& parent);

    /**
     * \brief Detect whether an alignment seems to be in multi-FASTA format.
     * \see Sequence::detect_fasta
     */
    static bool detect_mfa (const std::string& filename);

    /**
     * \brief Initialize from multi-FASTA alignment.
     *
     * Clears all alignment and sequence data.
     * Populates seq_db as appropriate,
     * being very careful to maintain ordering of sequences
     * for consistency when indexing seq_db and __rows.
     * \param filename alignment filename
     * \param strict Require that alignment be flush.
     */
    void read_mfa (const std::string& filename, const bool strict = true,
		   const bool verbose = true);

    /**
     * \brief Write in multi-FASTA format.
     * \param strict Require that alignment be flush.
     */
    void write_mfa (std::ostream& o, const bool strict = true) const;

    /**
     * \brief Assert all alignment rows of equal length.
     */
    void assert_flush() const;

    /**
     * \brief Map a sequence position to the corresponding alignment coordinate.
     * \param seq sequence name
     * \see Alignment_row::map_seq_to_align
     */
    inline size_t map_seq_to_align (const std::string& seq, const unsigned pos) const {
      return get_row (seq).map_seq_to_align (pos);
    }

    /**
     * \brief Map subalignment coordinates to the corresponding sequence interval.
     *
     * Maps to the sequence interval which is strictly contained by
     * the subalignment.
     * \param seq sequence name
     * \param start start coordinate of alignment row
     * \param end end coordinate of alignment row
     * \return sequence interval, or (1, 0) if the interval is empty
     * \see Alignment_row::map_align_to_seq
     */
    inline Interval map_align_to_seq (const std::string& seq, const unsigned start, const unsigned end) const {
      return get_row (seq).map_align_to_seq (start, end);
    }

    /**
     * \brief Map an interval in one sequence to the corresponding interval in another.
     *
     * Maps to the sequence interval in seq_to which is strictly contained by
     * the subalignment implied by the sequence interval in seq_from.
     * \return sequence interval, or (1, 0) if the interval is empty
     * \see map_seq_to_align
     * \see map_align_to_seq
     */
    inline Interval map_seq_to_seq (const std::string& seq_from, const unsigned start_from, const unsigned end_from,
				    const std::string& seq_to) const;

    /**
     * \brief Calculate average per-column percentage identity of alignment.
     *
     * Reports the average per-column percent id:
     *  sum_columns (# identical characters in column / total # non-gap characters in column) / num_columns
     * Gaps are ignored in both the numerator and denominator.
     */
    double percent_id() const;

    /**
     * \brief Calculate average per-column percentage identity of subalignment [start, end].
     *
     * Reports the average per-column percent id:
     *  sum_columns (# identical characters in column / total # non-gap characters in column) / num_columns
     * Gaps are ignored in both the numerator and denominator.
     * \param start start column of subalignment
     * \param end end column of subalignment
     */
    double percent_id (const unsigned start, const unsigned end) const;

    /**
     * \brief Calculate fraction of gaps in the alignment.
     */
    double gap_fraction() const;

    /**
     * \brief Calculate fraction of gaps in the subalignment [start, end].
     *
     * \param start start column of subalignment
     * \param end end column of subalignment
     */
    double gap_fraction (const unsigned start, const unsigned end) const;

    /**
     * \brief Get row of alignment.
     * \param name sequence name
     */
    const Alignment_row& get_row (const std::string& name) const;

    /**
     * \brief Get gapped row of alignment formatted for display.
     * \param name sequence name
     */
    const Sequence get_gapped_row (const std::string& name) const;

    /**
     * \brief Get gapped row of alignment formatted for display.
     * \param r index in __rows or seq_db
     */
    const Sequence get_gapped_row (const size_t r) const;

    /**
     * \brief Remove all gaps from the passed string.
     */
    static void remove_gaps (std::string& str);

    /**
     * \brief Get name of row of alignment.
     * \param r index in __rows or seq_db
     */
    const std::string& get_row_name (const size_t r) const {
      assert (r < __rows.size());
      return seq_db->get_seq (r).name;
    }

    /**
     * \brief Get index of row in alignment.
     * \param name sequence name
     */
    size_t get_row_index (const std::string& name) const {
      assert (exists_row (name));
      return __row_index.find (name)->second;
    }

    /**
     * \brief Get names of rows in the alignment.
     *
     * Row names are ordered as they are indexed in the alignment.
     */
    std::vector<std::string> get_row_names() const;

    /**
     * \brief Add a row to the alignment.
     * 
     * Helper for reading alignment files from disk.
     * Generally preferred over the other add_row for efficiency.
     * Add sequence information to Sequence_database as well
     * if nothing is present for the sequence name.
     * Beware: This function may not behave as you expect.
     * \param sequence sequence for alignment row
     * \param row_path Alignment_row::row_path for row
     */
    void add_row (Sequence* sequence,
		  Alignment_row::Row_path* row_path);

    /**
     * \brief Add a row to the alignment.
     * 
     * Helper for reading alignment files from disk.
     * Add sequence information to Sequence_database as well
     * if nothing is present for the sequence name.
     * Beware: This function may not behave as you expect.
     * \param sequence sequence for alignment row
     * \param alignment_row Alignment_row for row
     */
    void add_row (Sequence* sequence,
		  Alignment_row* alignment_row);

    /**
     * \brief Set a row in the alignment.
     * 
     * This function should be used when constructing a new alignment
     * of sequence data which has already been stored in 
     * the Sequence_database&; if you want to read in an alignment
     * and sequence data simultaneously (e.g., when reading an alignment
     * file from disk), then you should use add_row.
     * 
     * Assumes (and requires) that there is a sequence named 'name'
     * in the stored Sequence_database; dies if not.
     * If there is a corresponding row already stored in __rows,
     * then overwrite it; otherwise store 
     * Note that if improperly used, this function can cause
     * the sequence indices for the Alignment object and the 
     * stored Sequence_database& object to go out of sync!
     * Be very careful.
     * \param name name of seqence for alignment row
     * \param row_path Alignment_row::row_path for row
     */
    void set_row (const std::string& name,
		  Alignment_row::Row_path* row_path);

    /**
     * \brief Do we have a row for a particular sequence name?
     * \param name sequence name
     */
    bool exists_row (const std::string& name) const {
      return (__row_index.find (name) != __row_index.end());
    }

    /**
     * \brief Clear all alignment and sequence data.
     *
     * \see Sequence_database::clear
     */
    virtual void clear();

    /**
     * \brief Number of sequences in alignment.
     */
    inline size_t rows() const;

    /**
     * \brief Number of columns in alignment.
     */
    inline size_t columns() const;

    /**
     * \brief Is a character a gap?
     * 
     * Gaps can be: gap_char . _
     */
    static inline bool is_gap_char (char c) {
      return (c == gap_char || c == '.' || c == '_');
    }

    /**
     * \brief Is a particular (row, column) gapped?
     * \param row row index for sequence
     * \param col column of alignment
     */
    inline bool is_gapped (const size_t row, const size_t col) const;

  protected:

    /**
     * \brief Get row of alignment.
     * \param r index in __rows or seq_db
     */
    const Alignment_row& get_row (const size_t r) const {
      assert (r < __rows.size());
      return *__rows[r];
    }

    /**
     * \brief Get row of alignment.
     * \param r index in __rows or seq_db
     */
    Alignment_row& get_row (const size_t r) {
      assert (r < __rows.size());
      return *__rows[r];
    }

    /**
     * \brief Get the character at a particular (row, column).
     * \param row row index for sequence
     * \param col column of alignment
     * \see Alignment_row::get_char
     */
    char get_char (const size_t row, const size_t col) const;

    /**
     * \brief Destructor.
     * 
     * Note that the memory for seq_db is NOT freed,
     * since we assume that it was allocated elsewhere.
     */
    virtual ~Alignment();

    Sequence_database* seq_db;                   ///< sequences of the alignment
    std::vector<Alignment_row*> __rows;          ///< individual rows of the alignment
    std::map<std::string, size_t> __row_index;   ///< map from row names to position in __rows

  };

  /**
   * \brief Represent a Stockholm-format alignment.
   */
  struct Stockholm : public Alignment {

    /**
     * \brief Constructor.
     */
    Stockholm (Sequence_database& seq_db);

    /**
     * \brief Initialize from Stockholm alignment.
     *
     * Clears all alignment and sequence data.
     * Populates seq_db as appropriate,
     * being very careful to maintain ordering of sequences
     * for consistency when indexing seq_db and __rows.
     * \param filename alignment filename.
     * \param strict Require that alignment be flush.
     * \see Alignment::clear
     */
    void read_stockholm (const std::string& filename, const bool strict = true,
			 const bool verbose = true);

    /**
     * \brief Initialize from Stockholm or multi-FASTA alignment.
     *
     * Clears all alignment and sequence data.
     * \param filename alignment filename.
     * \param strict Require that alignment be flush.
     * \see Alignment::clear
     */
    void read_stockholm_or_mfa (const std::string& filename, const bool strict = true,
				const bool verbose = true);

    /**
     * \brief Write in Stockholm format.
     * \param strict Require that alignment be flush.
     */
    void write_stockholm (std::ostream& o, const bool strict = true) const;

    /**
     * \brief Extract subalignment.
     * 
     * \param seq_db_subalign Sequence_database to store sequence information for subalignment
     * \param start start coordinate of alignment
     * \param end end coordinate of alignment
     */
    Stockholm* subalignment (Sequence_database& seq_db_subalign,
			     const unsigned start, const unsigned end) const;

    /**
     * \brief Get all #=GF lines with a given key.
     */
    std::string get_gf_annot (const std::string& key) const;

    /**
     * \brief Get #=GC line with a given key.
     */
    std::string get_gc_annot (const std::string& key) const;

    /**
     * \brief Get #=GS line for a given sequence with a given key.
     */
    std::string get_gs_annot (const std::string& seq, const std::string& key) const;

    /**
     * \brief Get #=GR line for a given sequence with a given key.
     */
    std::string get_gr_annot (const std::string& seq, const std::string& key) const;

    /**
     * \brief Add a #=GF line.
     */
    void add_gf_annot (const std::string& key, const std::string& value);

    /**
     * \brief Set a #=GC line.
     */
    void set_gc_annot (const std::string& key, const std::string& value);

    /**
     * \brief Set a #=GS line.
     */
    void set_gs_annot (const std::string& seq, const std::string& key, const std::string& value);

    /**
     * \brief Set a #=GR line.
     */
    void set_gr_annot (const std::string& seq, const std::string& key, const std::string& value);

    /**
     * \brief Clear all alignment and sequence data as well as annotations.
     *
     * \see Sequence_database::clear
     */
    void clear();

    /**
     * \brief Clear all annotation lines.
     */
    void clear_annot();

    /**
     * \brief Assert all alignment rows and per-column annotations of equal length.
     */
    void assert_all_flush() const;

    /**
     * \brief Reverse-complement alignment and associated sequence data under the passed alphabet.
     *
     * Handles annotations, etc. properly.
     * (This is why it's a method of Stockholm rather than the parent Alignment class.)
     * \param alphabet Alphabet to reverse-complement under
     */
    void revcomp (const Alphabet& alphabet);

    /**
     * \brief Map amino acid alignment to corresponding codon alignment.
     *
     * Assumes that the sequences in seq_db_codon are ordered as in this->seq_db.
     * Handles annotations, etc. properly (overhangs are annotated as '.').
     * \param seq_db_codon sequence data in nucleotide space
     */
    Stockholm get_codon_from_aa_alignment (Sequence_database& seq_db_codon) const;

    /**
     * \brief Annotate alignment with alignment statistics.
     * 
     * Assumes that the (implicit) alphabet is NOT case-sensitive.
     * \see Alignment::percent_id
     * \see Alignment::gap_fraction
     */
    void annotate_with_statistics();


    static const std::string gff_annotation;                         ///< GFF annotation key
    static const std::string percent_id_annotation;                  ///< percent id annotation key
    static const std::string gap_fraction_annotation;                ///< gap fraction annotation key



  private:

    static const std::string format_identifier;
    static const std::string version_identifier;
    static const std::string alignment_header;
    static const std::string alignment_separator;

    static const std::string file_annotation;
    static const std::string column_annotation;
    static const std::string sequence_annotation;
    static const std::string sequence_column_annotation;

    static const char annotation_wildcard_char;                      ///< Annotation unknown character

    typedef std::map<std::string, std::string> Annotation;           ///< keyed annotations
    typedef std::map<std::string, Annotation> Row_annotation;        ///< Annotation keyed by sequence

    std::vector<std::pair<std::string, std::string> > __gf_annot;    ///< sorted list of #=GF annotations
    Annotation __gc_annot;                                           ///< per-column annotation
    Row_annotation __gs_annot;                                       ///< per-sequence annotation
    Row_annotation __gr_annot;                                       ///< per-sequence, per-column annotation
    std::map<std::string, std::set<size_t> > __gf_index;             ///< map from keys to #=GF lines

  };


  /****************************************
   * Function definitions.
   ****************************************/

  inline unsigned Alignment_row::map_seq_to_align (const unsigned pos) const {
    assert (pos < __seq_to_align_coords_map.size());
    return __seq_to_align_coords_map[pos];
  }

  inline Interval Alignment_row::map_align_to_seq (const unsigned start, const unsigned end) const {

    if (start > end)
      cerr << "start = " << start << ", end = " << end << endl;
    assert (start <= end);
    assert (end < __align_to_seq_coords_map.size());

    const unsigned s = __align_to_seq_coords_map[start]; // next ungapped position
    unsigned e = __align_to_seq_coords_map[end];

    // check for e == 0
    if (e == 0) {

      // if e == 0 and sequence is gapped at column end,
      // then [s, e] must map to a run of gaps at the start of the alignment
      if (is_gapped (end))
	return Interval (1, 0);

      // if sequence isn't gapped at column end,
      // then we must have a situation like this:
      // 000000123
      // -----AGGC
      //   s  e
      // where column end corresponds to the first ungapped position in the sequence
      // very very tricky...
      else {
	assert (s == e);
	assert (e == 0);
	return Interval (0, 0);
      }

    }

    // from here on we can safely assume that e > 0

    // check for end being gapped;
    // if it is, then set e to the previous ungapped position
    // (this also catches case of e running off the end of the sequence)
    if (is_gapped (end))
      --e;

    // catch cases giving empty sequence:
    // - no sequence data in alignment
    // - case of e < s (interval must be all gaps)
    if ((__seqlength == 0) || (e < s))
      return Interval (1, 0);

    return Interval (s, e);
  }

  inline void Alignment::remove_gaps (std::string& str) {

    std::string::iterator last_pos = std::remove_if (str.begin(), str.end(),
						     Alignment::is_gap_char);
    str.erase (last_pos, str.end());

  }

  inline size_t Alignment::rows() const {
    assert (seq_db->size() == __rows.size());
    assert (__rows.size() == __row_index.size());
    return __rows.size();
  }

  inline size_t Alignment::columns() const {
    if (__rows.size())
      return __rows[0]->length();
    return 0;
  }

  inline Interval Alignment::map_seq_to_seq (const std::string& seq_from, const unsigned start_from, const unsigned end_from,
					     const std::string& seq_to) const {

    assert (start_from <= end_from);

    const unsigned start_align = map_seq_to_align (seq_from, start_from);
    const unsigned end_align = map_seq_to_align (seq_from, end_from);

    return map_align_to_seq (seq_to, start_align, end_align);

  }

  inline bool Alignment::is_gapped (const size_t row, const size_t col) const {

    return get_row (row).is_gapped (col);

  }

  inline char Alignment_row::get_char (const size_t col,
				       const Sequence& ungapped, const char gap_char) const {

    assert (col < __align_to_seq_coords_map.size());

    if (is_gapped (col))
      return gap_char;

    return ungapped.seq[__align_to_seq_coords_map[col]];

  }

  inline char Alignment::get_char (const size_t row, const size_t col) const {

    return get_row (row).get_char (col,
				   seq_db->get_seq (row), Alignment::gap_char);
  }

}

#endif /* SEQ_ALIGNMENT_INCLUDED */

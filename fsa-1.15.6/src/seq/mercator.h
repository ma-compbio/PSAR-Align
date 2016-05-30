
/**
 * \file mercator.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef SEQ_MERCATOR_INCLUDED
#define SEQ_MERCATOR_INCLUDED

#include "config.h"
#include "util/regexp.h"
#include "seq/alignment.h"
#include "seq/gff.h"
#include "seq/interval.h"

namespace fsa {

  /**
   * \brief Represent a single genomic interval in a mapping.
   *
   * Uses 0-based, fully-closed coordinates.
   * Note that this is NOT the same as Mercator itself,
   * which uses 0-based, half-open [start, end) coordinates.
   * \see Mercator_map
   */
  struct Mercator_interval : public Genomic_interval {

  public:

    Mercator_interval() { }

    /**
     * \brief Constructor.
     */
    Mercator_interval (unsigned bin,
		       const std::string& genome, const std::string& chromosome,
		       const unsigned start, const unsigned end, char strand)
      : Genomic_interval (genome, chromosome,
			  start, end, strand),
      bin (bin)
    { }

    unsigned bin;   ///< Mercator bin number

    /**
     * \brief Write in Mercator format.
     */
    std::string to_string() const {
      std::string str;
      str = chromosome + '\t'
	+ Util::to_string (start) + '\t' + Util::to_string (end) + '\t'
	+ strand;
      return str;
    }

    /**
     * \brief Output operator.
     */
    friend std::ostream& operator<< (std::ostream& o, const Mercator_interval& interval) {
      o << interval.bin << '\t'
<< interval.genome << '\t'
	<< interval.to_string() << endl;
      return o;
    }

    // explicitly import less_than
    using Genomic_interval::operator<;

  };

  /**
   * \brief Represent a Mercator mapping between genomes.
   *
   * Uses 0-based, fully-closed coordinates.
   * Mercator uses 0-based, half-open [start, end)
   * coordinates; the conversion is done internally by this class.
   * If present, strips off the leading 'chr' frequently prepended 
   * to chromosome names (this helps to prevent frustrations due
   * to, e.g., having 'chr' prepended to the chromosomes in the map
   * file but not when calling slice.
   */
  struct Mercator_map {

  public:

    /**
     * \brief Constructor.
     * \param map_dir directory for Mercator genomes and map files
     */
    Mercator_map (const std::string& map_dir);

    /**
     * \brief Find the Mercator_interval for a particular genomic sequence.
     *
     * Require that the sequence is completely or partially contained
     * in interval.
     * Uses a binary search over all Mercator_interval objects.
     * \param strict interval must be completely contained within Mercator_interval
     * \return index of Mercator_interval in __intervals, or size() if it can't be found
     */
    size_t find_mercator_interval (const std::string& genome, const std::string& chromosome,
				   const unsigned start, const unsigned end,
				   const bool strict = true) const;

    /**
     * \brief Choose random Genomic_intervals from sequence in a Mercator_map.
     *
     * Chosen according to a uniform distribution over sequences:
     * - choose a Mercator_interval according to the length of sequence it has for genome
     * - choose a point within the Mercator_interval
     * - use this point as a centroid for a Genomic_interval
     *   (note that the Genomic_interval will NOT be truncated, even if it does not
     *   fit within the containing Mercator_interval)
     * \param genome genome to select intervals for
     * \param length length of desired intervals
     * \param num number of intervals desired
     */
    std::vector<Genomic_interval> choose_random_intervals (const std::string& genome,
							   const size_t length, const size_t num) const;

    /**
     * \brief Get a Mercator_interval by its index.
     * \see find_mercator_interval
     */
    const Mercator_interval& get_interval (const size_t idx) const {
      return __intervals[idx];
    }

    /**
     * \brief Number of intervals.
     * 
     * This is equal to (# of bins) * (# of genomes).
     */
    size_t size() const { return __intervals.size(); }

    /**
     * \brief Number of bins.
     */
    size_t bins() const { return __num_bins; }

    /**
     * \brief Get all Mercator_intervals for a particular genome.
     */
    std::vector<Mercator_interval> get_mercator_intervals (const std::string& genome) const;

    /**
     * \brief Do we have information for a particular genome?
     */
    bool exists_genome (const std::string& genome) const {
      return __genome_index.find (genome) != __genome_index.end();
    }

    /**
     * \brief Get the list of genomes in the Mercator file format.
     */
    std::string get_genomes_string() const;

    /**
     * \brief Get the Mercator_interval for a particular bin and genome.     
     */
    inline const Mercator_interval& get_interval_by_bin (const unsigned bin, const std::string& genome) const;


  protected:

    /**
     * \brief Get indices in __intervals for all Mercator_interval objects for a particular bin.
     */
    inline const std::vector<size_t>& get_intervals_by_bin (const unsigned bin) const;

    /**
     * \brief Get iterator to start of __intervals.
     */
    std::vector<Mercator_interval>::const_iterator begin() const {
      return __intervals.begin();
    }

    /**
     * \brief Get iterator to end of __intervals.
     */
    std::vector<Mercator_interval>::const_iterator end() const {
      return __intervals.end();
    }


  private:

    /**
     * \brief Initialize from Mercator file.
     * \see parse_line
     */
    void from_file (const std::string& genomes_filename,
		    const std::string& map_filename);

    /**
     * \brief Initialize map from bins to Mercator_interval objects.
     *
     * Must be called after from_file.
     */
    void init_bin_index();

    /**
     * \brief Parse single line of Mercator map file.
     */
    void parse_line (const std::string& line);


    std::vector<std::string> __genomes;             ///< map from numeric genome index to name (ordered list of genomes)
    std::map<std::string, size_t> __genome_index;   ///< map from a genome to its numeric index in the map file

    size_t __num_bins;                              ///< number of bins in the map
    std::vector<Mercator_interval> __intervals;     ///< Mercator_interval entries

    /**
     * \brief Map from a bin and genome to the index for the corresponding Mercator_interval objects in __intervals.
     * The appropriate index in __intervals is accessed as __bin_index[bin][__genome_index[genome]].
     * If there is no such index (i.e., for a particular bin there is no associated interval for some genome,
     * then the index will be the nonsense index size().
     */
    std::vector<std::vector<size_t> > __bin_index;

  private:

    static const std::string mercator_empty_interval; ///< string designating no information for a particular genome in a map file

  };

  /**
   * \brief 
   *
   * Assumes a DNA alphabet, but DOES NOT check or enforce this!
   */
  struct Mercator_alignment : public Mercator_map {

  public:

    /**
     * \brief Constructor.
     * \param map_dir directory for Mercator genomes and map files
     * \param alignments_dir directory for alignment files
     */
    Mercator_alignment (const std::string& map_dir,
			const std::string& alignments_dir);

    /**
     * \brief Find and read the alignment file for a particular bin.
     * \param bin Mercator bin
     * \param stockholm_bin read alignment into this
     */
    void read_bin_alignment (const unsigned bin,
			     Stockholm& stockholm_bin,
			     const bool verbose = true) const;

    /**
     * \brief Get the subalignment for a particular genomic sequence.
     *
     * If requested, will return a truncated subalignment if only
     * part of the sequence can be mapped (rather than returning nothing).
     * If the strand is -, then the returned subalignment
     * will be reverse-complemented.
     * \param seq_db_subalign Sequence_database hold subalignment sequences
     * \param lazy warn, instead of die, if the interval can't be mapped
     * \param truncate_ok if the sequence can't all be mapped, then return a truncated subalignment
     * \param annotate_homology_information mark up alignment with homology information for sequences (in GFF format)
     * \return subalignment, or empty alignment if can't be found
     * \see subalignment
     */
    Stockholm* slice (Sequence_database& seq_db_subalign,
		      const std::string& genome, const std::string& chromosome,
		      char strand,
		      const unsigned start, const unsigned end,
		      const bool lazy = false, const bool truncate_ok = false,
		      const bool annotate_homology_information = false,
		      const bool verbose = true) const;

    /**
     * \brief Get the subalignment for a particular genomic sequence.
     *
     * If requested, will return a truncated subalignment if only
     * part of the sequence can be mapped (rather than returning nothing).
     * If the strand is -, then the returned subalignment
     * will be reverse-complemented (i.e., if the requested strand differs from
     * the strand specified in the Mercator mapping, then the alignment will 
     * be reverse-complemented w.r.t. the Mercator multiple alignment).
     * If the strand is unknown, then the + strand is assumed.
     * Looks for an alignment of the form 'bin.{stock,stk,mfa,fasta,fa,fas}'
     * in __alignments_dir.
     * This method populates Sequence_database and Stockholm objects
     * for both the entire bin alignment as well as the requested subalignment.
     * Saving on overhead associated with reading in the bin alignment,
     * this allows for efficient extraction of multiple sequences
     * from the same bin.
     * \param subalign_seq_db hold subalignment sequences
     * \param bin Mercator bin which stockholm_bin holds the alignment for (set to 0 for dummy/unknown value)
     * \param stockholm_bin hold Stockholm alignment for bin
     * \param lazy warn, instead of die, if the interval can't be mapped
     * \param truncate if the sequence can't all be mapped, then return a truncated subalignment
     * \param annotate_homology_information mark up alignment with homology information for sequences (in GFF format)
     * \see map_genomic_to_align
     * \return subalignment, or empty alignment if can't be found
     */
    Stockholm* slice (Sequence_database& seq_db_subalign,
		      unsigned& bin,
		      Stockholm& stockholm_bin,
		      const std::string& genome, const std::string& chromosome,
		      char strand,
		      const unsigned start, const unsigned end,
		      const bool lazy = false, const bool truncate_ok = false,
		      const bool annotate_homology_information = false,
		      const bool verbose = true) const;

    /**
     * \brief Find the homologous interval according to the Mercator mapping and corresponding multiple alignment.
     *
     * The returned Mercator_interval DOES NOT correspond to an actual entry
     * in the Mercator map file, but rather a sub-interval thereof.
     * Maps genome_from -> genome_to.
     * As always, coordinates are 0-based and fully-closed.
     * \param genome_from genome with the specified sequence interval
     * \param genome_to genome to find the homologous interval in
     * \param lazy warn, instead of die, if the interval can't be mapped
     * \param truncate if the sequence can't all be mapped, then use a truncated subalignment
     * \return Mercator_interval specifying the homologous interval; returned value is empty if it can't be mapped or if empty subsequence (all gaps in alignment)
     */
    Genomic_interval find_homologous_interval
      (const std::string& genome_from, const std::string& chromosome_from,
       const unsigned start_from, const unsigned end_from,
       const std::string& genome_to,
       const bool lazy = false, const bool truncate_ok = false,
       const bool verbose = true) const;

    /**
     * \brief Find the homologous interval according to the Mercator mapping and corresponding multiple alignment.
     *
     * The returned Mercator_interval DOES NOT correspond to an actual entry
     * in the Mercator map file, but rather a sub-interval thereof.
     * Maps genome_from -> genome_to.
     * This method populates Sequence_database and Stockholm objects
     * for the entire bin alignment.
     * Saving on overhead associated with reading in the bin alignment,
     * this allows for efficient multiple queries from the same bin.
     * As always, coordinates are 0-based and fully-closed.
     * \param genome_from genome with the specified sequence interval
     * \param genome_to genome to find the homologous interval in
     * \param lazy warn, instead of die, if the interval can't be mapped
     * \param truncate if the sequence can't all be mapped, then use a truncated subalignment
     * \return Mercator_interval specifying the homologous interval; returned value is empty if it can't be mapped or if empty subsequence (all gaps in alignment)
     */
    Genomic_interval find_homologous_interval
      (unsigned& bin,
       Stockholm& stockholm_bin,
       const std::string& genome_from, const std::string& chromosome_from,
       const unsigned start_from, const unsigned end_from,
       const std::string& genome_to,
       const bool lazy = false, const bool truncate_ok = false,
       const bool verbose = true) const;

    /**
     * \brief Map the features in a GFF database to another genome.
     *
     * \param lazy warn, instead of die, if the interval can't be mapped
     * \param truncate_ok if the sequence can't all be mapped, then return a truncated subalignment
     * \param force_entry if a feature can't be mapped, then add an empty entry to the GFF file (rather than skipping it entirely)
     */
    GFF_database map_gff_database
      (const std::string& genome_from, const std::string& genome_to,
       const GFF_database& gff_db_from,
       const bool lazy = false, const bool truncate_ok = false,
       const bool verbose = true,
       const bool force_entry = false) const;

    /**
     * \brief Map alignment coordinates (w.r.t. the alignment of a particular Mercator bin)
     * to genomic coordinates for a genome.
     *
     * \param bin bin which we have the alignment for
     * \param stockholm_bin the alignment of the bin
     * \param genome genome to map coordinates for
     * \param start 0-based column start coordinate
     * \param end 0-based column end coordinate
     * \param rc has the passed alignment of the bin been reverse-complemented?  (if so,
     * then appropriately deals with the necessary coordinate conversions)
     * \return Genomic_interval representing the appropriate genomic sequence for genome,
     * with nonsense coordinates [1, 0] if the interval is empty (no ungapped sequence).
     */
    Genomic_interval map_align_to_genomic
      (const unsigned bin,
       const Stockholm& stockholm_bin,
       const std::string& genome,
       const unsigned start, const unsigned end,
       const bool rc = false) const;

    /**
     * \brief Map genomic coordinates to alignment coordinates for a particular bin.
     *
     * \param start 0-based genomic start coordinate
     * \param end 0-based genomic end coordinate
     * \param bin Mercator bin which stockholm holds the alignment for (set to 0 for dummy/unknown value)
     if reading in a new alignment is necessary, then this function sets bin appropriately
     * \param stockholm_bin hold Stockholm alignment for bin
     * \return alignment column coordinates or nonsense coordinates [1, 0] if no mapping
     */
    Interval map_genomic_to_align
      (const std::string& genome, const std::string& chromosome,
       const unsigned start, const unsigned end,
       unsigned& bin,
       Stockholm& stockholm_bin,
       const bool lazy = false, const bool truncate_ok = false,
       const bool verbose = true) const;

    /**
     * \brief Map genomic coordinates to sequence coordinates within a particular bin.
     * \param start 0-based genomic start coordinate
     * \param end 0-based genomic end coordinate
     * \param bin Mercator bin which stockholm holds the alignment for (set to 0 for dummy/unknown value)
     if reading in a new alignment is necessary, then this function sets bin appropriately
     * \param stockholm_bin hold Stockholm alignment for bin
     * \return alignment column coordinates or nonsense coordinates [1, 0] if no mapping
     */
    Interval map_genomic_to_bin_seq
      (const std::string& genome, const std::string& chromosome,
       const unsigned start, const unsigned end,
       unsigned& bin,
       Stockholm& stockholm_bin,
       const bool lazy = false, const bool truncate_ok = false,
       const bool verbose = true) const;
    

  private:

    std::string __alignments_dir;                   ///< directory to search for alignments
    static const Alphabet alphabet;                 ///< DNA alphabet
    std::vector<std::string> __suffixes_mfa;        ///< file suffixes to search for alignments
    std::vector<std::string> __suffixes_stockholm;  ///< file suffixes to search for alignments

  };

  inline const std::vector<size_t>& Mercator_map::get_intervals_by_bin (const unsigned bin) const {
    return __bin_index[bin];
  }

  inline const Mercator_interval& Mercator_map::get_interval_by_bin (const unsigned bin, const std::string& genome) const {

    // check sane
    assert (exists_genome (genome));
    assert (bin <= bins()); // remember that the numbering of Mercator bins is 1-based
                            // and that we preserve the 1-based indexing in __bin_index
                            // (see init_bin_index)

    const size_t index = __bin_index[bin][__genome_index.find (genome)->second];

    // make sure that we have a valid Mercator_interval for this (bin, genome) pair
    assert (index < size());

    // make sure that we found it
    assert ((get_interval (index).bin == bin) && (get_interval (index).genome == genome));

    return get_interval (index);

  }

}

#endif /* SEQ_MERCATOR_INCLUDED */

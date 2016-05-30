
/**
 * \file mercator.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <algorithm>
#include <numeric>

#include "seq/mercator.h"
#include "seq/gff.h"

using namespace fsa;

const std::string Mercator_map::mercator_empty_interval = "NA";
const Alphabet Mercator_alignment::alphabet = DNA_alphabet();


GFF_database Genomic_interval::convert_to_gff_db (const std::vector<Genomic_interval>& intervals) {

  GFF_database gff_db;
  for (std::vector<Genomic_interval>::const_iterator interval = intervals.begin(); interval != intervals.end(); ++interval)
    gff_db.store_entry (GFF (interval->chromosome, interval->start + 1, interval->end + 1));  // convert to 1-based coordinates

  return gff_db;

}

Mercator_map::Mercator_map (const std::string& map_dir)
  : __num_bins (0) {

  from_file (map_dir + "/genomes", map_dir + "/map");
  init_bin_index(); // this MUST only be done after sorting __intervals...

}

void Mercator_map::from_file (const std::string& genomes_filename,
			      const std::string& map_filename) {

  cerr << "Reading Mercator map '" << map_filename << "'...";

  std::ifstream filestream;

  // read genomes file
  filestream.open (genomes_filename.c_str(), std::ifstream::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << genomes_filename << "' for reading." << endl;
    exit (1);
  }

  std::string line;
  while (!filestream.eof()) {
    getline (filestream, line);

    std::stringstream ss (line);
    std::string buffer;
    std::vector<std::string> tokens; // hold whitespace-separated tokens
    while (ss >> buffer)
      tokens.push_back (buffer);
    
    for (size_t i = 0; i < tokens.size(); ++i) {
      __genomes.push_back (tokens[i]);
      __genome_index.insert (std::make_pair (tokens[i], i));
    }

  }
  filestream.close();

  // read map file
  filestream.open (map_filename.c_str(), std::ifstream::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << map_filename << "' for reading." << endl;
    exit (1);
  }

  while (!filestream.eof()) {
    getline (filestream, line);
    Util::chomp (line);
    if (!line.length()) { continue; }
    parse_line (line);
  }
  filestream.close();

  // sort entries!  otherwise binary search will fail
  std::sort (__intervals.begin(), __intervals.end(),
	     Genomic_interval::Genomic_interval_less());

  cerr << "done." << endl;

}

void Mercator_map::init_bin_index() {
  
  // __bin_index[bin][__genome_index[genome]] is the index in __intervals
  // for a particular (bin, genome) pair

  // initialize stuff
  // add a dummy entry for the 0th bin (hence bins() + 1)
  // so that we can index with the actual Mercator bin number
  // set everything to the nonsense index size() at first
  __bin_index.assign (bins() + 1, std::vector<size_t> (__genome_index.size(), size()));

  // store indices for each Mercator_interval
  for (size_t i = 0; i < size(); ++i) {
    const Mercator_interval& interval = get_interval (i);
    assert (interval.bin < __bin_index.size());
    assert (__bin_index[interval.bin].size() > __genome_index[interval.genome]);
    __bin_index[interval.bin][__genome_index[interval.genome]] = i;
  }
  
}

void Mercator_map::parse_line (const std::string& line) {

  std::stringstream ss (line);
  std::string buffer;
  std::vector<std::string> tokens; // hold whitespace-separated tokens
  while (ss >> buffer)
    tokens.push_back (buffer);

  // sanity checks
  if (tokens.size() != 4 * __genomes.size() + 1) {
    cerr << "WARNING: Can't parse line: '" << line << "'; skipping." << endl;
    return;
  }

  // extract bin number for this line of the map file  
  const unsigned bin = static_cast<unsigned> (atoi (tokens[0].c_str()));

  // keep track of how many bins there are
  // (remember that bins use 1-based numbering)
  if (bin > __num_bins)
    __num_bins = bin;

  for (size_t i = 1; i < tokens.size(); i += 4) {

    // if this is an empty interval,
    // meaning that this genome has no orthologous sequence in this bin,
    // then store nothing for this genome
    if (tokens[i] == mercator_empty_interval)
      continue;

    const std::string& genome = __genomes[i / 4];
    std::string& chromosome = tokens[i];
    Util::strip_leading_chr (chromosome);
    const unsigned start = atoi (tokens[i + 1].c_str());
    const unsigned end = atoi (tokens[i + 2].c_str()) - 1; // (remember that Mercator coordinates are half-open, [start, end); convert to fully-closed here)
    const char strand = (tokens[i + 3])[0];

    assert (end >= start);
    // strand sanity checks
    assert (tokens[i + 3].size() == 1);
    if (strand != '+' && strand != '-') {
      cerr << "ERROR: Unrecognized strand " << strand << "." << endl;
      exit (1);
    }

    // store interval plus indexing from bin -> interval
    __intervals.push_back (Mercator_interval (bin,
					      genome, chromosome,
					      start, end, strand));

  }
  

}

std::vector<Mercator_interval> Mercator_map::get_mercator_intervals (const std::string& genome) const {

  std::vector<Mercator_interval> gintervals;
  for (std::vector<Mercator_interval>::const_iterator interval = __intervals.begin(); interval != __intervals.end(); ++interval) {
    if (interval->genome == genome)
      gintervals.push_back (*interval);
  }

  return gintervals;

}

std::string Mercator_map::get_genomes_string() const {

  return Util::join (__genomes, "\t");

}

size_t Mercator_map::find_mercator_interval (const std::string& genome, const std::string& chromosome,
					     const unsigned start, const unsigned end,
					     const bool strict /* = true */) const {

  // initialize dummy object for comparison
  Mercator_interval m;
  m.genome = genome;
  m.chromosome = chromosome;
  Util::strip_leading_chr (m.chromosome);
  m.start = start;
  m.end = end;

  std::vector<Mercator_interval>::const_iterator mercator = std::upper_bound (__intervals.begin(),
									      __intervals.end(),
									      m,
									      Genomic_interval::Genomic_interval_less());

  // decrement element to try to get to the Mercator object
  // which (possibly) contains data for pos
  // (if element is at beginning, then we have no complete map for sequence,
  // but we may have a partial map)
  if (mercator != __intervals.begin())
    --mercator;

  // now check that we're on the correct genome & chromosome
  if (mercator->genome != m.genome || mercator->chromosome != m.chromosome)
    return __intervals.size();

  bool mapped = false;

  // does the complete interval [start, end] fall outside the found mapping?
  if (start > mercator->end || end < mercator->start)
    mapped = false;

  // does the complete interval [start, end] fall into the found mapping?
  else if (start >= mercator->start && end <= mercator->end)
    mapped = true;

  // if we don't require complete coverage, then look for partial coverage
  else if (!strict) {
    if ((start <= mercator->start && end >= mercator->end)
	|| (start >= mercator->start)
	|| (end <= mercator->end))
      mapped = true;
  }

  if (mapped)
    return mercator - __intervals.begin();

  // else couldn't map successfully
  return __intervals.size();

}

std::vector<Genomic_interval> Mercator_map::choose_random_intervals (const std::string& genome,
								     const size_t length, const size_t num) const {

  // create lists of all intervals for the requested genome
  std::vector<double> interval_lengths; // weight intervals for genome according to the sequence contained therein
  std::vector<size_t> interval_indices; // indices for intervals in *this
  for (size_t i = 0; i < size(); ++i) {
    if (get_interval (i).genome == genome) {
      interval_lengths.push_back (get_interval (i).length());
      interval_indices.push_back (i);
    }
  }

  // now normalize properly to get a probability distribution
  const double norm = accumulate (interval_lengths.begin(), interval_lengths.end(),
				  0.);
  for (size_t i = 0; i < interval_lengths.size(); ++i)
    interval_lengths[i] /= norm;

  // choose intervals!
  std::vector<Genomic_interval> random_intervals;
  while (random_intervals.size() < num) {

    // choose a Mercator_interval according to the length of sequence it has for genome
    const Mercator_interval& interval = get_interval (interval_indices[Util::choose_from_distribution (interval_lengths)]);
    
    // choose a point within the Mercator_interval
    const unsigned centroid = interval.start + Util::rand (interval.length() - 1);

    unsigned random_start = centroid > static_cast<size_t> (length / 2) ? centroid - static_cast<size_t> (length / 2) : 0;
    if (random_start < interval.start)
      random_start = interval.start;
    const unsigned random_end = centroid + static_cast<size_t> (length / 2);

    // use this point as a centroid for a Genomic_interval
    random_intervals.push_back (Genomic_interval (genome, interval.chromosome,
						  random_start, random_end, interval.strand));

  }

  return random_intervals;

}

Mercator_alignment::Mercator_alignment (const std::string& map_dir,
					const std::string& alignments_dir)
  : Mercator_map (map_dir),
    __alignments_dir (alignments_dir) {

  __suffixes_stockholm.push_back ("stock");
  __suffixes_stockholm.push_back ("stk");
  __suffixes_mfa.push_back ("mfa");
  __suffixes_mfa.push_back ("fasta");
  __suffixes_mfa.push_back ("fa");
  __suffixes_mfa.push_back ("fas");

}

void Mercator_alignment::read_bin_alignment (const unsigned bin,
					     Stockholm& stockholm_bin,
					     const bool verbose /* = true */) const {

  std::string filename;

  // look for Stockholm files
  for (std::vector<std::string>::const_iterator suffix = __suffixes_stockholm.begin(); suffix != __suffixes_stockholm.end(); ++suffix) {
    filename = __alignments_dir + '/' + Util::to_string (bin) + '.' + *suffix;
    if (Util::exists_file (filename)) {
      stockholm_bin.read_stockholm (filename, true, // strict = true
				    verbose);
      return;
    }
  }

  // look for MFA files
  for (std::vector<std::string>::const_iterator suffix = __suffixes_mfa.begin(); suffix != __suffixes_mfa.end(); ++suffix) {
    filename = __alignments_dir + '/' + Util::to_string (bin) + '.' + *suffix;
    if (Util::exists_file (filename)) {
      stockholm_bin.read_mfa (filename, true, // strict = true
			      verbose);
      return;
    }
  }

  // if we're here, then we didn't locate an alignment file
  cerr << "ERROR: Couldn't locate alignment file of the form '"
       << bin << ".{stock,stk,mfa,fasta,fa,fas}' in directory '" << __alignments_dir << "/'." << endl;
  exit (1);

}

Stockholm* Mercator_alignment::slice
(Sequence_database& seq_db_subalign,
 const std::string& genome, const std::string& chromosome,
 char strand,
 const unsigned start, const unsigned end,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool annotate_homology_information /* = false */,
 const bool verbose /* = true */) const {


  Sequence_database seq_db_bin;
  Stockholm stockholm_bin (seq_db_bin);
  unsigned bin;

  return slice (seq_db_subalign,
		bin,
		stockholm_bin,
		genome, chromosome,
		strand,
		start, end,
		lazy, truncate_ok,
		annotate_homology_information,
		verbose);

}

Stockholm* Mercator_alignment::slice
(Sequence_database& seq_db_subalign,
 unsigned& bin,
 Stockholm& stockholm_bin,
 const std::string& genome, const std::string& chromosome,
 char strand,
 const unsigned start, const unsigned end,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool annotate_homology_information /* = false */,
 const bool verbose /* = true */) const {

  // sanity check on strand
  if (strand != '+' && strand != '-' && verbose) {
    cerr << "WARNING: Unknown strand '" << strand << "'; assuming + strand." << endl;
    strand = '+';
  }

  // map from genomic to alignment coordinates
  const Interval align_interval = map_genomic_to_align (genome, chromosome,
							start, end,
							bin,
							stockholm_bin,
							lazy, truncate_ok,
							verbose);
  // catch the case of an empty interval, meaning that no mapping could be found
  if (align_interval.length() == 0) {
    seq_db_subalign.clear();
    return new Stockholm (seq_db_subalign);
  }  

  const unsigned start_align = align_interval.start;
  const unsigned end_align = align_interval.end;

  // pull out the subalignment
  Stockholm* subalignment = stockholm_bin.subalignment (seq_db_subalign,
							start_align, end_align);

  // revcomp if necessary
  const size_t idx = find_mercator_interval (genome, chromosome,
					     start, end,
					     !truncate_ok); // not strict if truncate_ok
  const Mercator_interval& interval = get_interval (idx);
  if (strand != interval.strand)
    subalignment->revcomp (Mercator_alignment::alphabet);

  // show homology information for all seqs if desired
  // (i.e., what subsequences are in the subalignment)
  if (annotate_homology_information) {

    // has the subaligment been reverse-complemented w.r.t. the raw Mercator mapping?
    const bool was_rc = (strand != interval.strand);

    for (size_t r = 0; r < stockholm_bin.rows(); ++r) {

      const std::string& genome2 = stockholm_bin.get_row_name (r);
      
      GFF gff;

      // if "source" genome, annotate directly
      if (genome2 == genome) {

	gff.seqid = chromosome;
	Util::strip_leading_chr (gff.seqid); // for consistency with rest of sequences
	gff.start = start + 1;               // convert to 1-based coordinates (b/c GFF)
	gff.end = end + 1;
	gff.strand = strand;

      }

      // else find the homologous sequence to annotate with
      else {

	// find homologous interval
	Genomic_interval interval2 = find_homologous_interval (bin,
							       stockholm_bin,
							       genome, chromosome,
							       start, end,
							       genome2,
							       lazy, truncate_ok,
							       verbose);

	// figure out what strand genome2 is on in the subalignment
	if (was_rc)
	  Sequence::complement_strand (interval2.strand);

	// annotate
	gff.seqid = interval2.chromosome;
	gff.start = interval2.start + 1; // convert to 1-based coordinates (b/c GFF)
	gff.end = interval2.end + 1;
	gff.strand = interval2.strand;

      }
      
      // annotate!
      subalignment->set_gs_annot (genome2, Stockholm::gff_annotation, gff.to_string());

    }

  }

  return subalignment;

}


Genomic_interval Mercator_alignment::find_homologous_interval
(const std::string& genome_from, const std::string& chromosome_from,
 const unsigned start_from, const unsigned end_from,
 const std::string& genome_to,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool verbose /* = true */) const {

  Sequence_database seq_db_bin;
  Stockholm stockholm_bin (seq_db_bin);
  unsigned bin;

  return find_homologous_interval (bin,
				   stockholm_bin,
				   genome_from, chromosome_from,
				   start_from, end_from,
				   genome_to,
				   lazy, truncate_ok,
				   verbose);  

}

Genomic_interval Mercator_alignment::find_homologous_interval 
(unsigned& bin,
 Stockholm& stockholm_bin,
 const std::string& genome_from, const std::string& chromosome_from,
 const unsigned start_from, const unsigned end_from,
 const std::string& genome_to,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool verbose /* = true */) const {

  // check that these are valid genomes
  if (!exists_genome (genome_from)) {
    cerr << "ERROR: No information for genome '" << genome_from << "'; valid genomes are:" << endl
	 << get_genomes_string() << endl;
    exit (1);
  }
  if (!exists_genome (genome_to)) {
    cerr << "ERROR: No information for genome '" << genome_to << "'; valid genomes are:" << endl
	 << get_genomes_string() << endl;
    exit (1);
  }

  // try to map the interval
  const size_t idx_from = find_mercator_interval (genome_from, chromosome_from,
						  start_from, end_from,
						  !truncate_ok);
  if (idx_from == size()) {
    if (!lazy) {
      cerr << "ERROR: No mapping exists for: " << genome_from << '\t' << chromosome_from << '\t'
	   << start_from << '\t' << end_from << endl;
      exit (1);
    } else {
      cerr << "WARNING: No mapping exists for: " << genome_from << '\t' << chromosome_from << '\t'
	   << start_from << '\t' << end_from << endl;
      return Genomic_interval();
    }
  }
  const Mercator_interval& interval_from = get_interval (idx_from);
  const unsigned this_bin = interval_from.bin;

  // find the bin for genome_map_to
  Mercator_interval mercator_interval_to = get_interval_by_bin (this_bin, genome_to);

  // now extract the homologous interval using the base-level alignment

  // if this_bin isn't the passed alignment,
  // then find and open the appropriate bin
  if (bin != this_bin) {
    // find & read alignment
    read_bin_alignment (this_bin, stockholm_bin, verbose);
    // record current bin
    bin = this_bin;
  }


  // map the sequence interval using the alignment
  // (remember to convert coordinates w.r.t. Mercator_interval)

  // first calculate the sequence coordinate offsets w.r.t. the interval_from alignment
  // the conditionals truncate the requested sequence as appropriate if it's not entirely
  // contained within the interval_from alignment
  unsigned start_from_seq, end_from_seq;
  // if Mercator_interval is on - strand, then need to flip the interval and calculate w.r.t.
  // the end of the mapped Mercator_interval (because the sequence in the alignment has
  // been reverse-complemented, so it "begins" at the end of the Mercator_interval)
  if (interval_from.strand == '-') {
    start_from_seq = static_cast<int> (interval_from.end - end_from) >= 0
      ? interval_from.end - end_from
      : 0;
    end_from_seq = static_cast<int> (interval_from.end - start_from) >= 0
      ? interval_from.end - start_from
      : interval_from.end - interval_from.start;
  }
  // if Mercator_interval is on + strand, then just subtract offsets
  else {
    start_from_seq = static_cast<int> (start_from - interval_from.start) >= 0
      ? start_from - interval_from.start
      : 0;
    end_from_seq = static_cast<int> (interval_from.end - end_from) >= 0
      ? end_from - interval_from.start
      : interval_from.end - interval_from.start;
  }

  // warn if we've truncated
  if (static_cast<int> (start_from - interval_from.start) < 0 || static_cast<int> (interval_from.end - end_from) < 0
      || static_cast<int> (interval_from.end - end_from) < 0 || static_cast<int> (interval_from.end - start_from) < 0)
    cerr << "WARNING: Truncating alignment slice to mapped sequence: " << genome_from << '\t' << chromosome_from << '\t'
	 << start_from << '\t' << end_from << '\t' << interval_from.strand << endl;

  // now map to coordinates in the other sequence
  const Interval interval_to = stockholm_bin.map_seq_to_seq (genome_from, start_from_seq, end_from_seq,
							     genome_to);

  // check for the case of no sequence (empty interval) in the "to" genome subsequence
  // Alignment::map_seq_to_seq returns (1, 0) if no sequence
  if (interval_to.end < interval_to.start)
    return Genomic_interval();

  // now convert mercator_interval_to to hold the base-level mapping
  // else add offsets w.r.t. end of the Mercator_interval
  if (mercator_interval_to.strand == '-') {
    mercator_interval_to.start = mercator_interval_to.end - interval_to.end;
    mercator_interval_to.end = mercator_interval_to.end - interval_to.start;
  }
  // if Mercator_interval is on + strand, then just add offsets back on
  else {
    mercator_interval_to.end = mercator_interval_to.start + interval_to.end;
    mercator_interval_to.start += interval_to.start;
  }
  // (NB: both of these transformations are derived by 
  // inverting the ones used above to set {start,end}_from_seq)

  return mercator_interval_to;

}

GFF_database Mercator_alignment::map_gff_database
(const std::string& genome_from, const std::string& genome_to,
 const GFF_database& gff_db_from,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool verbose /* = true */,
 const bool force_entry /* = false */) const {

  GFF_database gff_db_to;

  // initialize persistent sequence & alignment information for current bin
  Sequence_database seq_db_bin;
  Stockholm stockholm_bin (seq_db_bin);
  unsigned bin = 0; // initialize to dummy 0 value

  // iterate through GFF
  for (std::vector<GFF>::const_iterator gff_from = gff_db_from.begin(); gff_from != gff_db_from.end(); ++gff_from) {
    
    // use the alignment to find the homologous interval
    Genomic_interval region_target = find_homologous_interval (bin,
							       stockholm_bin,
							       genome_from, gff_from->seqid,
							       gff_from->start - 1, gff_from->end - 1, // convert to 0-based coordinates
							       genome_to,
							       lazy, truncate_ok,
							       verbose);

    // check that it was mapped successfully
    // (Mercator_alignment::find_homologous_interval returns an empty value
    // if it wasn't or if homologous subsequence is empty)
    // NB no need to enforce lazy here, b/c find_homologous_interval does this for us
    if (region_target.genome == "") {
      cerr << "WARNING: No homologous sequence found (all gaps) in " << genome_to << " for: " << genome_from << '\t' << gff_from->seqid << '\t'
	   << gff_from->start - 1 << '\t' << gff_from->end - 1 << '\t' << gff_from->strand << endl;
      // store an empty entry if requested
      if (force_entry) {
	GFF gff_to;
	gff_to.source = std::string (PACKAGE) + '\\' + gff_from->source;
	gff_to.type = gff_from->type;
	gff_to.start = 1; // nonsense values for start and end: [1, 0]
	gff_to.end = 0;
	gff_db_to.store_entry (gff_to);
      }
      continue;
    }

    // get the Mercator_interval for the source
    const size_t idx_source = find_mercator_interval (genome_from, gff_from->seqid,
						      gff_from->start - 1, gff_from->end - 1, // convert to 0-based coordinates
						      !truncate_ok);
    assert (idx_source < size());
    const Mercator_interval& interval_source = get_interval (idx_source);

    // use the source interval to determine whether the source feature was 
    // reverse-complemented w.r.t. the strand in the Mercator_interval
    // if it was, then the target feature will also be reverse-complemented
    // w.r.t. the strand in the Mercator_interval
    const bool is_rc = (interval_source.strand != gff_from->strand);

    GFF gff_to;
    gff_to.seqid = region_target.chromosome;
    gff_to.source = std::string (PACKAGE) + '\\' +  gff_from->source;
    gff_to.type = gff_from->type;
    gff_to.start = region_target.start + 1; // convert back to 1-based coordinates
    gff_to.end = region_target.end + 1;
    if (gff_from->strand == GFF::unknown_strand)   // if the source feature had no annotated strand, then do the same for the target feature
      gff_to.strand = GFF::unknown_strand;
    else {                           // otherwise set the target strand appropriately using the homology mapping implied by the alignment
      gff_to.strand = region_target.strand;
      if (is_rc)
	Sequence::complement_strand (gff_to.strand);
    }
    for (std::vector<std::string>::const_iterator key = gff_from->attributes_ordering.begin(); key != gff_from->attributes_ordering.end(); ++key) {
      const std::vector<std::string>& values = gff_from->attributes_map.find (*key)->second;
      for (std::vector<std::string>::const_iterator value = values.begin(); value != values.end(); ++value) {
	// watch for case of quoted fields
	if ((*value)[0] == '\"' && (*value)[value->length() - 1] == '\"') {
	  gff_to.add_value (*key, '\"' + genome_to + '\\'+ std::string (*value, 1, value->length() - 1));
	} else if ((*value)[0] == '\'' && (*value)[value->length() - 1] == '\'') {
	  gff_to.add_value (*key, '\'' + genome_to + '\\'+ std::string (*value, 1, value->length() - 1));
	} else {
	  gff_to.add_value (*key, genome_to + '\\'+ *value);
	}
      }
    }

    gff_db_to.store_entry (gff_to);

  }

  return gff_db_to;

}

Genomic_interval Mercator_alignment::map_align_to_genomic
(const unsigned bin,
 const Stockholm& stockholm_bin,
 const std::string& genome,
 const unsigned start, const unsigned end,
 const bool rc /* = false */) const {

  // check sane
  assert (start <= end);
  assert (end < stockholm_bin.columns());
  assert (exists_genome (genome));

  // find the Mercator_interval describing the sequence for genome in this bin
  // check whether the passed alignment of the bin sequence has been reverse-complemented:
  // if so, then we need to flip the strand of Mercator_interval for this bin
  Mercator_interval mercator_interval = get_interval_by_bin (bin, genome);
  if (rc)
    Sequence::complement_strand (mercator_interval.strand);

  // map column indices to sequence coordinates (w.r.t. stockholm_bin)
  const Interval align_seq_coords = stockholm_bin.map_align_to_seq (genome,
								    start, end);

  // now map from sequence coordinates within the alignment of the bin to genomic coordinates
  Genomic_interval genomic_interval;
  genomic_interval.genome = genome;
  genomic_interval.chromosome = mercator_interval.chromosome;
  genomic_interval.strand = mercator_interval.strand;

  // check for the case of an empty interval
  if (align_seq_coords.start > align_seq_coords.end) {
    genomic_interval.start = 1;
    genomic_interval.end = 0;
    return genomic_interval;
  }

  // tricky if on the - strand
  // (sequence has been reversed, so we need to count from right to left)
  if (genomic_interval.strand == '-') {
    genomic_interval.start = mercator_interval.end - align_seq_coords.end;
    genomic_interval.end = mercator_interval.end - align_seq_coords.start;
  }
  // ...easier if aligned sequence is on the + strand
  else {
    genomic_interval.start = mercator_interval.start + align_seq_coords.start;
    genomic_interval.end = mercator_interval.start + align_seq_coords.end;
  }
  // Note that the above calculations for the - strand are equivalent to:
  //	genomic_coords.start = (seqlength - 1) - align_seq_coords.end + mercator_interval.start;
  //	genomic_coords.end = (seqlength - 1) - align_seq_coords.start + mercator_interval.start;
  // where seqlength is the ungapped sequence length of the bin.

  return genomic_interval;

}

Interval Mercator_alignment::map_genomic_to_align
(const std::string& genome, const std::string& chromosome,
 const unsigned start, const unsigned end,
 unsigned& bin,
 Stockholm& stockholm_bin,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool verbose /* = true */) const {

  // map to sequence coordinates within the bin
  const Interval bin_seq_coords = map_genomic_to_bin_seq (genome, chromosome,
							  start, end,
							  bin,
							  stockholm_bin,
							  lazy, truncate_ok,
							  verbose);

  // catch the case of a 0-lenth interval
  if (bin_seq_coords.length() == 0)
    return Interval (1, 0);

  // then map to alignment coordinates
  const unsigned start_align = stockholm_bin.map_seq_to_align (genome, bin_seq_coords.start);
  const unsigned end_align = stockholm_bin.map_seq_to_align (genome, bin_seq_coords.end);

  return Interval (start_align, end_align);

}

Interval Mercator_alignment::map_genomic_to_bin_seq
(const std::string& genome, const std::string& chromosome,
 const unsigned start, const unsigned end,
 unsigned& bin,
 Stockholm& stockholm_bin,
 const bool lazy /* = false */, const bool truncate_ok /* = false */,
 const bool verbose /* = true */) const {

  // check that this is a valid genome
  if (!exists_genome (genome)) {
    cerr << "ERROR: No information for genome '" << genome << "'; valid genomes are:" << endl
	 << get_genomes_string() << endl;
    exit (1);
  }

  // try to map the interval
  const size_t idx = find_mercator_interval (genome, chromosome,
					     start, end,
					     !truncate_ok); // not strict if truncate_ok
  if (idx == size()) {
    if (!lazy) {
      cerr << "ERROR: No mapping exists for: " << genome << '\t' << chromosome << '\t'
	   << start << '\t' << end << endl;
      exit (1);
    } else {
      cerr << "WARNING: No mapping exists for: " << genome << '\t' << chromosome << '\t'
	   << start << '\t' << end << endl;
      return Interval (1, 0);
    }
  }
  const Mercator_interval& interval = get_interval (idx);
  const unsigned this_bin = interval.bin;

  // if this_bin isn't the passed alignment,
  // then find and open the appropriate bin
  if (bin != this_bin) {
    // find & read alignment
    read_bin_alignment (this_bin, stockholm_bin, verbose);
    // record current bin
    bin = this_bin;
  }

  // map sequence to alignment coordinates

  // first calculate the sequence coordinates within the interval alignment
  // (i.e., subtract off the offsets for the mapped Mercator_interval)
  // the conditionals truncate the requested sequence as appropriate if it's not entirely
  // contained within the interval alignment
  unsigned start_seq, end_seq;

  // if Mercator_interval is on - strand, then need to flip the interval and calculate w.r.t.
  // the end of the mapped Mercator_interval (because the sequence in the alignment has
  // been reverse-complemented, so it "begins" at the end of the Mercator_interval)
  // (very tricky...)
  if (interval.strand == '-') {
    start_seq = static_cast<int> (interval.end - end) >= 0
      ? interval.end - end
      : 0;
    end_seq = static_cast<int> (interval.end - start) >= 0
      ? interval.end - start
      : interval.end - interval.start;
  }
  // if Mercator_interval is on + strand, then just subtract offsets
  else {
    start_seq = static_cast<int> (start - interval.start) >= 0
      ? start - interval.start
      : 0;
    end_seq = static_cast<int> (interval.end - end) >= 0
      ? end - interval.start
      : interval.end - interval.start;
  }
  
  // warn if we've truncated
  if (static_cast<int> (start - interval.start) < 0 || static_cast<int> (interval.end - end) < 0
      || static_cast<int> (interval.end - end) < 0 || static_cast<int> (interval.end - start) < 0)
    cerr << "WARNING: Truncating alignment slice to mapped sequence: " << genome << '\t' << chromosome << '\t'
	 << start << '\t' << end << endl;
  
  return Interval (start_seq, end_seq);

}

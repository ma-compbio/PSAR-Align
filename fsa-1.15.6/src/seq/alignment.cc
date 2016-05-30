
/**
 * \file alignment.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <numeric>
#include <functional>
#include <algorithm>

#include "util/regexp.h"
#include "seq/alignment.h"

using namespace fsa;

// NB: Who knows why, but I get linker errors when these definitions
// are in the header file.  It works fine when they're not.  Very strange!

const char Alignment::gap_char = '-';

const std::string Stockholm::gff_annotation = "GFF";
const std::string Stockholm::percent_id_annotation = "Percent_id";
const std::string Stockholm::gap_fraction_annotation = "Gap_fraction";

const std::string Stockholm::format_identifier = "STOCKHOLM";
const std::string Stockholm::version_identifier = "1.0";
const std::string Stockholm::alignment_header = "# " + Stockholm::format_identifier + ' ' + Stockholm::version_identifier;
const std::string Stockholm::alignment_separator = "//";

const std::string Stockholm::file_annotation = "#=GF";
const std::string Stockholm::column_annotation = "#=GC";
const std::string Stockholm::sequence_annotation = "#=GS";
const std::string Stockholm::sequence_column_annotation = "#=GR";

const char Stockholm::annotation_wildcard_char = '.';



Alignment_row::Alignment_row (Row_path* row_path)
  : __row_path (row_path) {

  // calculate the total sequence length
  __seqlength = static_cast<size_t> (accumulate (__row_path->begin(), __row_path->end(), 0));

  build_coordinate_indices();

}

Alignment_row::Alignment_row (const Alignment_row& parent) {

  // deep copy
  __seq_to_align_coords_map = parent.__seq_to_align_coords_map;
  __align_to_seq_coords_map = parent.__align_to_seq_coords_map;
  __row_path = new Row_path (*parent.__row_path);
  __seqlength = parent.__seqlength;

}

Alignment_row& Alignment_row::operator= (const Alignment_row& parent) {

  // check for self-assignment
  if (this == &parent)
    return *this;

  // clear old memory
  delete __row_path;

  // deep copy
  __seq_to_align_coords_map = parent.__seq_to_align_coords_map;
  __align_to_seq_coords_map = parent.__align_to_seq_coords_map;
  __row_path = new Row_path (*parent.__row_path);
  __seqlength = parent.__seqlength;

  return *this;

}

Alignment_row::~Alignment_row() {

  delete __row_path;

}

void Alignment_row::build_coordinate_indices() {

  __seq_to_align_coords_map.clear();
  __align_to_seq_coords_map.clear();

  // build map from sequence to alignment coordinates, as well as
  // map from alignment to sequence coordinates
  // if an alignment position is gapped, then store the sequence
  // position corresponding to the next ungapped alignment position
  // note that this means that if the alignment ends with a run of gaps,
  // then the gapped alignment positions will be mapped to __seqlength,
  // so this end case must be explicitly checked for
  __seq_to_align_coords_map.reserve (__seqlength);
  __align_to_seq_coords_map.reserve (__row_path->size());
  size_t seqpos = 0; // this always holds the "next" sequence position
  for (size_t alignpos = 0; alignpos < __row_path->size(); ++alignpos) {

    if (is_gapped (alignpos)) {
      // do nothing for __seq_to_align_coords_map
      // record previous sequence position for __align_to_seq_coords_map
      __align_to_seq_coords_map.push_back (seqpos);
    }

    else {
      // record for __seq_to_align_coords_map
      __seq_to_align_coords_map.push_back (alignpos);
      // records for __align_to_seq_coords_map
      __align_to_seq_coords_map.push_back (seqpos);
      ++seqpos;
    }
    
  }
  assert (__seq_to_align_coords_map.size() == __seqlength);
  assert (__align_to_seq_coords_map.size() == __row_path->size());
  assert (seqpos == __seqlength); // sanity check

}

Alignment_row* Alignment_row::subalignment (const unsigned start, const unsigned end) const {

  assert (!__seqlength || (end < length()));

  Row_path* row_path = NULL;
  if (end < start)
    row_path = new Row_path();
  else
    row_path = new Row_path (__row_path->begin() + start, __row_path->begin() + end + 1);

  return new Alignment_row (row_path);
  
}

Alignment_row::Row_path* Alignment_row::subalignment_row_path (const unsigned start, const unsigned end) const {

  assert (!__seqlength || (end < length()));
  if (end < start)
    return new Row_path();
  return new Row_path (__row_path->begin() + start, __row_path->begin() + end + 1);
  
}

void Alignment_row::reverse() {

  // reverse the alignment
  std::reverse (__row_path->begin(), __row_path->end());

  // and then re-build coordinate indices (it is NOT sufficient to merely reverse these!)
  build_coordinate_indices();

}

Alignment::Alignment (Sequence_database& seq_db)
  : seq_db (&seq_db) {

}

Alignment::Alignment (const Alignment& parent) {

  // perform deep copy
  seq_db = parent.seq_db;
  for (std::vector<Alignment_row*>::const_iterator row = parent.__rows.begin(); row != parent.__rows.end(); ++row) {
    Alignment_row* a = new Alignment_row (**row);
    __rows.push_back (a);
  }
  __row_index = parent.__row_index;

}

Alignment& Alignment::Alignment::operator= (const Alignment& parent) {

  // check for self-assignment
  if (this == &parent)
    return *this;

  // perform deep copy
  seq_db = parent.seq_db;
  for (std::vector<Alignment_row*>::const_iterator row = parent.__rows.begin(); row != parent.__rows.end(); ++row) {
    Alignment_row* a = new Alignment_row (**row);
    __rows.push_back (a);
  }
  __row_index = parent.__row_index;

  return *this;

}

Alignment::~Alignment() {

  for (std::vector<Alignment_row*>::iterator row = __rows.begin(); row != __rows.end(); ++row)
    delete *row;

}

bool Alignment::detect_mfa (const std::string& filename) {

  return Sequence::detect_fasta (filename);

}

void Alignment::read_mfa (const std::string& filename, const bool strict /* = true */,
			  const bool verbose /* = true */) {

  clear();

  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << filename << "' for reading." << endl;
    exit (1);
  }

  if (verbose)
    cerr << "Reading multi-FASTA alignment '" << filename << "'...";

  Regexp re_name ("^[ \t]*" + Sequence::fasta_seq_start + "([^ \t]*)[ \t]*(.*)");

  std::string line;

  std::vector<std::string> gapped_seqs_names;
  std::vector<std::string> gapped_seqs;
  std::vector<std::string> gapped_seqs_info;

  bool in_seq_body = false;
  size_t curr_idx = 0; // dummy value to prevent compiler warnings
  while (!filestream.eof()) {

    getline (filestream, line);

    // are we at a name line?
    if (re_name.Match (line.c_str())) {
      if (re_name.SubStrings() < 1) { continue; }
      // store name
      const std::string name = re_name[1];
      curr_idx = gapped_seqs_names.size();
      gapped_seqs_names.push_back (name);
      // store info following name, if present
      std::string info;
      if (re_name.SubStrings() >= 2)
	info = re_name[2];
      gapped_seqs_info.push_back (info);
      // store empty row
      gapped_seqs.push_back ("");
      in_seq_body = true;
    }
    // if not, then we must be reading sequence data
    else if (in_seq_body) {
      // if this is the first line of data for this sequence
      if (gapped_seqs.size() < gapped_seqs_names.size())
	gapped_seqs.push_back (line);
      else
	gapped_seqs[curr_idx] += line;
    }

  }
  filestream.close();

  // assert sane
  if (gapped_seqs.size() != gapped_seqs_names.size()) {
    cerr << "ERROR: Sequence names and data don't match up.  Is your MFA file improperly formatted?" << endl;
    exit (1);
  }
  assert (gapped_seqs_names.size() == gapped_seqs_info.size());

  // now store
  for (size_t r = 0; r < gapped_seqs.size(); ++r) {

    const std::string& name = gapped_seqs_names[r];
    std::string& gapped_seq = gapped_seqs[r];
    std::string& info = gapped_seqs_info[r];

    // strip out whitespace
    std::string::iterator last_pos = std::remove_if (gapped_seq.begin(), gapped_seq.end(),
						     isspace);
    gapped_seq.erase (last_pos, gapped_seq.end());

    // remove newline from info, if present
    Util::chomp (info);

    // create ungapped sequence
    std::string ungapped_seq;
    ungapped_seq.reserve (gapped_seq.length());
    last_pos = std::remove_copy_if (gapped_seq.begin(), gapped_seq.end(),
				    ungapped_seq.begin(),
				    Alignment::is_gap_char);
    ungapped_seq.erase (last_pos, ungapped_seq.end()); // erase unnecessary positions

    // create alignment path
    Alignment_row::Row_path* row_path = new Alignment_row::Row_path (gapped_seq.length()); // allocate space
    std::transform (gapped_seq.begin(), gapped_seq.end(),
		    row_path->begin(),
		    std::not1 (std::ptr_fun (Alignment::is_gap_char)));

    // now store row
    Sequence* sequence = new Sequence (name, ungapped_seq, info);
    add_row (sequence,
	     row_path);
  }

  if (strict)
    assert_flush();

  if (verbose)
    cerr << "done." << endl;

}

const Alignment_row& Alignment::get_row (const std::string& name) const {

#ifndef NDEBUG
  if (!exists_row (name)) {
    cerr << "ERROR: No row named '" << name << "'." << endl;
    write_mfa (cerr, false);
  }
#endif
  assert (exists_row (name));

  return *__rows[__row_index.find (name)->second];

}

double Alignment::percent_id() const {

  return percent_id (0, columns() - 1);

}

double Alignment::percent_id (const unsigned start, const unsigned end) const {

  assert (start <= end);
  assert (end < columns());

  double id = 0.;   // sum of per-column percent IDs
  size_t cols = 0;  // number of columns used in calculation

  for (size_t c = start; c <= end; ++c) {

    // count the (maximal) fraction of identical characters in the column
    std::map<char, size_t> char_count;
    size_t colsize = 0; // # of non-gap characters in column
    for (size_t r = 0; r < rows(); ++r) {
      if (!is_gapped (r, c)) {
	++char_count[static_cast<char> (tolower (get_char (r, c)))];
	++colsize;
      }
    }

    // only count columns with > 1 non-gap character
    if (colsize < 2)
      continue;

    // find the most common character
    std::map<char, size_t>::const_iterator max = std::max_element (char_count.begin(), char_count.end(),
								   Util::Map_value_less<std::map<char, size_t> >());
    // if no character appears more than once,
    // then percent id = 0
    if (max->second > 1)
      id += static_cast<double> (max->second) / colsize;

    // increment columns counter
    ++cols;

  }

  if (cols == 0)
    return 0;

  return (id / cols);

}

double Alignment::gap_fraction() const {

  return gap_fraction (0, columns() - 1);

}

double Alignment::gap_fraction (const unsigned start, const unsigned end) const {

  assert (start <= end);
  assert (end < columns());

  if (columns() == 0 || rows() == 0)
    return 0;

  size_t gaps = 0;

  for (size_t c = start; c <= end; ++c) {
    for (size_t r = 0; r < rows(); ++r) {
      if (is_gapped (r, c))
	++gaps;
    }
  }

  return (static_cast<double> (gaps) / ((end - start + 1) * rows()));

}

const Sequence Alignment::get_gapped_row (const std::string& name) const {

  assert (exists_row (name));

  return get_gapped_row (__row_index.find (name)->second);

}

const Sequence Alignment::get_gapped_row (const size_t r) const {

  assert (r < rows());

  return get_row (r).get_gapped_seq (seq_db->get_seq (r), Alignment::gap_char);  

}

std::vector<std::string> Alignment::get_row_names() const {

  std::vector<std::string> names (rows());
  for (size_t r = 0; r < rows(); ++r)
    names[r] = get_row_name (r);

  return names;

}

void Alignment::add_row (Sequence* sequence,
			 Alignment_row::Row_path* row_path) {

#ifndef NDEBUG
  // check that sequence length matches (implied) alignment sequence length
  if (sequence->seq.length() != static_cast<size_t> (accumulate (row_path->begin(), row_path->end(), 0))) {
    cerr << "ERROR: Sequence length mismatch for " << sequence->name << "." << endl
	 << "Sequence:                " << sequence->seq << endl
	 << "Alignment_row::Row_path: " << Util::join (*row_path, "") << endl;
    exit (1);
  }
#endif
  assert (sequence->seq.length() == static_cast<size_t> (accumulate (row_path->begin(), row_path->end(), 0)));

  // store in Sequence_database
  if (!seq_db->exists_seq (sequence->name))
    seq_db->add_seq (sequence);
  // now store alignment information
  __rows.push_back (new Alignment_row (row_path));
  // only create a new entry in __row_index if there isn't one already
  if (!exists_row (sequence->name))
    __row_index.insert (make_pair (sequence->name, __rows.size() - 1));

#ifndef NDEBUG
  // check that the sequence indices for the rows in *this
  // and sequences in the Sequence_database& are still in sync
  // after this operation
  if (get_row_index (sequence->name) != seq_db->get_seq_index (sequence->name)) {
    cerr << "ERROR: Sequence indices are out of sync between the Alignment (index = " << get_row_index (sequence->name) << ") and corresponding Sequence_database& (index = " << seq_db->get_seq_index (sequence->name) << ") after storing '" << sequence->name << "'." << endl;
    exit (1);
  }
#endif

}

void Alignment::add_row (Sequence* sequence,
			 Alignment_row* alignment_row) {

#ifndef NDEBUG
  // check that sequence length matches alignment sequence length
  if (sequence->seq.length() != alignment_row->seqlength()) {
    cerr << "ERROR: Sequence length mismatch for " << sequence->name << "." << endl
	 << "Sequence:                " << sequence->seq << endl
	 << "Alignment_row::Row_path: " << *alignment_row << endl;
  }
#endif
  assert (sequence->seq.length() == alignment_row->seqlength());  

  // store in Sequence_database
  if (!seq_db->exists_seq (sequence->name))
    seq_db->add_seq (sequence);
  // now store alignment information
  __rows.push_back (alignment_row);
  // only create a new entry in __row_index if there isn't one already
  if (!exists_row (sequence->name))
    __row_index.insert (make_pair (sequence->name, __rows.size() - 1));

#ifndef NDEBUG
  // check that the sequence indices for the rows in *this
  // and sequences in the Sequence_database& are still in sync
  // after this operation
  if (get_row_index (sequence->name) != seq_db->get_seq_index (sequence->name)) {
    cerr << "ERROR: Sequence indices are out of sync between the Alignment (index = " << get_row_index (sequence->name) << ") and corresponding Sequence_database& (index = " << seq_db->get_seq_index (sequence->name) << ") after storing '" << sequence->name << "'." << endl;
    exit (1);
  }
#endif

}

void Alignment::set_row (const std::string& name,
			 Alignment_row::Row_path* row_path) {

  // confirm that such a sequence already exists in the Sequence_database
  if (!seq_db->exists_seq (name)) {
    cerr << "ERROR: No sequence 'name' found in the stored Sequence_database." << endl;
    exit (1);
  }

#ifndef NDEBUG
  const Sequence& sequence = seq_db->get_seq (name);
  // check that sequence length matches (implied) alignment sequence length
  if (sequence.seq.length() != static_cast<size_t> (accumulate (row_path->begin(), row_path->end(), 0))) {
    cerr << "ERROR: Sequence length mismatch for " << sequence.name << "." << endl
	 << "Sequence:                " << sequence.seq << endl
	 << "Alignment_row::Row_path: " << Util::join (*row_path, "") << endl;
    exit (1);
  }
  assert (sequence.seq.length() == static_cast<size_t> (accumulate (row_path->begin(), row_path->end(), 0)));
#endif

  // store alignment information

  // only create a new entry in __row_index if there isn't one already
  if (!exists_row (name)) {
    __rows.push_back (new Alignment_row (row_path));
    __row_index.insert (make_pair (name, __rows.size() - 1));
  } else {
    __rows[get_row_index (name)] = new Alignment_row (row_path);
  }

#ifndef NDEBUG
  // check that the sequence indices for the rows in *this
  // and sequences in the Sequence_database& are still in sync
  // after this operation
  if (get_row_index (name) != seq_db->get_seq_index (name)) {
    cerr << "ERROR: Sequence indices are out of sync between the Alignment (index = " << get_row_index (name) << ") and corresponding Sequence_database& (index = " << seq_db->get_seq_index (name) << ") after storing '" << name << "'." << endl;
    exit (1);
  }
#endif

}

void Alignment::clear() {

  seq_db->clear();

  for (std::vector<Alignment_row*>::iterator row = __rows.begin(); row != __rows.end(); ++row)
    delete *row;
  __rows.clear();

  __row_index.clear();

}

Sequence Alignment_row::get_gapped_seq (const Sequence& ungapped, const char gap_char) const {

#ifndef NDEBUG
  // check lengths ok
  if (ungapped.length() != static_cast<size_t> (accumulate (__row_path->begin(), __row_path->end(), 0))) {
    cerr << "ERROR: Sequence length mismatch for " << ungapped.name << "." << endl
	 << "Sequence:                " << ungapped.seq << endl
	 << "Alignment_row::Row_path: " << *this << endl;
  }
#endif
  assert (ungapped.length() == static_cast<size_t> (accumulate (__row_path->begin(), __row_path->end(), 0)));

  std::string gapped (__row_path->size(), gap_char);
  size_t pos = 0;                                        // sequence position
  for (size_t col = 0; col < gapped.length(); ++col) {   // column (alignment position)
    if (!is_gapped (col))
      gapped[col] = ungapped.seq[pos++];
  }

  return Sequence (ungapped.name, gapped, ungapped.info);

}

void Alignment::assert_flush() const {

#ifndef NDEBUG

  for (std::vector<Alignment_row*>::const_iterator row = __rows.begin(); row != __rows.end(); ++row) {
    if (columns() != (*row)->length()) {
      write_mfa (cerr, false);
      cerr << "ERROR: Alignment not flush." << endl;
      exit (1);
    }
  }

#endif

}

void Alignment::write_mfa (std::ostream& o, const bool strict /* = true */) const {

  if (strict)
    assert_flush();

  // align to left
  o.setf (std::ios_base::left);

  for (size_t r = 0; r < rows(); ++r) {

    const Sequence& sequence = get_gapped_row (r);
    o << Sequence::fasta_seq_start << sequence.name
      << (sequence.info.length() ? (' ' + sequence.info) : "") << endl
      << sequence.seq << endl;
  }

}

Stockholm::Stockholm (Sequence_database& seq_db)
  : Alignment (seq_db) {

}

void Stockholm::read_stockholm (const std::string& filename, const bool strict /* = true */,
				const bool verbose /* = true */) {

  clear();

  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << filename << "' for reading." << endl;
    exit (1);
  }

  if (verbose)
    cerr << "Reading Stockholm alignment '" << filename << "'...";

  Regexp re_stock ("^[ \t]*#[ \t]*" + format_identifier + "[ \t]*" + version_identifier + "[ \t]*$");               // format & version identifiers
  Regexp re_sep ("^[ \t]*" + alignment_separator + "[ \t]*$");                                                      // alignment separator lines, "//"
  Regexp re_gf ("^[ \t]*" + file_annotation + "[ \t]+([^ \t]+)[ \t]+(.*)$");                                        // #=GF lines
  Regexp re_gc ("^[ \t]*" + column_annotation + "[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]*$");                           // #=GC lines
  Regexp re_gs ("^[ \t]*" + sequence_annotation + "[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]+(.*)$");                     // #=GS lines
  Regexp re_gr ("^[ \t]*" + sequence_column_annotation + "[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]*$");   // #=GR lines
  Regexp re_row ("^[ \t]*([^ \t#]+)[ \t]+([^ \t]*)[ \t]*$");                                                        // alignment row data
  Regexp re_nonwhite ("[^ \t]");                                                                                    // whitespace

  // first store gapped sequences in database,
  // then parse them later to initialize the alignment

  // hold raw (gapped) sequence data read in from file
  std::vector<std::string> gapped_seqs;
  std::vector<std::string> gapped_seqs_names;

  std::string line;
  while (!filestream.eof()) {

    getline (filestream, line);
    Util::chomp (line);

    // ignore format-version identifiers
    if (re_stock.Match (line.c_str()))
      continue;

    // break if we encounter an alignment separator
    else if (re_sep.Match (line.c_str()))
      break;

    // #=GF
    else if (re_gf.Match (line.c_str()) && re_gf.SubStrings() == 2) {
      const std::string& gf_key = re_gf[1];
      const std::string& gf_data = re_gf[2];
      __gf_annot.push_back (std::make_pair (gf_key, gf_data));
      __gf_index[gf_key].insert (__gf_annot.size() - 1);
    }

    // #=GC
    else if (re_gc.Match (line.c_str()) && re_gc.SubStrings() == 2)
      __gc_annot[re_gc[1]].append (re_gc[2]);

    // #=GS
    else if (re_gs.Match (line.c_str()) && re_gs.SubStrings() == 3)
      __gs_annot[re_gs[1]][re_gs[2]] = re_gs[3];

    // #=GR
    else if (re_gr.Match (line.c_str()) && re_gr.SubStrings() == 3)
      __gr_annot[re_gr[1]][re_gr[2]].append (re_gr[3]);

    // row data
    else if (re_row.Match (line.c_str()) && re_row.SubStrings() == 2) {

      const std::string name = re_row[1];
      const std::string gapped_seq = re_row[2];

      // is it the first time that we've seen this sequence?
      if (!exists_row (name)) {

	gapped_seqs_names.push_back (name);
	__row_index.insert (std::make_pair (name, __row_index.size()));
	gapped_seqs.push_back (gapped_seq);

	if (__gs_annot.find (name) == __gs_annot.end())
	  __gs_annot[name] = Annotation();

	if (__gr_annot.find (name) == __gr_annot.end())
	  __gr_annot[name] = Annotation();

      }
      // if not, then just append sequence data
      else {
	gapped_seqs[__row_index.find (name)->second] += gapped_seq;
      }

    }
      
    else if (line.size() && re_nonwhite.Match (line.c_str())) {
      cerr << "WARNING: couldn't parse the following alignment line:" << endl << line << endl;
      continue;
    }

  }
  filestream.close();

  // now parse this gapped sequence data into an alignment
  for (size_t r = 0; r < gapped_seqs.size(); ++r) {

    const std::string& name = gapped_seqs_names[r];
    const std::string& gapped_seq = gapped_seqs[r];

    // create ungapped sequence
    std::string ungapped_seq;
    ungapped_seq.reserve (gapped_seq.length());
    std::string::iterator last_pos = remove_copy_if (gapped_seq.begin(), gapped_seq.end(),
						     ungapped_seq.begin(),
						     Alignment::is_gap_char);
    ungapped_seq.erase (last_pos, ungapped_seq.end()); // erase unnecessary positions

    // create alignment path
    Alignment_row::Row_path* row_path = new Alignment_row::Row_path (gapped_seq.length()); // allocate space
    transform (gapped_seq.begin(), gapped_seq.end(),
	       row_path->begin(),
	       std::not1 (std::ptr_fun (Alignment::is_gap_char)));

    // now store row
    Sequence* sequence = new Sequence (name, ungapped_seq);
    add_row (sequence,
	     row_path);

  }

  if (strict)
    assert_all_flush();

  if (verbose)
    cerr << "done." << endl;

}

void Stockholm::write_stockholm (std::ostream& o, const bool strict /* = true */) const {

  if (strict)
    assert_all_flush();

  o << alignment_header << endl;

  // align to left
  o.setf (std::ios_base::left);

  // get the field width for the key columns
  // max_name_width takes into account sequence names, GF/GS/GC/GR annotations,
  // and GS annotations for internal nodes...
  // scroll down to see where it gets its final value
  size_t max_name_width = 0;
  for (size_t r = 0; r < rows(); r++) {
    // sequence names
    const std::string& name = get_row_name (r);
    max_name_width = std::max (max_name_width, name.length());
    // 6 chars for "#=GS name "
    const Row_annotation::const_iterator gs_annot_row = __gs_annot.find (name);
    if (gs_annot_row != __gs_annot.end()) {
      for (Annotation::const_iterator gs = gs_annot_row->second.begin();gs != gs_annot_row->second.end(); ++gs)
	max_name_width = std::max (max_name_width, gs->first.size() + name.length() + 6);
    }
    // 6 chars for "#=GR name "
    const Row_annotation::const_iterator gr_annot_row = __gr_annot.find (name);
    if (gr_annot_row != __gr_annot.end()) {
      for (Annotation::const_iterator gr = gr_annot_row->second.begin();gr != gr_annot_row->second.end(); ++gr)
	max_name_width = std::max (max_name_width, gr->first.size() + name.length() + 6);
    }
  }
  // 5 chars for "#=GF "
  for (std::vector<std::pair<std::string, std::string> >::const_iterator gf = __gf_annot.begin(); gf != __gf_annot.end(); ++gf)
    max_name_width = std::max (max_name_width, gf->first.length() + 5);
  // 5 chars for "#=GC "
  for (Annotation::const_iterator gc = __gc_annot.begin();gc != __gc_annot.end(); ++gc)
    max_name_width = std::max (max_name_width, gc->first.size() + 5);

  // #=GF lines
  for (std::vector<std::pair<std::string, std::string> >::const_iterator gf = __gf_annot.begin(); gf != __gf_annot.end(); ++gf) {
    std::string key = file_annotation + ' ' + gf->first;
    o.width (max_name_width + 1);
    o << key << gf->second << endl;
  }

  // #=GS lines
  for (size_t r = 0; r < rows(); ++r) {
    const std::string& name = get_row_name (r);
    const Row_annotation::const_iterator gs_annot_row = __gs_annot.find (name);
    if (gs_annot_row != __gs_annot.end()) {
      for (Annotation::const_iterator gs = gs_annot_row->second.begin(); gs != gs_annot_row->second.end(); ++gs) {
	std::string key = sequence_annotation + ' ' + name + ' ' + gs->first;
	o.width (max_name_width + 1);
	o << key << gs->second << endl;
      }
    }

  }

  // main body of alignment
  for (size_t r = 0; r < rows(); r++) {
    const std::string& name = get_row_name (r);

    // sequence data
    o.width (max_name_width + 1);
    const Sequence& sequence = get_gapped_row (r);
    o << sequence.name
      << sequence.seq << endl;

    // #=GR lines
    const Row_annotation::const_iterator gr_annot_row = __gr_annot.find (name);
    if (gr_annot_row != __gr_annot.end()) {
      for (Annotation::const_iterator gr = gr_annot_row->second.begin(); gr != gr_annot_row->second.end(); ++gr) {
	std::string key = sequence_column_annotation + ' ' + name + ' ' + gr->first;
	o.width (max_name_width + 1);
	o << key << gr->second << endl;
      }
    }

  }

  // #=GC lines
  for (Annotation::const_iterator gc = __gc_annot.begin(); gc != __gc_annot.end(); ++gc) {
    std::string key;
    key = column_annotation + ' ' + gc->first;
    o.width (max_name_width + 1);
    o << key << gc->second << endl;
  }

  o << alignment_separator << endl;

}

void Stockholm::read_stockholm_or_mfa (const std::string& filename, const bool strict /* = true */,
				       const bool verbose /* = true */) {

  if (!detect_mfa (filename))
    read_stockholm (filename, strict,
		    verbose);
  else
    read_mfa (filename, strict,
	      verbose);

}

std::string Stockholm::get_gf_annot (const std::string& key) const {

  std::string s;
  const std::map<std::string, std::set<size_t> >::const_iterator index = __gf_index.find (key);
  if (index != __gf_index.end()) {
    for (std::set<size_t>::const_iterator i = index->second.begin(); i != index->second.end(); ++i)
      s.append (__gf_annot[*i].second);
  }

  return s;

}

std::string Stockholm::get_gc_annot (const std::string& key) const {

  std::string s;
  const Annotation::const_iterator gc_iter = __gc_annot.find (key);
  if (gc_iter != __gc_annot.end())
    s = gc_iter->second;

  return s;

}

std::string Stockholm::get_gs_annot (const std::string& seq, const std::string& key) const {

  std::string s;
  const Row_annotation::const_iterator gs_row_annot = __gs_annot.find (seq);
  if (gs_row_annot != __gs_annot.end()) {
    const Annotation::const_iterator gs = gs_row_annot->second.find (key);
    if (gs != gs_row_annot->second.end())
      s = gs->second;
  }

  return s;

}

std::string Stockholm::get_gr_annot (const std::string& seq, const std::string& key) const {

  std::string s;
  const Row_annotation::const_iterator gr_row_annot = __gr_annot.find (seq);
  if (gr_row_annot != __gr_annot.end()) {
    const Annotation::const_iterator gr = gr_row_annot->second.find (key);
    if (gr != gr_row_annot->second.end())
      s = gr->second;
  }

  return s;

}

void Stockholm::clear() {

  this->Alignment::clear();
  clear_annot();

}

void Stockholm::clear_annot() {

  __gf_annot.clear();
  __gc_annot.clear();
  __gs_annot.clear();
  __gr_annot.clear();
  __gf_index.clear();

}

void Stockholm::add_gf_annot (const std::string& key, const std::string& value) {
  __gf_annot.push_back (std::make_pair (key, value));
  __gf_index[key].insert (__gf_annot.size() - 1);
}

void Stockholm::set_gc_annot (const std::string& key, const std::string& value) {
  __gc_annot[key] = value;
}

void Stockholm::set_gs_annot (const std::string& seq, const std::string& key, const std::string& value) {
  __gs_annot[seq].insert (std::make_pair (key, value));
}


void Stockholm::set_gr_annot (const std::string& seq, const std::string& key, const std::string& value) {
  __gr_annot[seq].insert (std::make_pair (key, value));
}

Stockholm* Stockholm::subalignment (Sequence_database& seq_db_subalign,
				    const unsigned start, const unsigned end) const {

  Stockholm* subalign = new Stockholm (seq_db_subalign);
  subalign->clear();

  // #=GF lines
  for (std::vector<std::pair<std::string, std::string> >::const_iterator gf = __gf_annot.begin(); gf != __gf_annot.end(); ++gf)
    subalign->add_gf_annot (gf->first, gf->second);

  // #=GS lines
  for (size_t r = 0; r < rows(); ++r) {
    const std::string& name = get_row_name (r);
    const Row_annotation::const_iterator gs_annot_row = __gs_annot.find (name);
    if (gs_annot_row != __gs_annot.end()) {
      for (Annotation::const_iterator gs = gs_annot_row->second.begin(); gs != gs_annot_row->second.end(); ++gs)
	subalign->set_gs_annot (name, gs->first, gs->second);
    }
  }

  // main body of alignment
  for (size_t r = 0; r < rows(); ++r) {

    const std::string& name = get_row_name (r);
    const Alignment_row& row = get_row (r);

    // sequence data
    const Interval seq_coords = row.map_align_to_seq (start, end);
    subalign->add_row (seq_db->get_seq (r).subsequence (seq_coords.start, seq_coords.end),
		       row.subalignment_row_path (start, end));
    // (note that the called functions properly handle the case of the 0-length interval (1, 0)
    // which map_align_to_seq may return)

    // #=GR lines
    const Row_annotation::const_iterator gr_annot_row = __gr_annot.find (name);
    if (gr_annot_row != __gr_annot.end()) {
      for (Annotation::const_iterator gr = gr_annot_row->second.begin(); gr != gr_annot_row->second.end(); ++gr)
	subalign->set_gr_annot (name, gr->first, gr->second.substr (start, start < end ? end - start + 1 : 0)); // catch case of the empty alignment
    }

  }

  // #=GC lines
  for (Annotation::const_iterator gc = __gc_annot.begin(); gc != __gc_annot.end(); ++gc)
    subalign->set_gc_annot (gc->first, gc->second.substr (start, start < end ? end - start + 1 : 0));

  return subalign;

}

void Stockholm::assert_all_flush() const {

#ifndef NDEBUG

  assert_flush();

  // now check annotations
  for (size_t r = 0; r < rows(); ++r) {

    const std::string& name = get_row_name (r);

    // #=GR lines
    const Row_annotation::const_iterator gr_annot_row = __gr_annot.find (name);
    if (gr_annot_row != __gr_annot.end()) {
      for (Annotation::const_iterator gr = gr_annot_row->second.begin(); gr != gr_annot_row->second.end(); ++gr) {
	if (columns() != gr->second.length()) {
	  write_stockholm (cerr, false);
	  cerr << "ERROR: Alignment annotations not flush." << endl;
	  exit (1);
	}
      }
    }

  }

  // #=GC lines
  for (Annotation::const_iterator gc = __gc_annot.begin(); gc != __gc_annot.end(); ++gc) {
    if (columns() != gc->second.length()) {
      write_stockholm (cerr, false);
      cerr << "ERROR: Alignment annotations not flush." << endl;
      exit (1);
    }
  }

#endif

}

void Stockholm::revcomp (const Alphabet& alphabet) {

  // revcomp sequences
  seq_db->revcomp (alphabet);

  // now reverse per-column annotations
  for (size_t r = 0; r < rows(); ++r) {

    const std::string& name = get_row_name (r);

    // reverse alignment row
    // note that this method properly re-builds the coordinate indices with Alignment_row::build_coordinate_indices
    get_row (r).reverse();

    // #=GR lines
    const Row_annotation::iterator gr_annot_row = __gr_annot.find (name);
    if (gr_annot_row != __gr_annot.end()) {
      for (Annotation::iterator gr = gr_annot_row->second.begin(); gr != gr_annot_row->second.end(); ++gr)
	std::reverse (gr->second.begin(), gr->second.end());
    }

  }

  // #=GC lines
  for (Annotation::iterator gc = __gc_annot.begin(); gc != __gc_annot.end(); ++gc)
    std::reverse (gc->second.begin(), gc->second.end());

}

Stockholm Stockholm::get_codon_from_aa_alignment (Sequence_database& seq_db_codon) const {

  // sanity check on number of seqs
  assert (this->rows() == seq_db_codon.size());

  // calculate the appropriate number of columns for the codon alignment
  size_t cols = 3 * this->columns(); // columns in aa alignment
  for (size_t i = 0; i < seq_db_codon.size(); ++i) {
    const Sequence& sequence = seq_db_codon.get_seq (i);
    const unsigned rem = sequence.length() % 3;        // add columns for overhanging nt
    cols += rem;
  }

  // initialize codon alignment
  Stockholm stockholm_codon (seq_db_codon);

  // #=GF lines
  stockholm_codon.__gf_annot = this->__gf_annot;
  stockholm_codon.__gf_index = this->__gf_index;

  // #=GS lines
  for (size_t r = 0; r < rows(); ++r) {
    const std::string& name = get_row_name (r);
    const Row_annotation::const_iterator gs_annot_row = __gs_annot.find (name);
    if (gs_annot_row != __gs_annot.end()) {
      for (Annotation::const_iterator gs = gs_annot_row->second.begin(); gs != gs_annot_row->second.end(); ++gs)
	stockholm_codon.set_gs_annot (name, gs->first, gs->second);
    }
  }

  // store nt alignment
  unsigned indent = 0;
  for (size_t i = 0; i < seq_db_codon.size(); ++i) {

    const std::string& name = seq_db_codon.get_seq (i).name;

    // #=GR lines
    const Row_annotation::const_iterator gr_annot_row = __gr_annot.find (name);
    if (gr_annot_row != __gr_annot.end()) {
      for (Annotation::const_iterator gr = gr_annot_row->second.begin(); gr != gr_annot_row->second.end(); ++gr) {
	std::string gr_codon (cols, Stockholm::annotation_wildcard_char);
	for (size_t s = 0; s < gr->second.length(); ++s)
	  gr_codon[(3 * s)] = gr_codon[(3 * s) + 1] = gr_codon[(3 * s) + 2] = gr->second[s];
	stockholm_codon.set_gr_annot (name, gr->first, gr_codon);
      }
    }

    // store this row of alignment as vector of booleans
    const Alignment_row& row_aa = this->get_row (i);
    Alignment_row::Row_path* row_path_nt = new Alignment_row::Row_path (cols, false);
    size_t ii;
    for (ii = 0; ii < row_aa.length(); ++ii) {
      // if a character (aa) here, mark the corresponding codon as aligned
      if (!row_aa.is_gapped (ii))
	(*row_path_nt)[(3 * ii)] = (*row_path_nt)[(3 * ii) + 1] = (*row_path_nt)[(3 * ii) + 2] = true;
    }

    // store overhangs if necessary as unaligned nt
    const unsigned len_aa = this->get_row (i).seqlength();
    const unsigned len_nt = seq_db_codon.get_seq (i).length();
    if (static_cast<unsigned> (len_nt / 3) != len_aa) { // sanity check on sequence lengths
      cerr << "ERROR: Sequence lengths in protein and codon space don't agree (even when rounded to the lowest multiple of 3)." << endl;
      exit (1);
    }

    if (len_nt > 3 * len_aa) {
      for (size_t overhang = len_nt - (3 * len_aa); overhang > 0; --overhang)
	(*row_path_nt)[(3 * ii) + overhang + indent - 1] = true;
      indent += len_nt - (3 * len_aa);
    }

    // add the row to the alignment
    stockholm_codon.set_row (seq_db_codon.get_seq (i).name, row_path_nt);

  }

  // #=GC lines
  for (Annotation::const_iterator gc = __gc_annot.begin(); gc != __gc_annot.end(); ++gc) {
    std::string gc_codon (cols, Stockholm::annotation_wildcard_char);
    for (size_t s = 0; s < gc->second.length(); ++s)
      gc_codon[(3 * s)] = gc_codon[(3 * s) + 1] = gc_codon[(3 * s) + 2] = gc->second[s];
    stockholm_codon.set_gc_annot (gc->first, gc_codon);
  }

  stockholm_codon.assert_all_flush();

  return stockholm_codon;

}

void Stockholm::annotate_with_statistics() {

  add_gf_annot (Stockholm::percent_id_annotation, Util::to_string (percent_id()));
  add_gf_annot (Stockholm::gap_fraction_annotation, Util::to_string (gap_fraction()));

}

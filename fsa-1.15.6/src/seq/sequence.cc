
/**
 * \file sequence.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <algorithm>

#include "util/regexp.h"
#include "math/mathematics.h"
#include "seq/sequence.h"

using namespace fsa;

const std::string Sequence::fasta_seq_start = ">";

// Mapping from codons to amino acids.
// This really should be a static member variable,
// but I kept getting annoying linker errors, so here it is.
// Nucleotide->integer mapping is e.g. AGC => [0][3][2].
// Stop codons are 'X'.
const char Translated_sequence::__codon_map[4][4][4] = {
  { { 'K', 'N', 'K', 'N' },
    { 'T', 'T', 'T', 'T' },
    { 'R', 'S', 'R', 'S' },
    { 'I', 'I', 'M', 'I' } },
  { { 'Q', 'H', 'Q', 'H' },
    { 'P', 'P', 'P', 'P' },
    { 'R', 'R', 'R', 'R' },
    { 'L', 'L', 'L', 'L' } },
  { { 'E', 'D', 'E', 'D' },
    { 'A', 'A', 'A', 'A' },
    { 'G', 'G', 'G', 'G' },
    { 'V', 'V', 'V', 'V' } },
  { { 'X', 'Y', 'X', 'Y' },
    { 'S', 'S', 'S', 'S' },
    { 'X', 'C', 'W', 'C' },
    { 'L', 'F', 'L', 'F' } },
};



Sequence::Sequence (const std::string& name, const std::string& seq)
  : name (name), seq (seq)
{ }

Sequence::Sequence (const std::string& name, const std::string& seq,
		    const std::string& info)
  : name (name), seq (seq), info (info)
{ }

bool Sequence::detect_fasta (const std::string& filename) {

  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << filename << "' for reading." << endl;
    exit (1);
  }

  Regexp re_name ("^[ \t]*" + fasta_seq_start + "([^ \t]*)[ \t]*(.*)");

  std::string line;
  while (!filestream.eof()) {

    getline (filestream, line);
    if (re_name.Match (line.c_str())) {
      filestream.close();
      return true;
    }

  }
  filestream.close();

  return false;  

}

Sequence* Sequence::subsequence (const unsigned start, const unsigned end) const {

  assert (!length() || (end < length()));
  Sequence* sequence = NULL;
  if (end < start)
    sequence = new Sequence (name, "");
  else
    sequence = new Sequence (name, seq.substr (start, end - start + 1), info);

  return sequence;

}

void Sequence::revcomp (const Alphabet& alphabet, bool (*is_gap_char) (char) /* = NULL */) {

  alphabet.revcomp (seq, is_gap_char);

}

Sequence Sequence::revcomp (const Alphabet& alphabet, bool (*is_gap_char) (char) /* = NULL */) const {

  Sequence sequence = *this;
  sequence.revcomp (alphabet, is_gap_char);
  return sequence;

}

Translated_sequence::Translated_sequence (const Sequence& orig_seq)
  : orig_seq (orig_seq),
    __alphabet (DNA_alphabet()) {             // impose DNA alphabet

  __forward.resize (3);
  __reverse.resize (3);

  // size strings for all reading frames for speed
  const unsigned rem = orig_seq.length() % 3;
  const unsigned len = (orig_seq.length() - rem) / 3;
  for (size_t frame = 0; frame < 3; ++frame) {
    // calculate length for this frame
    const size_t l = len - ((rem - frame) >= 0 ? 0 : 1); // remember possibility of overhang at end
    __forward[frame].reserve (l);
    __reverse[frame].reserve (l);
  }
  
  // forward strand
  unsigned zero;
  unsigned one = __alphabet.get_char_index (orig_seq.seq[0]);  // initialize to "previous" values from fake last frame
  unsigned two = __alphabet.get_char_index (orig_seq.seq[1]);
  for (size_t i = 2; i < orig_seq.length(); ++i) {

    // get amino acid
    zero = one;                                            // increment frame
    one = two;
    two = __alphabet.get_char_index (orig_seq.seq[i]);     // and get new character
    const char ch = __codon_map[zero][one][two];

    // get frame and position in translated sequence
    const unsigned frame = (i - 2) % 3; // because we begin at frame 0 but i = 2

    // store character
    __forward[frame] += ch;

  }

  // reverse strand
  std::string orig_seq_revcomp = __alphabet.revcomp (orig_seq.seq);

  one = __alphabet.get_char_index (orig_seq_revcomp[0]);  // initialize to "previous" values from fake last frame
  two = __alphabet.get_char_index (orig_seq_revcomp[1]);
  for (size_t i = 2; i < orig_seq.length(); ++i) {

    // get amino acid
    zero = one;                                                // increment frame
    one = two;
    two = __alphabet.get_char_index (orig_seq_revcomp[i]);     // and get new character
    const char ch = __codon_map[zero][one][two];

    // get frame and position in translated sequence
    const unsigned frame = (i - 2) % 3; // because we begin at frame 0 but i = 2

    // store character
    __reverse[frame] += ch;

  }

}

void Sequence::init_hardmasking (size_t min_hardmask_length, bool (*is_hardmask_char) (char)) {

  // mark sequence as hardmasked
  __is_hardmasked = true;

  // find hardmasked intervals
  bool in_run = false;
  size_t start = 0;
  size_t len = 0;
  for (size_t i = 0; i < seq.length(); ++i) {
    if (is_hardmask_char (seq[i])) {
      // if just started a run of hardmasked characters
      if (!in_run) {
	start = i;
	len = 0;
	in_run = true;
      }
      // increment length of run
      ++len;
    }
    else {
      // if just ended a run
      if (in_run) {
	// store if meets length criteria
	if (len >= min_hardmask_length)
	  __masked_intervals.push_back (Interval (start, i - 1)); // -1 b/c now outside of run
	in_run = false;
      }
    }
  }

  // check for case of ending in a run
  if (in_run) {
    // store if meets length criteria
    if (len >= min_hardmask_length)
      __masked_intervals.push_back (Interval (start, seq.length()));
  }  

}

Sequence Sequence::get_stripped_sequence() const {

  if (!__is_hardmasked) {
    cerr << "ERROR: Sequence isn't hardmasked." << endl;
    exit (1);
  }

  // handle case of no masked intervals
  if (!__masked_intervals.size())
    return *this;

  std::string stripped;

  // iterate over masked intervals, avoiding them as we go  
  size_t start = 0;
  for (std::vector<Interval>::const_iterator iter = __masked_intervals.begin(); iter != __masked_intervals.end(); ++iter) {
    for (size_t i = start; i < iter->start; ++i)
      stripped += seq[i];
    start = iter->end + 1;
  }
  // catch the sequence past the last interval
  for (size_t i = __masked_intervals.back().end + 1; i < seq.length(); ++i)
    stripped += seq[i];

  return Sequence (name, stripped, info);

}

bool Sequence::is_pos_hardmasked (const unsigned orig_pos) const {

  if (!__is_hardmasked) {
    cerr << "ERROR: Sequence isn't hardmasked." << endl;
    exit (1);
  }

  for (std::vector<Interval>::const_iterator iter = __masked_intervals.begin(); iter != __masked_intervals.end(); ++iter) {
    if ((orig_pos >= iter->start) && (orig_pos <= iter->end))
      return true;
  }

  return false;
}

unsigned Sequence::map_stripped_to_orig (const unsigned pos) const {

  if (!__is_hardmasked) {
    cerr << "ERROR: Sequence isn't hardmasked." << endl;
    exit (1);
  }

  size_t orig_pos = pos;
  for (std::vector<Interval>::const_iterator iter = __masked_intervals.begin(); iter != __masked_intervals.end(); ++iter) {
    
    // have we found it?
    // (is it before this interval begins?)
    if (orig_pos < iter->start)
      return orig_pos;
    
    // if position is past this masked interval,
    // increment offset appropriately
    if (orig_pos >= iter->start)
      orig_pos += iter->end - iter->start + 1;
    
  }

  // check sane
  if (orig_pos >= seq.length()) {
    cerr << "map_stripped_to_orig: mapped " << pos << " -> " << orig_pos << "; seqlen = " << seq.length() << endl;
  }
  assert (orig_pos < seq.length());

  return orig_pos;

}

unsigned Sequence::map_orig_to_stripped (const unsigned orig_pos) const {

  if (!__is_hardmasked) {
    cerr << "ERROR: Sequence isn't hardmasked." << endl;
    exit (1);
  }

  // sanity check  
  assert (orig_pos < (*this).length());

  // loop over the masked intervals,
  // subtracting their length as we go
  size_t pos = orig_pos;
  for (std::vector<Interval>::const_iterator iter = __masked_intervals.begin(); iter != __masked_intervals.end(); ++iter) {

    // return out-of-bounds answer if original coordinate
    // falls within a masked interval
    if ((orig_pos >= iter->start) && (orig_pos <= iter->end))
      return seq.length();

    // have we reached the original coordinate?
    if (orig_pos < iter->start)
      return pos;

    pos -= iter->length();
    
  }

  return pos;

}

void Sequence::write_fasta (std::ostream& o) const {

  o << fasta_seq_start << name
    << (info.length() ? (' ' + info) : "") << endl
    << seq << endl;

}

Sequence_database::Sequence_database (const Sequence_database& parent) {

  // deep copy
  for (std::vector<Sequence*>::const_iterator sequence = parent.__sequences.begin(); sequence != parent.__sequences.end(); ++sequence) {
    Sequence* s = new Sequence (**sequence);
    __sequences.push_back (s);
  }

  __sequence_index = parent.__sequence_index;

}

Sequence_database& Sequence_database::operator= (const Sequence_database& parent) {

  // check for self-assignment
  if (this == &parent)
    return *this;

  // clear old memory
  this->clear();

  // perform deep copy
  for (std::vector<Sequence*>::const_iterator sequence = parent.__sequences.begin(); sequence != parent.__sequences.end(); ++sequence) {
    Sequence* s = new Sequence (**sequence);
    __sequences.push_back (s);
  }

  __sequence_index = parent.__sequence_index;

  return *this;

}

Sequence_database::~Sequence_database() {

  for (std::vector<Sequence*>::iterator sequence = __sequences.begin(); sequence != __sequences.end(); ++sequence)
    delete *sequence;

}

void Sequence_database::read_fasta (const std::string& filename, bool (*is_gap_char) (char) /* = NULL */,
				    const bool strip_leading_chr /* = true */,
				    const bool verbose /* = true */) {

  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open()) {
    cerr << "ERROR: Couldn't open file '" << filename << "' for reading." << endl;
    exit (1);
  }

  if (verbose)
    cerr << "Reading FASTA-format sequence file '" << filename << "'...";

  Regexp re_name ("^[ \t]*" + Sequence::fasta_seq_start + "([^ \t]*)[ \t]*(.*)");

  std::string line;

  std::vector<std::string> seqs_names;
  std::vector<std::string> seqs;
  std::vector<std::string> seqs_info;

  bool in_seq_body = false;
  size_t curr_idx = 0; // dummy value to prevent compiler warnings
  while (!filestream.eof()) {

    getline (filestream, line);

    // are we at a name line?
    if (re_name.Match (line.c_str())) {
      if (re_name.SubStrings() < 1) { continue; }
      // store name
      const std::string name = re_name[1];
      curr_idx = seqs_names.size();
      seqs_names.push_back (name);
      // store info following name, if present
      std::string info;
      if (re_name.SubStrings() >= 2)
	info = re_name[2];
      seqs_info.push_back (info);
      // store empty sequence
      seqs.push_back ("");
      in_seq_body = true;
    }
    // if not, then we must be reading sequence data
    else if (in_seq_body) {
      // if this is the first line of data for this sequence
      if (seqs.size() < seqs_names.size())
	seqs.push_back (line);
      else
	seqs[curr_idx] += line;
    }

  }
  filestream.close();

  // assert sane
  if (seqs.size() != seqs_names.size()) {
    cerr << "ERROR: Sequence names and data don't match up.  Is your FASTA file improperly formatted?" << endl;
    exit (1);
  }
  assert (seqs_names.size() == seqs_info.size());

  // now store
  for (size_t i = 0; i < seqs_names.size(); ++i) {

    std::string& name = seqs_names[i];
    std::string& seq = seqs[i];
    std::string& info = seqs_info[i];

    // remove leading chr if requested and present
    if (strip_leading_chr)
      Util::strip_leading_chr (name);

    // strip out whitespace
    std::string::iterator last_pos = std::remove_if (seq.begin(), seq.end(),
						     isspace);
    seq.erase (last_pos, seq.end());

    // strip out gaps if requested
    if (is_gap_char != 0) {
      last_pos = std::remove_if (seq.begin(), seq.end(),
				 is_gap_char);
      seq.erase (last_pos, seq.end());
    }

    // remove newline from info, if present
    Util::chomp (info);

    // store sequence
    Sequence* sequence = new Sequence (name, seq, info);
    add_seq (sequence);
  }

  if (verbose)
    cerr << "done." << endl;

}

void Sequence_database::write_fasta (std::ostream& o) const {

  for (std::vector<Sequence*>::const_iterator seq = this->begin(); seq != this->end(); ++seq)
    (*seq)->write_fasta (o);

}

void Sequence_database::add_seq (Sequence* sequence) {

  if (exists_seq (sequence->name)) {
    cerr << "ERROR: Duplicate sequence name '" << sequence->name << "'; this is forbidden." << endl
	 << "Note that descriptions following the sequence names are not considered part of the names themselves." << endl;
    exit (1);
  }
  __sequences.push_back (sequence);
  __sequence_index.insert (std::make_pair (sequence->name, __sequences.size() - 1));

}

void Sequence_database::add_seq (const Sequence& sequence) {

  if (exists_seq (sequence.name)) {
    cerr << "ERROR: Duplicate sequence name '" << sequence.name << "'; this is forbidden." << endl
	 << "Note that descriptions following the sequence names are not considered part of the names themselves." << endl;
    exit (1);
  }

  Sequence* sequence_p = new Sequence (sequence);
  __sequences.push_back (sequence_p);
  __sequence_index.insert (std::make_pair (sequence.name, __sequences.size() - 1));

}

void Sequence_database::revcomp (const Alphabet& alphabet, bool (*is_gap_char) (char) /* = NULL */) {

  for (std::vector<Sequence*>::iterator seq = this->begin(); seq != this->end(); ++seq)
    (*seq)->revcomp (alphabet, is_gap_char);

}

bool Sequence_database::matches_alphabet (const Alphabet& alphabet, const double threshold /* = 0.95 */) const {

  size_t matches = 0;
  size_t total = 0;

  for (size_t i = 0; i < size(); ++i) {
    const Sequence& sequence = get_seq (i);
    total += sequence.length();

    for (std::string::const_iterator s = sequence.seq.begin(); s != sequence.seq.end(); ++s) {
      if (alphabet.is_nondegen_char (*s) || alphabet.is_unknown_char (*s))
	++matches;
    }

  }

  // catch case of no sequence
  if (total == 0)
    return true;

  const double frac =  static_cast<double> (matches) / total;

  if (frac >= threshold) {
    cerr << "Sequences match a " << alphabet.name() << " alphabet." << endl;
    return true;
  }

  return false;
}

void Sequence_database::clear() {

  for (std::vector<Sequence*>::iterator sequence = __sequences.begin(); sequence != __sequences.end(); ++sequence)
    delete *sequence;
  __sequences.clear();

  __sequence_index.clear();

}

size_t Sequence_database::meanlength() const {

  size_t total = 0;
  for (size_t i = 0; i < size(); ++i) {
    const Sequence& sequence = get_seq (i);
    total += sequence.length();
  }
  
  return static_cast<size_t> (total / size());

}

size_t Sequence_database::median_length() const {

  std::vector<size_t> lengths (size());
  for (size_t i = 0; i < size(); ++i)
    lengths[i] = get_seq (i).length();
  
  return Mathematics::median (lengths);

}

std::vector<std::string> Sequence_database::get_sequence_list() const {

  std::vector<std::string> names;
  for (std::vector<Sequence*>::const_iterator sequence = begin(); sequence != end(); ++sequence)
    names.push_back ((*sequence)->name);

  return names;

}

Sequence_database Sequence_database::translate() const {

  Sequence_database tr_db;

  // store first reading frame of forward strand
  for (size_t i = 0; i < this->size(); ++i) {
    const Sequence& sequence = get_seq (i);
    Translated_sequence tr (sequence);
    Sequence* forward = new Sequence (tr.get_forward (0));
    tr_db.add_seq (forward);
  }

  return tr_db;

}

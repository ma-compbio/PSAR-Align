
/**
 * \file sequence.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <climits>
#include <algorithm>

#include "seq/alphabet.h"

using namespace fsa;

const std::string DNA_alphabet::DNA_alphabet_name = "DNA";
const std::string RNA_alphabet::RNA_alphabet_name = "RNA";
const std::string Protein_alphabet::Protein_alphabet_name = "Protein";



Alphabet::Alphabet (std::string name, size_t size, bool case_sensitive /* = false */)
  : __name (name), __size (size), __case_sensitive (case_sensitive) {

}

void Alphabet::init_chars (const std::string& chars, const std::string complement /* = "" */) {

  __has_complement = (complement.length() > 0);

  // check sane
  assert (__size == chars.length());
  if (__has_complement) {
    assert (__size == complement.length());
    __char_complement_list.resize (__size);
  }

  // now create the __char_index data structure:
  // this holds a mapping from (non-degenerate) characters in the alphabet
  // to their numerical indices in __char_list
  // for speed, we use a vector rather than a map;
  // elements are accessed by casting the character to int
  // if a character isn't in the alphabet then it's given an index of __size;
  // if it is, then it's given the corresponding index in __char_list
  // (which is smaller than __size by definition)

  // assign dummy value to all entries in __char_index
  // (indicating not present in alphabet) of __size
  // NB: CHAR_MAX is defined in climits
  __char_index.assign (CHAR_MAX + 1, __size);

  // store the alphabet
  __char_list.resize (__size);  
  for (size_t i = 0; i < chars.length(); ++i) {
    char ch = chars[i];
    if (!__case_sensitive)
      ch = tolower (ch);
    __char_list[i] = ch;
    if (__has_complement)
      __char_complement_list[i] = complement[i];
    assert (__char_index.size() > static_cast<size_t> (ch));
    __char_index[ch] = i;
  }

}

void Alphabet::add_degen_char (char ch, const std::string& nondegen) {

  // initialize __degen_char_map by storing a
  // dummy empty vector for all characters
  // (similar in spirit to how __char_index works; see init_chars)
  // if this is the first time that this function was called
  if (!__degen_char_map.size())
  __degen_char_map.assign (CHAR_MAX + 1, std::vector<char>());

  if (!__case_sensitive)
    ch = tolower (ch);

  std::vector<char> nondegen_chars (nondegen.length());
  for (size_t i = 0; i < nondegen.length(); ++i)
    nondegen_chars[i] = nondegen[i];
  __degen_char_map[ch] = nondegen_chars;
  
}

void Alphabet::set_unknown_char (const char ch) {

  __unknown_char = ch;

}

void Alphabet::make_nondegen (std::string& seq) const {

  for (std::string::iterator c = seq.begin(); c != seq.end(); ++c)
    *c = get_nondegen_char (*c);

}

std::string Alphabet::get_nondegen (const std::string& seq) const {

  std::string seq_nondegen;
  seq_nondegen.resize (seq.length());
  for (size_t i = 0; i < seq.length(); ++i)
    seq_nondegen[i] = get_nondegen_char (seq[i]);

  return seq_nondegen;

}

void Alphabet::revcomp (std::string& seq, bool (*is_gap_char) (char) /* = NULL */) const {

  if (!__has_complement) {
    cerr << "ERROR: Tried to reverse-complement under an alphabet without a complementary alphabet." << endl;
    exit (1);
  }

  // first reverse
  std::reverse (seq.begin(), seq.end());

  // then complement
  for (std::string::iterator ch = seq.begin(); ch != seq.end(); ++ch) {
    // if requested, ignore gap characters
    if ((is_gap_char != 0) && is_gap_char (*ch))
      continue;
    // complement nondegen chars
    if (is_nondegen_char (*ch)) {
      if (isupper (*ch))
	*ch = toupper (__char_complement_list[__char_index[tolower (*ch)]]);
      else
	*ch = tolower (__char_complement_list[__char_index[*ch]]);
    }
    // and set degen & unknown chars to __unknown_char
    else {
      *ch = isupper (*ch) ? toupper (__unknown_char) : __unknown_char;
    }
  }

}

std::string Alphabet::revcomp (const std::string& seq, bool (*is_gap_char) (char) /* = NULL */) const {
  
  std::string seq_revcomp = seq;
  revcomp (seq_revcomp, is_gap_char);
  return seq_revcomp;

}

DNA_alphabet::DNA_alphabet()
  : Alphabet (DNA_alphabet_name, 4, false) {

  init_chars ("acgt", "tgca");
  set_unknown_char ('n');
  add_degen_char ('u', "t");
  add_degen_char ('r', "ag");
  add_degen_char ('y', "ct");
  add_degen_char ('m', "ac");
  add_degen_char ('k', "gt");
  add_degen_char ('s', "cg");
  add_degen_char ('w', "at");
  add_degen_char ('h', "act");
  add_degen_char ('b', "cgt");
  add_degen_char ('v', "acg");
  add_degen_char ('d', "agt");

}

RNA_alphabet::RNA_alphabet()
  : Alphabet (RNA_alphabet_name, 4, false) {

  init_chars ("acgu", "ugca");
  set_unknown_char ('n');
  add_degen_char ('t', "u");
  add_degen_char ('r', "ag");
  add_degen_char ('y', "ct");
  add_degen_char ('m', "ac");
  add_degen_char ('k', "gt");
  add_degen_char ('s', "cg");
  add_degen_char ('w', "at");
  add_degen_char ('h', "act");
  add_degen_char ('b', "cgt");
  add_degen_char ('v', "acg");
  add_degen_char ('d', "agt");

}

Protein_alphabet::Protein_alphabet()
  : Alphabet (Protein_alphabet_name, 20, false) {

  init_chars ("acdefghiklmnpqrstvwy");
  set_unknown_char ('x');
  add_degen_char ('b', "nd");
  add_degen_char ('z', "qe");

}

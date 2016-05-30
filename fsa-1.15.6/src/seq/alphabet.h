
/**
 * \file alphabet.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 * The alphabet representation code is loosely based on Ian Holmes's
 * Alphabet class.
 */

#ifndef SEQ_ALPHABET_INCLUDED
#define SEQ_ALPHABET_INCLUDED

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "util/misc.h"

namespace fsa {

  /**
   * \brief Represent a sequence alphabet.
   */
  struct Alphabet {

  public:

    /**
     * \brief Constructor.
     *
     * Unless the alphabet is designated as case-sensitive,
     * alphabet information is stored as lower-case internally.
     * \param name alphabet name
     * \param size alphabet size
     * \param case_sensitive is the alphabet case-sensitive?
     */
    Alphabet (std::string name, size_t size, bool case_sensitive = false);

    /**
     * \brief Get the name of the alphabet.
     */
    const std::string& name() const { return __name; }

    /**
     * \brief Make a sequence nondegenerate.
     */
    void make_nondegen (std::string& sequence) const;

    /**
     * \brief Get nondegenerate sequence.
     */
    std::string get_nondegen (const std::string& sequence) const;

    /**
     * \brief Does the alphabet have a complement defined?
     */
    bool has_complement() const { return __has_complement; }

    /**
     * \brief Reverse-complement sequence under alphabet.
     *
     * Degenerate characters are set to __unknown_char after reverse-complementing.
     * If requested, ignores gap characters as defined by the passed function.
     * \param seq sequence to be reverse-complemented
     * \param is_gap_char function pointer defining a gap character
     */
    void revcomp (std::string& seq, bool (*is_gap_char) (char) = NULL) const;

    /**
     * \brief Reverse-complement sequence under alphabet.
     * \see revcomp
     */
    std::string revcomp (const std::string& seq, bool (*is_gap_char) (char) = NULL) const;

    /**
     * \brief Get alphabet size.
     */
    inline size_t size() const { return __size; };

    /**
     * \brief Is this a non-degenerate character in the alphabet?
     */
    inline bool is_nondegen_char (char ch) const;

    /**
     * \brief Is this a degenerate character in the alphabet?
     */
    inline bool is_degen_char (char ch) const;

    /**
     * \brief Does the alphabet contain a possibly-degenerate character?
     * \see is_nondegen_char
     * \see is_degen_char
     */
    inline bool contains_char (const char ch) const;

    /**
     * \brief Is this an unknown character?
     */
    inline bool is_unknown_char (const char ch) const;

    /**
     * \brief Convert a possibly-degenerate character to a non-degenerate character.
     *
     * Randomizes unknown characters across the entire alphabet.
     * Preserves character case.
     */
    inline char get_nondegen_char (const char ch) const;

    /**
     * \brief Get the numerical index for a character.
     *
     * If the character is degenerate or unknown,
     * randomizes the character prior to getting the index.
     */
    inline size_t get_char_index (char ch) const;

    /**
     * \brief Get a character from its numerical index.
     */
    inline char get_char_from_index (const size_t index) const;

  protected:

    /**
     * \brief Initialize characters in alphabet.
     */
    void init_chars (const std::string& chars, const std::string complement = "");

    /**
     * \brief Add a degenerate character to the alphabet.
     */
    void add_degen_char (char ch, const std::string& nondegen);

    /**
     * \brief Set the unknown character for the alphabet.
     */
    void set_unknown_char (const char ch);

  private:

    std::string __name;          ///< alphabet name
    size_t __size;               ///< alphabet size (# of non-degenerate characters)
    bool __case_sensitive;       ///< is the alphabet case-sensitive?
    bool __has_complement;       ///< does the alphabet have a complement?

    std::vector<char> __char_list;                          ///< ordered list of characters
    std::vector<char> __char_complement_list;               ///< ordered list of complementary characters
    std::vector<size_t> __char_index;                       ///< map from characters (cast to int) to index in __char_list
    char __unknown_char;                                    ///< character representing unknown character
    std::vector<std::vector<char> > __degen_char_map;       ///< map from degenerate to non-degenerate characters (similar in spirit to __char_index)

  };

  /**
   * \brief Represent a DNA alphabet.
   */
  struct DNA_alphabet : public Alphabet {

  public:

    /**
     * \brief Constructor.
     */
    DNA_alphabet();

    /**
     * \brief Define the hardmasked character (ignores case!).
     */
    static bool is_hardmask_char (const char ch) {
      return toupper (ch) == 'N';
    }

    /**
     * \brief Call all lower-case characters softmasked.
     */
    static bool is_softmasked (const char ch) {
      return islower (ch);
    }

  private:

    static const std::string DNA_alphabet_name;   ///< alphabet name

  };

  /**
   * \brief Represent a RNA alphabet.
   */
  struct RNA_alphabet : public Alphabet {

  public:

    /**
     * \brief Constructor.
     */
    RNA_alphabet();

  private:

    static const std::string RNA_alphabet_name;   ///< alphabet name

  };

  /**
   * \brief Represent a protein alphabet.
   *
   * Characters are in IUPAC alphabetical order, 'acdefghiklmnpqrstvwy'.
   */
  struct Protein_alphabet : public Alphabet {

  public:

    /**
     * \brief Constructor.
     */
    Protein_alphabet();

  private:

    static const std::string Protein_alphabet_name;   ///< alphabet name

  };

  inline bool Alphabet::is_nondegen_char (char ch) const {
    if (!__case_sensitive)
      ch = tolower (ch);
    // test for membership in the alphabet
    if (__char_index[ch] < __size)
      return true;
    return false;
  }

  inline bool Alphabet::is_degen_char (char ch) const {
    if (!__case_sensitive)
      ch = tolower (ch);
    if (__degen_char_map[ch].size())
      return true;
    return false;
  }

  inline bool Alphabet::contains_char (const char ch) const {
    return (is_nondegen_char (ch) || is_degen_char (ch));
  }

  inline bool Alphabet::is_unknown_char (const char ch) const {
    return ((__case_sensitive ? ch : tolower (ch)) == __unknown_char);
  }

  inline char Alphabet::get_nondegen_char (const char ch) const {

    // if upper-case
    if (isupper (ch)) {

      // if nondegenerate
      if (is_nondegen_char (ch))
	return ch;

      // else if degenerate
      else if (is_degen_char (ch)) {
	const std::vector<char>& nondegen = __degen_char_map[__case_sensitive ? ch : tolower (ch)];
	return toupper (nondegen[Util::rand (nondegen.size() - 1)]);
      }
      // else we don't know anything about it, so treat as unknown
      else
	return toupper (__char_list[Util::rand (__size - 1)]);

    }

    // if lower-case
    else {

      // if nondegenerate
      if (is_nondegen_char (ch))
	return ch;

      // else if degenerate
      else if (is_degen_char (ch)) {
	const std::vector<char>& nondegen = __degen_char_map[__case_sensitive ? ch : tolower (ch)];
	return tolower (nondegen[Util::rand (nondegen.size() - 1)]);
      }
      // else we don't know anything about it, so treat as unknown
      else {
	return tolower (__char_list[Util::rand (__size - 1)]);
      }

    }

  }

  inline size_t Alphabet::get_char_index (char ch) const {

    if (!__case_sensitive)
      ch = tolower (ch);

    // if this is a non-degenerate character, then return the index
    if (__char_index[ch] < __size)
      return __char_index[ch];

    // else randomize and return the corresponding index
    else {
      ch = get_nondegen_char (ch);
      return __char_index[ch];
    }

  }
  
  inline char Alphabet::get_char_from_index (const size_t index) const {
    
    assert (index < __char_list.size());
    return __char_list[index];

  }

}

#endif /* SEQ_ALPHABET_INCLUDED */

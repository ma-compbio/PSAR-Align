
/**
 * \file gff.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef SEQ_GFF_INCLUDED
#define SEQ_GFF_INCLUDED

#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "util/misc.h"

namespace fsa {

  /**
   * \brief Bare-bones representation of a single GFF entry.
   *
   * This class repeatedly assumes that the seqid field holds the
   * chromosome which the feature is on.
   * Note that GFF coordinates are always 1-based and fully-closed.
   * This class stores coordinates for GFF features accordingly
   * (in contrast to the 0-based indexing used throughout other
   * parts of this code).
   * See http://www.sequenceontology.org/gff3.shtml.
   */
  struct GFF {

    typedef std::map<std::string, std::vector<std::string>, fsa::Util::String_equal_ci> Attribute_table;  ///< (case-insensitive) map of key, value pairs

  public:

    std::string seqid;        ///< field 0
    std::string source;       ///< field 1
    std::string type;
    unsigned start;
    unsigned end;
    float score;
    char strand;
    unsigned phase;
    std::vector<std::string> attributes_ordering;     ///< sorted list of keys in attributes field
    Attribute_table attributes_map;                   ///< map of key, value pairs

    /**
     * \brief Default constructor.
     *
     * Initializes a feature of length 0 (nonsense coordinates [2, 1]).
     */
  GFF()
  : seqid (""), source (""), type (""),
      start (2), end (1),
      score (-1.), strand (GFF::unknown_strand), phase (3)
    { }

    /**
     * \brief Constructor.
     *
     * Initialize from a (possibly-scored) interval.
     */
    GFF (std::string seqid, unsigned start, unsigned end, float score = -1.)
    : seqid (seqid), source (""), type (""),
      start (start), end (end),
      score (score), strand (GFF::unknown_strand), phase (3)
    { }

    /**
     * \brief Set start coordinate.
     *
     * This deserves an accessor method in order to catch the case of negative start coordinates,
     * which are unfortunately often found in real annotation files.
     */
    void set_start (const int s);

    /**
     * \brief Length of feature.
     */
    inline size_t length() const;

    /**
     * \brief 5' end of feature.
     */
    inline unsigned five_prime_end() const;

    /**
     * \brief 3' end of feature.
     */
    inline unsigned three_prime_end() const;

    /**
     * \brief Add a (possibly-new) key, value pair to the attributes field.
     *
     * If values for the key already exist, then the new value is appended.
     */
    template<typename K, typename V>
      inline void add_value (const K& key, const V& value);

    /**
     * \brief Set a (possibly-new) key, value pair in the attributes field.
     *
     * If values for the key already exist, then all old values are cleared
     * and replaced with the single new value.
     */
    template<typename K, typename V>
      inline void set_value (const K& key, const V& value);

    /**
     * \brief Get the values for a key.
     * \return values, or empty vector if no such key
     * \see get_value
     */
    inline std::vector<std::string> get_values (const std::string& key) const;

    /**
     * \brief Get the first value for a key.
     * \return value, or empty string if no such key
     * \see get_values
     */
    inline std::string get_value (const std::string& key) const;

    /**
     * \brief Set name in attributes field.
     */
    inline void set_name (const std::string& name);

    /**
     * \brief Get name from attributes field.
     */
    inline std::string get_name() const;

    /**
     * \brief Set ID in attributes field.
     */
    inline void set_id (const std::string& id);

    /**
     * \brief Get ID from attributes field.
     */
    inline std::string get_id() const;

    /**
     * \brief Initialize from a GFF-formatted std::string.
     * \param strip_leading_chr strips the leading 'chr', if present, from the seqid field
     * possibly holding the chromosome name
     */
    void from_string (const std::string& str, const bool strip_leading_chr = true);

    /**
     * \brief Write to GFF-formatted string.
     */
    std::string to_string() const;

    /**
     * \brief Output operator.
     *
     * Prints undef_char for undefined fields.
     */
    friend std::ostream& operator<< (std::ostream& o, const GFF& gff) {
      o << gff.to_string() << endl;
      return o;
    }


    static const std::string comment_char;                    ///< comment character
    static const std::string undef_char;                      ///< character for undefined fields

    static const char unknown_strand;                         ///< character for unknown strand

    static const std::string attributes_split_char;           ///< split key, value pairs in attributes field
    static const std::string attributes_assign_char;          ///< assign value to key in attributes field
    static const std::string attributes_list_char;            ///< separate values for a single key in attributes field

    static const std::string key_id;                          ///< GFF key for ID
    static const std::string key_name;                        ///< GFF key for Name

    static const size_t default_flanking;                     ///< \see find_closest_feature

    /**
     * \brief Function object for binary comparison of GFF objects (sort by seqid, start, end).
     */
    struct GFF_less : std::binary_function<GFF, GFF, bool> {
    public:
      bool operator() (const GFF& l, const GFF& r) const {
	if (l.seqid == r.seqid) {
	  if (l.start == r.start)
	    return l.end < r.end;
	  return l.start < r.start;
	}
	return l.seqid < r.seqid;
      }
    };


  protected:
  
    /**
     * \brief Get string for the attributes field.
     */
    std::string get_attributes_string() const;

  private:

    /**
     * \brief Set start coordinate.
     */
    inline void set_start (const unsigned s);

    /**
     * \brief Parse attributes string into key, value pairs.
     */
    void parse_attributes_string (const std::string& str);

    static Regexp re_key_value;                                  ///< match key, value pairs in attributes field

  };

  /**
   * \brief Represent a GFF file.
   */
  struct GFF_database {

    /**
     * \brief Constructor.
     */
    GFF_database()
    : __maxlen (0) { }

    /**
     * \brief Load from a file.
     * \param strip_leading_chr strips the leading 'chr', if present, from the seqid field
     * possibly holding the chromosome name
     */
    void from_file (const std::string& filename, const bool strip_leading_chr = true);

    /**
     * \brief Write.
     */
    void write (std::ostream& o) const;

    /**
     * \brief Store an entry.
     */
    inline void store_entry (const GFF& gff);

    /**
     * \brief Append a GFF_database to this one.
     */
    void append (const GFF_database& gff_db);

    /**
     * \brief Set a (possibly-new) key, value pair in the attributes field for every entry.
     *
     * If values for the key already exist, then all old values are cleared
     * and replaced with the single new value.
     * \see GFF::set_value
     */
    template<typename K, typename V>
      inline void set_value (const K& key, const V& value);

    /**
     * \brief Create nicely-formatted numerical IDs for all entries.
     */
    void create_unique_ids();

    /**
     * \brief Sort entries in database.
     *
     * This must be called in order for intersect_genomic_interval to work properly.
     */
    void sort_entries();

    /**
     * \brief Find the GFF feature whose 5' end is closest to the passed interval.
     * \see find_closest_feature
     */
    GFF find_closest_feature_five_prime
    (const std::string& chromosome,
     unsigned start, unsigned end,
     const size_t flanking = GFF::default_flanking) const;

    /**
     * \brief Find the GFF feature whose 3' end is closest to the passed interval.
     * \see find_closest_feature
     */
    GFF find_closest_feature_three_prime
    (const std::string& chromosome,
     unsigned start, unsigned end,
     const size_t flanking = GFF::default_flanking) const;

    /**
     * \brief Find GFF features which intersect the passed interval.
     * 
     * Assumes that the seqid field holds the chromosome
     * and that the passed interval coordinates are 0-based
     * and fully closed.  Note the contrast with the 1-based
     * coordinates used for the GFF features themselves.
     * Note that the entries MUST be sorted in order for this to function properly.
     * \see GFF_database::sort_entries
     */
    GFF_database intersect_genomic_interval (const std::string& chromosome,
					     const unsigned start, const unsigned end) const;

    /**
     * \brief Maximum length of entry.
     */
    size_t maxlen() const { return __maxlen; }

    /**
     * \brief Average length of entry.
     */
    size_t meanlen() const;

    /**
     * \brief Median length of entry.
     */
    size_t medianlen() const;

    /**
     * \brief Number of entries.
     */
    size_t size() const { return __entries.size(); }

    /**
     * \brief Data access operator.
     */
    const GFF& operator[] (const size_t i) const {
      assert (i < __entries.size());
      return __entries[i];
    }


    /**
     * \brief Get iterator to start of __entries.
     */
    std::vector<GFF>::iterator begin() {
      return __entries.begin();
    }

    /**
     * \brief Get const_iterator to start of __entries.
     */
    std::vector<GFF>::const_iterator begin() const {
      return __entries.begin();
    }

    /**
     * \brief Get iterator to end of __entries.
     */
    std::vector<GFF>::iterator end() {
      return __entries.end();
    }

    /**
     * \brief Get const_iterator to end of __entries.
     */
    std::vector<GFF>::const_iterator end() const {
      return __entries.end();
    }

  private:

    /**
     * \brief Find the GFF feature whose 3' or 5' end is closest to the passed interval.
     *
     * The closest feature is defined as the feature whose 3' or 5' end
     * (specified by use_feature_five_prime) is closest to the centroid of the passed interval.
     * \param flanking search for features within this distance of the requested interval
     * \param use_feature_five_prime look for features whose 5' end is closest to the passed interval (false for 3')
     * \return closest feature, or empty feature if no feature within flanking distance
     * \see intersect_genomic_interval
     */
    GFF find_closest_feature (const std::string& chromosome,
			      unsigned start, unsigned end,
			      const size_t flanking = GFF::default_flanking,
			      const bool use_feature_five_prime = true) const;

    std::vector<GFF> __entries; ///< individual GFF entries in database

    size_t __maxlen;            ///< maximum length of an entry in the database
    bool __is_sorted;           ///< have the entries been sorted?

  };



  /****************************************
   * Function definitions.
   ****************************************/

  inline size_t GFF::length() const {
    // catch case of 0-length feature
    if (end < start)
      return 0;
    return end - start + 1;
  }

  inline unsigned GFF::five_prime_end() const {
    if (strand == '-')
      return end;
    return start;
  }

  inline unsigned GFF::three_prime_end() const {
    if (strand == '-')
      return start;
    return end;
  }

  inline void GFF::set_start (const unsigned s) {
    if (s >= 1)
      start = s;
    else {
      cerr << "Setting start coordinate " << s << " to 0." << endl;
      start = 0;
    }
  }

  template<typename K, typename V>
    inline void GFF::add_value (const K& key, const V& value) {

    // use stringstream to convert key and value to string
    std::string key_str = Util::to_string (key);
    std::string value_str = Util::to_string (value);

    // now store
    // if we already have this key, then append value
    std::map<std::string, std::vector<std::string> >::iterator value_vector = attributes_map.find (key_str);
    if (value_vector != attributes_map.end())
      value_vector->second.push_back (value_str);
    // else add new key and append value
    else {
      attributes_ordering.push_back (key_str);
      attributes_map[key_str].push_back (value_str);
    }

  }

  template<typename K, typename V>
    inline void GFF::set_value (const K& key, const V& value) {

    // use stringstream to convert key and value to string
    std::string key_str = Util::to_string (key);
    std::string value_str = Util::to_string (value);

    // now store
    // if we already have this key, then clear old value and record new value
    std::map<std::string, std::vector<std::string> >::iterator value_vector = attributes_map.find (key_str);
    if (value_vector != attributes_map.end()) {
      value_vector->second.clear();
      value_vector->second.push_back (value_str);
    }
    else {
      attributes_ordering.push_back (key_str);
      attributes_map[key_str].push_back (value_str);
    }

  }

  inline std::vector<std::string> GFF::get_values (const std::string& key) const {

    std::map<std::string, std::vector<std::string> >::const_iterator value_vector = attributes_map.find (key);
    if (value_vector != attributes_map.end())
      return value_vector->second;
    else
      return std::vector<std::string>();

  }

  inline std::string GFF::get_value (const std::string& key) const {

    std::map<std::string, std::vector<std::string> >::const_iterator value_vector = attributes_map.find (key);
    if (value_vector != attributes_map.end() && value_vector->second.size())
      return (value_vector->second)[0];
    else
      return "";

  }
  
  inline void GFF::set_name (const std::string& name) {
    set_value (key_name, name);
  }

  inline std::string GFF::get_name() const {
    const std::vector<std::string> names = get_values (key_name);
    if (!names.size())
      return "";
    else if (names.size() > 1) {
      cerr << "Warning: More than one name detected in GFF entry!" << endl
	   << *this << endl;
    }

    return names[0];
  }

  inline void GFF::set_id (const std::string& id) {
    set_value (key_id, id);
  }

  inline std::string GFF::get_id() const {
    const std::vector<std::string> ids = get_values (key_id);
    if (!ids.size())
      return "";
    else if (ids.size() > 1) {
      cerr << "Warning: More than one ID detected in GFF entry!" << endl
	   << *this << endl;
    }

    return ids[0];
  }

  inline void GFF_database::store_entry (const GFF& gff) {

    __entries.push_back (gff);
    __maxlen = (gff.end - gff.start + 1 > __maxlen) ? gff.end - gff.start + 1 : __maxlen;

  }

  template<typename K, typename V>
    inline void GFF_database::set_value (const K& key, const V& value) {

    for (std::vector<GFF>::iterator gff = begin(); gff != end(); ++gff)
      gff->set_value (key, value);
  }

}

#endif /* SEQ_GFF_INCLUDED */

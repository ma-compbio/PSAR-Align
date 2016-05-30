
/**
 * \file hash_fcn.h
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file implements Paul Hsieh's hash function,
 *  found at http://www.azillionmonkeys.com/qed/hash.html.
 */

#ifndef HASH_FCN_INCLUDED
#define HASH_FCN_INCLUDED

#include <cstring>
#include <config.h>
#include <stdint.h>

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__)	\
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)	\
		      +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

namespace fsa {

  typedef unsigned short bit16_t;
  typedef unsigned int bit32_t;
  typedef unsigned long long bit64_t;

  /// Assorted hash functions.
  // NB: Who knows why, but I get "duplicate symbol" errors
  // if these functions aren't wrapped up in a struct.
  struct Hash_functions {

    /// Hash for 16-bit integer.
    /*
     * \param key int to hash
     */
    static inline bit32_t bit16_t_hash (const bit16_t key) {
      return __lh3_Wang_hash_int (bit32_t (key));
    }

    /// Hash for 32-bit integer.
    /*
     * \param key int to hash
     */
    static inline bit32_t bit32_t_hash (const bit32_t key) {
      return __lh3_Wang_hash_int (key);
    }

    /// Hash for 64-bit integer.
    /*
     * \param key int to hash
     */
    static inline bit32_t bit64_t_hash (const bit64_t key) {
      return bit32_t (__lh3_Jenkins_hash_64 (key));
    }

    /// Hash for two 64-bit integers.
    /*
     * \param first first int to hash
     * \param second second int to hash
     *
     * XORs the two integers, then uses the Jenkins hash.
     */
    static inline bit32_t bit64_t_pair_hash (const bit64_t first, const bit64_t second) {
      return bit32_t (__lh3_Jenkins_hash_64 (first ^ second));
    }

    /// Wang's hash function for 32-bit integers.
    static inline bit32_t __lh3_Wang_hash_int (bit32_t key) {
      key += ~(key << 15);
      key ^=  (key >> 10);
      key +=  (key << 3);
      key ^=  (key >> 6);
      key += ~(key << 11);
      key ^=  (key >> 16);
      return key;
    }

    /// Jenkins' hash function for 64-bit integers.
    static inline bit64_t __lh3_Jenkins_hash_64 (bit64_t key) {
      key += ~(key << 32);
      key ^= (key >> 22);
      key += ~(key << 13);
      key ^= (key >> 8);
      key += (key << 3);
      key ^= (key >> 15);
      key += ~(key << 27);
      key ^= (key >> 31);
      return key;
    }

    /// Hsieh's hash.
    /*
     * \param data char array to hash
     * \param len length of char array
     * \see hsieh_hash_incr
     */
    static inline uint32_t hsieh_hash (const char* data, int len) {
      return hsieh_hash_incr (data, len, len);
    }

    /// Hsieh's hash.
    /*
     * \param data char array to hash
     * \see hsieh_hash_incr
     */
    static inline uint32_t hsieh_hash (const char* data) {
      int len = std::strlen (data);
      return hsieh_hash_incr (data, len, len);
    }

    /// Hsieh's hash for pairs.
    /*
     * \param first first element (char array) to hash
     * \param second second element (char array) to hash
     * \param len length of char array
     * Hash first, then use that to hash second.
     * Assumes that both first and second have identical lengths.
     * \see hsieh_hash_incr
     */
    static inline uint32_t hsieh_hash_pair (const char* first, const char* second, int len) {
      return hsieh_hash_incr (second, len, hsieh_hash_incr (first, len, 0));
    }

    /// Hsieh's hash for pairs.
    /*
     * \param first first element (char array) to hash
     * \param second second element (char array) to hash
     * Hash first, then use that to hash second.
     * Assumes that both first and second have identical lengths.
     * \see hsieh_hash_incr
     */
    static inline uint32_t hsieh_hash_pair (const char* first, const char* second) {
      return hsieh_hash_incr (second, std::strlen (second), hsieh_hash_incr (first, std::strlen (first), 0));
    }

    /// Incremental version of Hsieh's hash.
    /*
     * \param data char array to hash
     * \param len length of char array
     * \param hash some constant
     */
    static inline uint32_t hsieh_hash_incr (const char* data, int len, uint32_t hash) {
      uint32_t tmp;
      int rem;

      // the below previously was
      //   if (len <= 0 || data == NULL) return 0;
      // -- RKB 9/27/08
      if (len <= 0 || data == 0) return 0;

      rem = len & 3;
      len >>= 2;

      /* Main loop */
      for (;len > 0; len--) {
	hash  += get16bits (data);
	tmp    = (get16bits (data+2) << 11) ^ hash;
	hash   = (hash << 16) ^ tmp;
	data  += 2*sizeof (uint16_t);
	hash  += hash >> 11;
      }

      /* Handle end cases */
      switch (rem) {
      case 3: hash += get16bits (data);
	hash ^= hash << 16;
	hash ^= data[sizeof (uint16_t)] << 18;
	hash += hash >> 11;
	break;
      case 2: hash += get16bits (data);
	hash ^= hash << 11;
	hash += hash >> 17;
	break;
      case 1: hash += *data;
	hash ^= hash << 10;
	hash += hash >> 1;
      }

      /* Force "avalanching" of final 127 bits */
      hash ^= hash << 3;
      hash += hash >> 5;
      hash ^= hash << 4;
      hash += hash >> 17;
      hash ^= hash << 25;
      hash += hash >> 6;

      return hash;
    }

  };

}

#endif /* HASH_FCN_INCLUDED */

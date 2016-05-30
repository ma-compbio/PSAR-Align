/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef DB_MISC_INCLUDED
#define DB_MISC_INCLUDED

#include <string>
#include "util/sstring.h"

#define BUFFER_RATIO 0.9
#define MAX_MEM_BUFFER_SIZE 1024*1024*10 // MAX SOCKET SIZE = 10MB
#define MAX_DB_BUFFER_SIZE 1024*1024*20  // MAX DB SIZE = 1GB
#define MAX_DB_ENTRY_LENGTH 100

#define DB_FSA_SUFIX           "fsa" 
#define DB_FSA_TABLE_PREFIX    "root"

#define DB_SEQUENCE_INFO       "sequence_info"
#define DB_SPARSE_MATRIX_INFO  "sparse_matrix_info"
#define DB_HEAP_INFO           "heap_info"
#define DB_PARAMETER_INFO      "parameter_info"

#define DB_SPARSE_MATRIX_SUFIX "s"
#define DB_HEAP_SUFIX          "h"
#define DB_MERGED_HEAP_SUFIX   "heap"
#define DB_NUM_CELLS_SUFIX     "ncells"

#define DB_MEMORY_RATIO        1
#define DEFAULT_DB_BUFFERSIZE  ULONG_MAX

namespace fsa {

  class FSA;

  /// Three types of response
  typedef enum {
    DB_NOT_AVAILABLE = -1,
    DB_OK  = 1,
    DB_BAD = 0
  } DB_Response_Type;

  /// Yes or No type returns
  typedef enum {
    DB_YES = 1,
    DB_NO  = 0
  } DB_Yes_No_Type;

  /// the structure of num cells buffer
  /*
   * the sparse matrix of the sequence pair (seq1, seq2)
   * requires size data cells.
   */
  typedef struct Num_Cells_Buffer {
    int seq1;
    int seq2;
    int size;
  } Num_Cells_Buffer;

  /// the structure of sparse matrix buffer
  typedef struct Sparse_Matrix_Buffer {
    int pos1;
    int pos2;
    float prob;
  } Sparse_Matrix_Buffer;

  /// the structure of heap buffer
  typedef struct Heap_Buffer {
    int seq1;
    int pos1;
    int seq2;
    int pos2;
    double weight;
    float  delta;
  } Heap_Buffer;

  /// the structure of MEM_Buffer 
  /*  this is used when workers send the data including pairwise posterior probabilities and 
   *  candidate edges to the master directly.
   */
  typedef struct MEM_Buffer {
    Sparse_Matrix_Buffer *sparse_matrix;
    Num_Cells_Buffer *num_cells;
    Heap_Buffer *heap;
  } MEM_Buffer;

  /// the structure of DB_Buffer
  /*
   *  this is used when workers put the data into the database
   */
  typedef struct DB_Buffer {
    sstring	    sparse_matrix;
    sstring	    num_cells;
    sstring	    heap;
  } DB_Buffer;

  /// type define pairs
  typedef std::pair <int, Sparse_Matrix_Buffer *> Sparse_pair;
  typedef std::pair <int, Num_Cells_Buffer *> Cells_pair;
  typedef std::pair <int, Heap_Buffer *> Heap_pair;

  /// type define vectors
  typedef std::vector <Sparse_pair> Sparse_vector;
  typedef std::vector <Cells_pair> Cells_vector;
  typedef std::vector <Heap_pair> Heap_vector;

  /// the structure of MEM_Buffers
  typedef struct MEM_Buffers {
    Sparse_vector sparse_matrix_buffers;
    Cells_vector num_cells_buffers;
    Heap_vector heap_buffers;
  } MEM_Buffers;


  /// class DB_opts
  /*
   * The main purpose of this class is to store states of options that have been used when running FSA
   */
  struct DB_opts {

    const FSA *fsa;

    // alignment speedup options
    int num_alignment_pairs;             /// total number of all (n choose 2) pairs to consider during alignment inference

    // sequence annealing options
    int num_refinement_steps;            /// number of iterative refinement steps

    // anchoring options
    bool anchored;                       /// use anchor annealing

    // parameter estimation options
    bool learn_gap;                           /// learn indel parameters
    bool regularize;                          /// regularize learned parameters with Dirichlet distribution specified by model
    bool learn_emit_all;                      /// learn emit parameters over all sequences

    // parallelization options
    int  num_parallelized_jobs;       /// num of jobs to be simutaneously run

    // database connection options
    sstring db_hostname;           	  /// database server host name 
    sstring db_hostaddr;              /// database server host IP address
    sstring db_name;                  /// database name
    int 	db_port;                  /// database server port
    sstring db_user;                  /// database user name
    sstring db_password;              /// database password
    int  db_max_ram;                  /// database maximum ram

    // output options
    bool write_db;           /// write post. prob. matrices and cand. edges to database

    /// constructor
    DB_opts() {}

    /// Copy all parameters.
    void copy_opts (const FSA *from); 
  };

}

#endif

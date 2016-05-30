/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef MANAGER_INCLUDED 
#define MANAGER_INCLUDED 
#include <queue>
#include <float.h>

#include "seq/sequence.h"
#include "annealing/SparseMatrix.h"
#include "manager/db_adapter.h"
#include "manager/mw_adapter.h"

#define MAX_SPARSE_MATRIX_COUNT 10

#ifdef HAVE_POSTGRES
#include "manager/db_postgres.h"
#endif

namespace fsa {

  class Edge;
  class Column;
  class smaller_weight;

  /// Point a cell (which is SparseMatrix*) in a SparseMatrix 
  typedef SparseMatrix** Sparse_matrix_entry;

  /// Map a position in a sequence to the containing Column*.
  typedef std::vector<Column*> Seq_pos_col_map;

  // Class Manager
  /*
   * This class is needed to glue the fsa class together with some classes for the parallelization and the database
   */
  class Manager {

  public:

    /// Constructor
    Manager();

    /// Constructor
    /*
     * Detects whether the database is running, and we can get the appropriate data from the database
     */
    Manager (const Sequence_database& seq_db_internal, DB_opts &db_opts);

    /// Destructor
    ~Manager(); 				


    /// Get methods

    /// Get a sparse matrix from the database
    /*
     * Create a new sparse matrix for the sequence pair (i,j)
     */
    bool get_sparse_matrix(std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const int i, const int j); 

    /// Get numbers of data cells that are required for the sparse matrices
    bool get_num_cells (); 

    /// Create all sparse matrices.
    /*
     * If the database mode is being used, pairwise posterior probabilities are gotten from the database. 
     * If the parallizaition mode is being used, the data are transfered from the workers
     */
    bool get_all_sparse_matrices (std::vector<std::vector<SparseMatrix*> >& sparse_matrices); 


    /// Construct the priority queue with the candidate edges of the null alignment
    /*
     * If the database mode is being used, candidate edges are gotten from the database. 
     * If the parallizaition mode is being used, the data are transfered from the workers
     */
    bool get_all_edges (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges, 
			std::vector<Seq_pos_col_map> &seq_pos_col_maps); 

    /// Get the size of edges
    /*
     * In the case that the --db-maxram option is used:
     * The priority queue, edges, contains only some portion of the total edges so that 
     * the number of tuples in the merged heap table is needed to be considered. 
     * Otherwise, just return the size of the priority queue 
     */
    int get_edges_size(std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges); 

    /// Get edges
    /*
     * In the case that the --db-maxram option is used:
     * instead of getting all candidate edges of the null alignment at once, some portion of 
     * candidate edges are gotten from the merged heap table. 
     * Otherwise, get all candidate edges at once from either the database or the workers
     */
    void get_edges (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
		    std::vector<Seq_pos_col_map> &seq_pos_col_maps);

    /// Get some of candidate edges from the database
    bool get_next_edges (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
			 std::vector<Seq_pos_col_map> &seq_pos_col_maps);

    /// Get sparse matrices
    /*
     * In the case that the --db-maxram option is used:
     * Create the merged heap table and copy all tuples in the heap tables into the merged heap table. 
     * Note that sparse matrix is gotten from the database on-the-fly during the sequence annealing
     */	
    void get_sparse_matrices (std::vector<std::vector<SparseMatrix*> >& sparse_matrices); 

    /// Get the top element in the priority queue
    /*
     * In the case that the --db-maxram option is used:
     * once next some edges have been gotten from the database, return the top element in the priority queue.
     * Otherwise, just return the top element in the priority queue
     */
    Edge* get_next_top_edge (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
			     std::vector<Seq_pos_col_map> &seq_pos_col_maps); 

    /// Push edge
    /*
     * In the case that the --db-maxram option is used:
     * If the recomputed weight of the edge is greater than the m_min_edge_weight, the smallest
     * edge weight among edge weights in the priority queue, push it into the priority queue. 
     * Otherwise, put this edge into the merged heap table. This is because if the recomputed 
     * weight is less than the m_min_edge_weight, there might be an edge whose initial weight is
     * greater then the edge. 
     *
     * In the case that the --db-maxram option is NOT used:
     * push it into the priority queue
     */
    void push_edge(std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges, Edge *edge); 

    /// Is the appropriate data in the database?
    bool is_data_available (); 

    /// Can we use the sparse matrix of the sequence pair (i,j)?
    bool is_sparse_matrix_available (const int i, const int j); 

    /// Are sparse matrices avialable?
    bool is_sparse_matrices_available (); 

    /// Are edges available?
    bool is_edges_available (); 

    /// Get the list of available sparse matrices
    bool check_available_sparse_matrices (); 

    /// run the master instance
    bool mw_master_run (int argc, char** argv, const Params &params_seed, const Params &pseudocounts);

    /// run the worker instance
    bool mw_worker_run(int argc, char** argv);

    bool mw_single_worker_run (Params &params_seed, Params &pseudocounts); 

  private:

    const Sequence_database *m_seq_db_internal;        // hold all input sequence data (FSA's internal format)
    DB_opts *m_db_opts;                    // hold all options that have been used when running FSA

    DB_adapter m_db_adapter;               // database adapter class
    MW_adapter m_mw_adapter;               // master-worker adapter class

#ifdef HAVE_POSTGRES
  public:

    /// Find the worker id that has generated the sparse matrix of the sequence pair (i,j)
    int look_up_sparse_matrix_table_id (const int i, const int j);                     

    /// Re-set the maximum number of sparse matrices to be stored in the memory
    /*
     * Note that this method is called when the --db-maxram option is used
     */
    bool update_size (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges);

  private:

    bool m_last_edge_chunk;              // is this the last chunk of edges?

    int m_num_inserted_edges;            // the number of inserted re-weighted edges to the merged heap table duing the sequence anenaling
    int m_merged_heap_offset;            // merged heap offet to be used to get some portion of edges
    double m_min_edge_weight;            // the smallest edge weight among edges in the priority queue

    std::vector<std::vector<int> > *m_num_cells;                      // 2-dim vector to keep numbers of data cells
    std::vector<std::vector<int> > *m_available_sparse_matrices;      // 2-dim vector for the available sparse matrices: 
    //   The entry [i][j] stores the worker id that has generated the sparse matrix[i][j]

    std::queue <Sparse_matrix_entry> *m_sparse_matrix_scheduler; // queue to keep the maximum number of sparse matrices that are on memory
    int m_cnt_sparse_matrices;                              // how many sparse matrices are on the memory?

    std::vector<int> m_sparse_matrix_ref;		                // 1-dim vector to be used in look_up_sparse_matrix_table_id method

    int m_orig_edges_size;                                  // the number of candidate edges of the null alignment

    int HEAP_WINDOW_SIZE;                      // how many edges can we hold at a time?
    int SPARSE_MATRIX_MAX_SIZE;                // how many sparse matrices can we hold at a time? 

    sstring	    m_db_merged_heap_buffer;              // the merged_heap_buffer

    double avg_sparse_matrix_size;                    // the average sparse matrix size
#endif

#ifdef HAVE_CONDOR
  private:
    MEM_Buffers m_mem_buffers;                        // MEM buffers
#endif

  };

}

#endif



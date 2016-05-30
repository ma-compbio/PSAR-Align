/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef MW_MASTER_INCLUDED
#define MW_MASTER_INCLUDED

#include "MWDriver.h"
#include "seq/sequence.h"
#include "manager/mw_task.h"
#include "manager/db_misc.h"

#ifdef HAVE_POSTGRES
#include "manager/db_postgres.h"
#endif

namespace fsa {

  class Params;
  class DB_opts;

  /// Application Driver subclass derived from MWDriver 
  class MW_master : public MWDriver 
  {
  public:

    /// constructor
    MW_master();

    /// constructor
    MW_master(const Params &params_seed, const Params &pseudocounts, MEM_Buffers &mem_buffers,
	      const int num_seqs, const int num_seqs_pairs, const int num_jobs,
	      const int seqs_schema_id, const int params_table_id); 


    /// Get the info from the user.  Don't forget to get the worker_executable! 
    MWReturn get_userinfo(int argc, char **argv);

    /// Set up an array of tasks here 
    MWReturn setup_initial_tasks( int *, MWTask *** );

    /// What to do when a task finishes 
    MWReturn act_on_completed_task( MWTask * );

    /// Put things in the send buffer here that go to a worker 
    MWReturn pack_worker_init_data( void );

    /// Print the results	
    void printresults();

    /// Write out the state of the master to a file. 
    void write_master_state( FILE *fp );

    /// Read in the state from a file. 
    void read_master_state( FILE *fp );

    /// Just return a newly constructed application task 
    MWTask* gimme_a_task();

    /// Compute next candidate starting pair 
    void compute_next_starting_pair(int& seq_i, int& seq_j, int length, int total_seq_num);

    /// Get cpuinfo 
    bool get_proc_cpuinfo(char *arch); 

    bool do_raw_unpack (MW_task *tf, double&, double& );

    void pack_params(const Params &params); 

    void AddTask( MWTask* t);

  private:
		
    int m_num_parallelized_jobs;
    int m_num_seq_pairs;
    int m_num_of_pairs; 
    int m_num_remains;
    int m_num_seqs;
    int m_seqs_schema_id; 
    int m_params_table_id;

    const Params *m_params_seed;
    const Params *m_pseudocounts;
    int    m_argc;
    char** m_argv;

    MEM_Buffers *m_mem_buffers;
  };

}

#endif

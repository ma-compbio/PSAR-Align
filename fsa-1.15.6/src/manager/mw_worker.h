/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef MW_WORKER_INCLUDED
#define MW_WORKER_INCLUDED

#include "MWWorker.h"

#include "manager/db_adapter.h"
#include "manager/mw_task.h"

namespace fsa {

  class Params;

  class MW_worker : public MWWorker 
  {
  public:

    /// constructor
    MW_worker();

    /// constructor
    MW_worker(DB_adapter &db_adapter);

    /// destructor
    ~MW_worker();

    /// Benchmarking 
    double benchmark( MWTask *t );

    /// return a task
    MWTask* gimme_a_task();

    /// unpack init data 
    MWReturn unpack_init_data( void );

    /// unpack parameters
    void unpack_params(Params &params);

    /// do the real work for each task 
    void execute_task( MWTask * );

    void init_transfer(const int worker_id, int &raw_flag); 

  private:

    int m_seqs_schema_id;
    int m_params_table_id;
		
    // command arguments
    int    m_argc;                
    char **m_argv;



  };

}

#endif

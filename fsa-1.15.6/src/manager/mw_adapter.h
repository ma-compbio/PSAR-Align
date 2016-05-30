/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef MW_ADAPTER_INCLUDED
#define MW_ADAPTER_INCLUDED

#include "seq/sequence.h"

#include "manager/db_misc.h"
#include "manager/db_adapter.h"

#ifdef HAVE_CONDOR
#include "MW.h"
#include "MWSocketRC.h"
#include "MWDriver.h"

#include "manager/mw_master.h"
#include "manager/mw_worker.h"
#endif

namespace fsa {

  class Params;

  /// MW_adapter class
  /*
   * This class manages all stuffs related to the Master-Worker framework
   */
  class MW_adapter {
  private:

  public:	
    /// constructor
    MW_adapter() {}
		
    /// destructor
    ~MW_adapter() {}

    /// run the worker instance
    int worker_run(int argc, char** argv, DB_adapter &db_adapter); 

    /// run the master instance
    int master_run(int argc, char** argv, 
		   const Params &params_seed, const Params &pseudocounts, MEM_Buffers &mem_buffers, 
		   const int num_seqs, const int num_seqs_pairs, const int num_jobs,
		   const int seqs_schema_id, const int params_table_id); 

  };

}

#endif

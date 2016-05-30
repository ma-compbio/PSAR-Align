/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef MW_TASK_INCLUDED
#define MW_TASK_INCLUDED

#include "MWTask.h"

namespace fsa {

  /// MW_task class
  /*
   * The master and the workers communicate by passing this class
   */
  class MW_task : public MWTask 
  {
  public:
    /// constructor
    MW_task();

    /// Destructor 
    ~MW_task();


    /// App is required to implement the following functions. 
		
    /// The driver packs the input data via RMC, the data which will be sent to a worker. 
    void pack_work( void );

    /// The worker unpacks input data via RMC, need to allocate space for data 
    void unpack_work( void );

    /// The worker packs result data via RMC, the result will be sent back to driver 
    void pack_results( void );

    /// The driver unpacks result data via RMC 
    void unpack_results( void );

    /// The following functions have default implementation. 

    /// Print the task to stdout 
    void printself( int level = 70 );

    /// write checkpoint info per task, for each task haven't been finished 
    void write_ckpt_info( FILE *fp );

    /// read checkpoint info, in the order written into the file 
    void read_ckpt_info( FILE *fp );

    int m_start_seq_i;
    int m_start_seq_j;
    int m_prev_length_sum;
    int m_length;
    int m_worker_id;

    int m_raw_flag;
  };

}

#endif

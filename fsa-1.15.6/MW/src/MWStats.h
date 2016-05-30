/***************************Copyright-DO-NOT-REMOVE-THIS-LINE**
  *
  * Condor Software Copyright Notice
  * Copyright (C) 1990-2004, Condor Team, Computer Sciences Department,
  * University of Wisconsin-Madison, WI.
  *
  * This source code is covered by the Condor Public License, which can
  * be found in the accompanying LICENSE.TXT file, or online at
  * www.condorproject.org.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * AND THE UNIVERSITY OF WISCONSIN-MADISON "AS IS" AND ANY EXPRESS OR
  * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY, OF SATISFACTORY QUALITY, AND FITNESS
  * FOR A PARTICULAR PURPOSE OR USE ARE DISCLAIMED. THE COPYRIGHT
  * HOLDERS AND CONTRIBUTORS AND THE UNIVERSITY OF WISCONSIN-MADISON
  * MAKE NO MAKE NO REPRESENTATION THAT THE SOFTWARE, MODIFICATIONS,
  * ENHANCEMENTS OR DERIVATIVE WORKS THEREOF, WILL NOT INFRINGE ANY
  * PATENT, COPYRIGHT, TRADEMARK, TRADE SECRET OR OTHER PROPRIETARY
  * RIGHT.
  *
  ****************************Copyright-DO-NOT-REMOVE-THIS-LINE**/
/* stats.h

   Header file for statistics class.
*/

#ifndef STATISTICS_H
#define STATISTICS_H

#include "MWWorkerID.h"
#include "MW.h"
#include "MWList.h"

/** This class serves as a holder for worker classes before they go to die.
    It holds the workers until a makestats() call is made, at which
    point it generates some statistics, given the workers.  
*/
class MWStatistics {

 public:
    /// Default constructor
    MWStatistics();
    /// Destructor...
    ~MWStatistics();

  /** Literally, take the stats from this worker and store 
      them in the stats class.  This is done when the worker
      dies and during a checkpoint.
  */
    void gather ( MWWorkerID * );

  /** This method is called at the end of a run.  It prints out some
      statistical information on that run.
  */
    void makestats();

  /** This function is just like makestats(), except it return the relative
      statistics information to the caller. 
  */
  void get_stats(double *average_bench,
		 double *equivalent_bench,
		 double *min_bench,
		 double *max_bench,
		 double *av_present_workers,
		 double *av_nonsusp_workers,
		 double *av_active_workers,
		 double *equi_pool_performance,
		 double *equi_run_time,
		 double *parallel_performance,
		 double *wall_time,
		 MWList<MWWorkerID> *workerlist );

 /** Everytime we get a new bench factor, we update the max and min ones if necessary */
  void update_best_bench_results( double bres ); 

  /** Each time a task finishes, we call this method, to record the 
    * best of times, the worst of times.
    */
  void update_best_task_results(double cpu_time);
  /** We find that we should save stats across checkpoints.  Here
      we're given the ckpt FILE* and the list of running workers.
      Basically, we write out the sum of the old stats we've got
      plus the stats of the running workers. */
  void write_checkpoint( FILE *cfp, MWList<MWWorkerID> *wlist );
  
  /** Read in stats information across the checkpoint */
  void read_checkpoint( FILE *cfp );

 private:

  double duration;
  double starttime;
  double previoustime; // the total run time of the checkpoints before...

  /// Array of times that the workerIDs have been "up"
  double *uptimes;
  /// Array of times that the workerIDs have been working
  double *workingtimes;
  /// Array of times holding the CPU times of the workerIDs
  double *cputimes;
  /// Array of times that the workerIDs have been suspended
  double *susptimes;

  /**@name Normalized Statistics.  MW computes "normalized" statistics
     about the runs performance based on the relative speed of the machines
     that the user obtained during the course of the run.
  */
     
     
  //@{

  double *workingtimesnorm; // the total working time normalized by benchmark factor

  /** Min benchmark result -- the "slowest" machine we ever had. */
  double max_bench_result;

  /** Max benchmark result -- the "faster" machine we ever had. */
  double min_bench_result;

  /** Min task result -- the "slowest" machine we ever had. */
  double max_task_result;

  /** Max benchmark result -- the "faster" machine we ever had. */
  double min_task_result;

  /// The sum of benchmark factors for this vid
  double *sumbenchmark;     
  /// The number of benchmar factors for this vid
  int *numbenchmark;        

  //@}

  /// Computes the mean of an array of number.
  double mean ( double *array, int length );

  /** Computes the variance of an array of numbers (given the mean, 
      to avoid "double processing" the list -- which you don't need
      to do anyway, but Mike wrote this function).
  */
  double var ( double *array, int length, double mean );

  /** What is the maximum "virtuid ID" that we ever assigned to a worker.
      This amounts to asking how many workers did we ever have in ALL states
      at any one time during the course of the run.
  */
  int get_max_vid();

};

#endif



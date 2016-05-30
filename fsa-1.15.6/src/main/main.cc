
/**
 * \file main.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley and Jaeyoung Do.
 */

#include "fsa/fsa.h"

using namespace fsa;

/**
 * \brief Check whether this appears to be a worker instance.
 *
 * Examines command-line arguments and checks whether they conform
 * to the expected format for worker instances.
 */
bool is_mw_worker_format (int argc, char** argv); 

/**
 * \brief Create & run an FSA object from command-line arguments
 */
int main (int argc, char** argv) {

  if (is_mw_worker_format (argc, argv))  {
    Manager manager;
    return manager.mw_worker_run (argc, argv);
  }
  else { 
    FSA fsa (argc, argv);
    return fsa.run();
  }
}

/**
 * \brief Simple method to check whether we seem to be running as a MW_worker.
 */
bool is_mw_worker_format(int argc, char** argv) {
	
  if (argc == 5) {
    bool chk_arg2 = (atoi(argv[2]) != 0 ) ? true : false; // 2nd parameter must be an integer value
    bool chk_arg3 = (atoi(argv[3]) != 0 ) ? true : false; // 3rd parameter must be an integer value

    int dot_cnt = 0;
    char tmp_argv[512],*dot_pch;

    strcpy(tmp_argv, argv[4]);
    dot_pch = strtok (tmp_argv, ".");

    while (dot_pch != NULL) {    // check whether 4th parameter is a valid IPv4
      ++dot_cnt;
      dot_pch = strtok(NULL, ".");
    }

    if (chk_arg2 && chk_arg3 && dot_cnt == 4) 
      return true;
  }	
  return false;
}

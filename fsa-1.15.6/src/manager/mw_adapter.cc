/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */
#include "manager/mw_adapter.h"

using namespace fsa;

int MW_adapter::worker_run(int argc, char** argv, DB_adapter &db_adapter) {
#ifdef HAVE_CONDOR
	MWWorker::RMC = new MWSocketRC ( FALSE, 0 );
	MWTask::RMC = MWWorker::RMC;
	MWDriver::RMC = NULL;

	MW_worker mw_worker (db_adapter);
	mw_worker.go (argc, argv); //start worker 
#endif
	return 0;
}

int MW_adapter::master_run(int argc, char** argv, const Params &params_seed, const Params &pseudocounts, MEM_Buffers &mem_buffers, 
		const int num_seqs, const int num_seqs_pairs, const int num_jobs,
		const int seqs_schema_id, const int params_table_id) {
#ifdef HAVE_CONDOR
	MW_master mw_master (params_seed, pseudocounts, mem_buffers, num_seqs, num_seqs_pairs, num_jobs, seqs_schema_id, params_table_id);
	mw_master.go (argc, argv); //start master
#endif
	return 0;
}



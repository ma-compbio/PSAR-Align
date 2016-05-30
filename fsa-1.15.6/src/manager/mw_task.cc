/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include "manager/mw_task.h"

using namespace fsa;

/// this funtion is called inside pack_results
extern void transfer_data(MW_task *);

MW_task::MW_task() {
	m_raw_flag = 0;
}

MW_task::~MW_task() {}

void MW_task::printself( int level ) {
	// do nothing
}

void MW_task::pack_work( void ) {

	RMC->pack(&m_worker_id,1,1);                // the number of sequence pairs
	RMC->pack(&m_prev_length_sum,1,1);
	RMC->pack(&m_length,1,1);                   // the number of sequence pairs
	RMC->pack(&m_start_seq_i,1,1);              // the first sequence of the starting sequence pair
	RMC->pack(&m_start_seq_j,1,1);              // the secnod sequence of the starting sequence pair
}

void MW_task::unpack_work( void ) {
	RMC->unpack(&m_worker_id,1,1);              // the number of sequence pairs
	RMC->unpack(&m_prev_length_sum,1,1);
	RMC->unpack(&m_length,1,1);                 // the number of sequence pairs
	RMC->unpack(&m_start_seq_i,1,1);            // the first sequence of the starting sequence pair
	RMC->unpack(&m_start_seq_j,1,1);            // the second sequence of the starting sequence pair

}

void MW_task::pack_results( void ) {
	RMC->pack(&m_raw_flag, 1, 1);

	if (m_raw_flag) {
		RMC->pre_send(RESULTS); 
	}

	transfer_data(this);
}

void MW_task::unpack_results( void ) {
	m_raw_flag = 0;
	RMC->unpack(&m_raw_flag, 1, 1);
}

void MW_task::write_ckpt_info( FILE *fp )  {
	// Nothing in this app, will lose data if it crashes. 
}

void MW_task::read_ckpt_info( FILE *fp ) {
	// Nothing to be read since nothing is written 
}

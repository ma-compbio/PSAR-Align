/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include "manager/db_misc.h"
#include "fsa/fsa.h"

using namespace fsa;

void DB_opts::copy_opts (const FSA *from) {

	// copy all options 

	fsa = from;	
	
	num_refinement_steps = from->num_refinement_steps;
	anchored = from->anchored;
	
	learn_emit_all = from->learn_emit_all;
	learn_gap = from->learn_gap;                         
	regularize = from->regularize;                      
	num_parallelized_jobs = from->num_parallelized_jobs;       
	num_alignment_pairs = from->num_alignment_pairs;
	
	db_hostname = from->db_hostname;
	db_hostaddr = from->db_hostaddr;
	db_name = from->db_name;
	db_port = from->db_port;
	db_user = from->db_user;
	db_password = from->db_password;
	db_max_ram = from->db_max_ram;

	write_db = from->write_db;
}


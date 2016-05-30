/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include <arpa/inet.h>
#include <pthread.h>
#include <string.h>

#include "manager/mw_worker.h"
#include "fsa/fsa.h"

using namespace fsa;

extern FSA *m_fsa;
extern DB_adapter *m_db_adapter;
extern Params m_params_seed;
extern Params m_pseudocounts;
extern MEM_Buffer m_mem_buffer;
extern DB_Buffer m_db_buffer;
extern bool final_buffer;
extern bool send_thread_waiting;



MW_worker::MW_worker() {}

MW_worker::MW_worker(DB_adapter &db_adapter) 
{
	// how much info you wants the driver to print 
	set_MWprintf_level( 1 );

	m_db_adapter = &db_adapter;

	workingTask = new MW_task;

	m_mem_buffer.sparse_matrix = NULL;
	m_mem_buffer.num_cells = NULL;
	m_mem_buffer.heap = NULL;
}

MW_worker::~MW_worker() {
	for (int i=0; i<m_argc; i++)
		free (m_argv[i]);
	free (m_argv);

	if (workingTask)
		delete workingTask;
}

double MW_worker::benchmark( MWTask *t ) {
	// As we don't have benchmark, we do nothing
	return (double)3.14;
}

void MW_worker::unpack_params (Params &params) {
	// bandwidth 
	RMC->unpack(&(params.bandwidth), 1,1);

	// is_indel2 
	char is_indel2;
	RMC->unpack(&is_indel2, 1, 1);
	params.is_indel2 = (is_indel2 == 1) ? true : false;

	// gap_open1 
	RMC->unpack(&(params.gap_open1), 1,1);

	// gap_open2 
	RMC->unpack(&(params.gap_open2), 1,1);

	// gap_extend1 
	RMC->unpack(&(params.gap_extend1), 1,1);

	// gap extend2 
	RMC->unpack(&(params.gap_extend2), 1,1);

	// to_end 
	RMC->unpack(&(params.to_end), 1,1);

	// time 
	RMC->unpack(&(params.time), 1,1);

	// alphabet_string 
	char *temp_char =  (char *) malloc (1024 * sizeof (char) );
	RMC->unpack(temp_char);
	params.alphabet_string = temp_char;
	free (temp_char);

	// single_dist: size 
	int size;
	RMC->unpack(&size, 1,1);
	params.single_dist.resize(size);

	// single_dist: values 
	for (int i=0;i<size;i++) 
		RMC->unpack(&(params.single_dist[i]), 1,1);

	// pair_dist: size(n) -> n by n array 
	RMC->unpack(&size, 1,1);
	params.pair_dist.resize(size);

	// pair_dist: values 
	for (int i=0;i<size;i++) { 
		params.pair_dist[i].resize(size);
		for (int j=0;j<size;j++)
			RMC->unpack(&(params.pair_dist[i][j]), 1,1);
	}

	// transaction_matrix: column size 
	RMC->unpack(&size, 1,1);
	params.transition_matrix.resize(size);

	for (int i=0;i<size;i++)  {
		params.transition_matrix[i].resize(size);
		for (int j=0;j<size;j++) 
			RMC->unpack(&(params.transition_matrix[i][j]), 1,1);
	}
}

/// unpack the init data from the driver 
MWReturn MW_worker::unpack_init_data( void )  {

	// Read the total number of arguments 
	RMC->unpack(&m_argc,1,1);

	// Read arguments 
	m_argv = (char**) malloc(m_argc * sizeof(char*));

	for (int i = 0; i < m_argc; i++) {
		m_argv[i] = (char *) malloc (512 * sizeof(char));
		RMC->unpack(m_argv[i]);
	}

	// Read schema id of DB 
	RMC->unpack(&m_seqs_schema_id, 1,1);

	// Read data id of DB 
	RMC->unpack(&m_params_table_id, 1,1);

	unpack_params (m_params_seed);
	unpack_params (m_pseudocounts);

	return OK;
}

void MW_worker::init_transfer(const int worker_id, int &raw_flag) {

	final_buffer = false;
	send_thread_waiting = false;

	if (m_fsa->write_db) {
		bool is_db_running = false;
		m_db_buffer.sparse_matrix.clear();
		m_db_buffer.num_cells.clear();
		m_db_buffer.heap.clear();

		// connect DB
		is_db_running = m_db_adapter->connect_db(m_fsa->db_hostname.c_str(), m_fsa->db_hostaddr.c_str(), m_fsa->db_name.c_str(), 
				m_fsa->db_port, m_fsa->db_user.c_str(), m_fsa->db_password.c_str());

		if (!is_db_running)
			THROWEXPR ("ERROR: We can not access database server.");

		m_db_adapter->set_up_ids (m_seqs_schema_id, m_params_table_id, worker_id);

		raw_flag = 0;

	}
	else {
		if (!m_mem_buffer.sparse_matrix) 
			m_mem_buffer.sparse_matrix = (Sparse_Matrix_Buffer*) malloc (
					sizeof (Sparse_Matrix_Buffer) * ((int) (MAX_DB_BUFFER_SIZE / sizeof (Sparse_Matrix_Buffer))));
		

		if (!m_mem_buffer.num_cells)	
			m_mem_buffer.num_cells = (Num_Cells_Buffer *) malloc ( 
					sizeof (Num_Cells_Buffer) * ((int) (MAX_DB_BUFFER_SIZE / sizeof (Num_Cells_Buffer))));

		if (!m_mem_buffer.heap)
			m_mem_buffer.heap = (Heap_Buffer *) malloc ( 
					sizeof (Heap_Buffer) * ((int) (MAX_DB_BUFFER_SIZE / sizeof (Heap_Buffer))));

		raw_flag = 1;

	}
}

/// Execute each task 
void MW_worker::execute_task( MWTask *t ) {

	MW_task *tf = dynamic_cast<MW_task *>(t);
	int num_seq_pairs = tf->m_length;
	int prev_seq_pairs_sum = tf->m_prev_length_sum;
	std::pair <int,int> start_seq_pair (tf->m_start_seq_i, tf->m_start_seq_j);

	// create fsa
	m_fsa = new FSA (m_argc, m_argv);

	// create buffers and init
	m_fsa->init_for_mw_worker (start_seq_pair, prev_seq_pairs_sum, num_seq_pairs);

	try
	{
		// It's important that we input the data /before/ parsing the command-line 
		// This allows us to set preset values for DNA, RNA and proteins and then let
		// them be overriden by the command-line 
		m_fsa->init_opts();
		m_fsa->input_data();
		m_fsa->set_up_defaults();
		m_fsa->parse_opts();
		m_fsa->assemble_sequence_data();
		m_fsa->choose_seq_pairs();

		init_transfer(tf->m_worker_id, tf->m_raw_flag);
	}
	catch (const Dart_exception& e)
	{
		CLOGERR << e.what();
		exit(1);
	}

}

MWTask* MW_worker::gimme_a_task() {
	return new MW_task;
}

// Just return a newly created application worker object 
MWWorker* gimme_a_worker () {
	return new MW_worker;
}


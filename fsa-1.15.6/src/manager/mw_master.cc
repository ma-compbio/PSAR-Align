/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

// TODO int meanlen = (int) floor (seq_db.mean_length());

#include <arpa/inet.h>
#include <pthread.h>

#include "manager/mw_master.h"
#include "manager/db_adapter.h"
#include "fsa/fsa.h"

using namespace fsa;

static pthread_mutex_t worker_thread_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct ThreadData {
	MW_task* m_task;
	MW_master* m_driver;
} ThreadData;

MW_master::MW_master(const Params &params_seed, const Params &pseudocounts, MEM_Buffers &mem_buffers,
		const int num_seqs, const int num_seqs_pairs, const int num_jobs,
		const int seqs_schema_id, const int params_table_id) {

	set_MWprintf_level (1);

	m_num_parallelized_jobs = num_jobs;
	m_params_seed = &params_seed;
	m_pseudocounts = &pseudocounts;
	m_num_seqs = num_seqs;
	m_num_seq_pairs = num_seqs_pairs;

	m_seqs_schema_id = seqs_schema_id;
	m_params_table_id = params_table_id;

	m_mem_buffers = &mem_buffers;
}


static void *__do_raw_unpack(void *ptr) 
{
	ThreadData* tdata = (ThreadData*)ptr;
	MW_task* t = tdata->m_task;
	MW_master* driver = tdata->m_driver;

	int worker_id = t->worker->get_id1();
	free(tdata);

	double wall_time = 0.0;
	double cpu_time = 0.0;
	bool ret = driver->do_raw_unpack(t, wall_time, cpu_time);

	// Now we have finished reading data from socket
	// So we transfer the control over the socket to main process
	t->RMC->thread_stop(worker_id);

	if( !ret ) {
		fprintf(stderr, "do_raw_unpack failed from worker(id=%d,name=%s)\n", 
				worker_id, t->worker->machine_name());
		fprintf(stderr, "Creating a new task for this failed task\n");
		fflush(stderr);
		// Fail
		//in case of error, need to put this task into job queue again..
		
		// Create a new Task
		MW_task *new_task = new MW_task;
		assert(new_task);

		new_task->m_worker_id = t->m_worker_id;
		new_task->m_prev_length_sum = t->m_prev_length_sum;
		new_task->m_length = t->m_length;
		new_task->m_start_seq_i = t->m_start_seq_i;
		new_task->m_start_seq_j = t->m_start_seq_j;

		driver->getGlobalLock();

		// Delete old one
		MWWorkerID* w = t->worker;
		wall_time = 0.0;
		cpu_time = 0.0;
		t->completedTask(wall_time, cpu_time);
		delete t;

		// Add new one
		driver->AddTask((MWTask*)new_task);

		driver->getGlobalUnLock();
	}else {
		// success 
		driver->getGlobalLock();
		MWWorkerID* w = t->worker;
		t->completedTask(wall_time, cpu_time);
		delete t;
		driver->getGlobalUnLock();
	}

	driver->decreasePendingNum();
	return NULL;
}

/// get_userinfo 
MWReturn MW_master::get_userinfo(int argc, char **argv) 
{

	char arch[10];
	char requirements[100];

	char *exec_ptr = NULL;

	m_argc = argc;
	m_argv = argv;

	/* exec classes */
	RMC->set_num_exec_classes(1);
	/* arch classes */
	RMC->set_num_arch_classes(1);
	/* set architecture and linux type */

	if (get_proc_cpuinfo(arch)) {
		sprintf(requirements,"((Arch==\"%s\") && (Opsys==\"LINUX\"))",arch);
		RMC->set_arch_class_attributes (0, requirements); 
	}
	else 
		RMC->set_arch_class_attributes (0, "((Arch==\"INTEL\") && (Opsys==\"LINUX\"))");  // default setting

	/* we have only one executable - fsa */
	RMC->set_num_executables(1);

	RMC->add_executable(0, 0, (((exec_ptr = strchr(argv[0], '/')) != NULL)? (exec_ptr + 1) : argv[0]), "");  

	/* checkpoint requirement */
	set_checkpoint_frequency(0);

	/* Set the number of jobs */

	RMC->set_target_num_workers(m_num_parallelized_jobs);
	RMC->set_worker_increment(m_num_parallelized_jobs);

	m_num_of_pairs = m_num_seq_pairs / m_num_parallelized_jobs; 
	m_num_remains = m_num_seq_pairs % m_num_parallelized_jobs;

	return OK;
}

/// setup (generate and push) the first batch of tasks in the beginning 
MWReturn MW_master::setup_initial_tasks(int *n_init , MWTask ***init_tasks) 
{
	// Create initial tasks
	*n_init = m_num_parallelized_jobs;
	*init_tasks = new MWTask* [m_num_parallelized_jobs];

	CTAG(9, MW) << "Creating " << m_num_parallelized_jobs << " job(s)." << endl;

	// Starting pair of sequence
	int seq_i = 0;
	int seq_j = 1;
	int i = 0;
	int prev_length_sum = 0;

	for (i = 0; i < m_num_parallelized_jobs; i++) {

		MW_task *new_task = new MW_task;

		int prev_length = (i == 0) ? 0 : (((i - 1) < m_num_remains) ? (1 + m_num_of_pairs) : (m_num_of_pairs));

		int length = (i < m_num_remains) ? m_num_of_pairs + 1 : m_num_of_pairs;
		prev_length_sum += prev_length;

		compute_next_starting_pair(seq_i, seq_j, prev_length, m_num_seqs);

		new_task->m_worker_id = i;
		new_task->m_prev_length_sum = prev_length_sum;
		new_task->m_length = length;
		new_task->m_start_seq_i = seq_i;
		new_task->m_start_seq_j = seq_j;

		(*init_tasks)[i] = new_task;
	}

	return OK;
}

void MW_master::pack_params(const Params &params) {

	// bandwidth 
	RMC->pack(&(params.bandwidth), 1,1);

	// is_indel2 
	char is_indel2 = (params.is_indel2)	? 1 : 0;
	RMC->pack(&is_indel2, 1,1);

	// gap_open1 
	RMC->pack(&(params.gap_open1), 1,1);

	// gap_open2 
	RMC->pack(&(params.gap_open2), 1,1);

	// gap_extend1 
	RMC->pack(&(params.gap_extend1), 1,1);

	// gap extend2 
	RMC->pack(&(params.gap_extend2), 1,1);

	// to_end 
	RMC->pack(&(params.to_end), 1,1);

	// time 
	RMC->pack(&(params.time), 1,1);

	// alphabet_string 
	RMC->pack(params.alphabet_string.c_str());

	// single_dist: size 
	int size = (int) params.single_dist.size();
	RMC->pack(&size, 1,1);

	// single_dist: values 
	for (int i=0;i<(int) params.single_dist.size();i++) 
		RMC->pack(&(params.single_dist[i]), 1,1);

	// pair_dist: size(n) -> n by n array 
	size = (int) params.pair_dist.size();
	RMC->pack(&size, 1,1);

	// pair_dist: values 
	for (int i=0;i<(int) params.pair_dist.size();i++) 
		for (int j=0;j<(int) params.pair_dist.size();j++)
			RMC->pack(&(params.pair_dist[i][j]), 1,1);

	// transition_matrix: column size 
	size = (int) params.transition_matrix.size();
	RMC->pack(&size, 1,1);

	for (int i=0;i<(int) params.transition_matrix.size();i++) 
		for (int j=0;j<(int) params.transition_matrix[i].size();j++) 
			RMC->pack(&(params.transition_matrix[i][j]), 1,1);
}

/// The first batch of data for a newly spawned worker, e.g. init data 
MWReturn MW_master::pack_worker_init_data() {

	int argc=0;
	char argv[100][512];

	// Send arguments without the arugment having information of the number of jobs
	for (int i=0; i<m_argc; i++)  {
		if (strstr(m_argv[i], "--parallelize") != NULL) 
			++i;
		else {
			strcpy(argv[argc], m_argv[i]);
			++argc;
		}
	}

	strcpy(argv[argc++], "--noannealing");


	// Send the total number of arguments 
	RMC->pack(&argc,1,1);

	for (int i=0; i<argc; i++) 
		RMC->pack(argv[i]); 

	// Send the schema id of DB 
	RMC->pack(&m_seqs_schema_id, 1,1);

	// Send the data id of DB 
	RMC->pack(&m_params_table_id, 1,1);

	pack_params (*m_params_seed);
	pack_params (*m_pseudocounts);

	return OK;
}

static bool raw_unpack(MWRMComm* _RMC, int worker_id, char* buf, int len)
{
	static int cnt = 0;
	cnt++;

	int got = _RMC->raw_unpack(worker_id, buf, len);
	if( got != len ) {
		fprintf(stderr, "raw_unpack[%d] failed to read the size of data.\n", cnt);
		fprintf(stderr, "need to read(%d) but got(%d)\n", len, got);
		fflush(stderr);
		return false;
	}
	return true;
}

bool MW_master::do_raw_unpack(MW_task *tf, double& wall_time, double& cpu_time)
{
	int tmpSize = 0;
	int hostSize = 0;

	int worker_id = tf->worker->get_id1();
	MWRMComm* _RMC = tf->RMC;

	Sparse_Matrix_Buffer *sparse_matrix_buffer = NULL;
	Num_Cells_Buffer *num_cells_buffer = NULL;
	Heap_Buffer *heap_buffer = NULL; 

	int what_data = 0;
	int num_edges = 0;
	int num_probs = 0;
	int num_logs = 0;

	MEM_Buffers mem_buffers;	

	do {
		if (!raw_unpack( _RMC, worker_id, (char*)&what_data, sizeof(int)))
			return false;

		if ( what_data == 1 ) {

			if (!raw_unpack( _RMC, worker_id, (char*)&num_logs, sizeof(int))) {
				return false;
			}

			num_cells_buffer = (Num_Cells_Buffer *) malloc (sizeof (Num_Cells_Buffer) * num_logs);
			assert(num_cells_buffer);


			if (!raw_unpack( _RMC, worker_id, (char*) num_cells_buffer, num_logs * sizeof(Num_Cells_Buffer))) {
				return false;
			}

			if (!raw_unpack( _RMC, worker_id, (char*)&num_probs, sizeof(int))) {
				return false;
			}

			sparse_matrix_buffer = (Sparse_Matrix_Buffer *) malloc (sizeof (Sparse_Matrix_Buffer) * num_probs);
			assert(sparse_matrix_buffer);

			if (!raw_unpack( _RMC, worker_id, (char*) sparse_matrix_buffer, num_probs * sizeof(Sparse_Matrix_Buffer))) {
				return false;
			}
			
			mem_buffers.num_cells_buffers.push_back (Cells_pair (num_logs, num_cells_buffer) );
			
			mem_buffers.sparse_matrix_buffers.push_back (Sparse_pair (num_probs, sparse_matrix_buffer) );


		}
		else if ( what_data == 2) {

			if (!raw_unpack( _RMC, worker_id, (char*)&num_edges, sizeof(int))) {
				return false;
			}

			heap_buffer = (Heap_Buffer *) malloc (sizeof (Heap_Buffer) * num_edges);
			assert(heap_buffer);

			if (!raw_unpack( _RMC, worker_id, (char*) heap_buffer, num_edges * sizeof(Heap_Buffer))) {
				return false;
			}

			mem_buffers.heap_buffers.push_back (Heap_pair (num_edges, heap_buffer) );
		}

	} while ( what_data > 0 );

	// We finished to read all data
	// Now, we need to read two doubles
	wall_time = 0.0;
	cpu_time = 0.0;
	if (!raw_unpack( _RMC, worker_id, (char*)&wall_time, sizeof(double))) {
		return false;
	}
	if (!raw_unpack( _RMC, worker_id, (char*)&cpu_time, sizeof(double))) {
		return false;
	}

	pthread_mutex_lock(&worker_thread_mutex);	

	for (int i = 0; i < mem_buffers.num_cells_buffers.size(); i++) {
		m_mem_buffers->num_cells_buffers.push_back (mem_buffers.num_cells_buffers[i]);
	}

	for (int i = 0; i < mem_buffers.sparse_matrix_buffers.size(); i++) {
		m_mem_buffers->sparse_matrix_buffers.push_back (mem_buffers.sparse_matrix_buffers[i]);
	}

	for (int i = 0; i < mem_buffers.heap_buffers.size(); i++) {
		m_mem_buffers->heap_buffers.push_back (mem_buffers.heap_buffers[i]);
	}
	
	pthread_mutex_unlock(&worker_thread_mutex);	


	return true;
}

/// Implement application behavior to process a just completed task */
MWReturn MW_master::act_on_completed_task( MWTask *t ) {
	MW_task *tf = dynamic_cast<MW_task *>(t);

	if( !tf->m_raw_flag ) {
		return OK;
	}

	/// It means that we need to call raw_unpack
	int worker_id = tf->worker->get_id1();

	// Set the timeout of Main process select/poll to 1 sec
	tf->RMC->setPollTimeOut(5);
		
	// This will indicate 
	// "this socket will be taken over by new thread"
	// "So main thread MUST do nothing on it"
	tf->RMC->thread_start(worker_id);

	// Create pthread
	increasePendingNum();

	ThreadData* tdata = (ThreadData*)malloc(sizeof(ThreadData));
	assert(tdata);

	tdata->m_task = tf;
	tdata->m_driver = this;

	pthread_t thread;
	if( pthread_create(&thread, NULL, __do_raw_unpack, (void *)tdata) != 0 ) {
		fprintf(stderr, "Failed to create a new thread\n");	
		free(tdata);

		tf->RMC->thread_stop(worker_id);

		decreasePendingNum();

		double wall_time = 0.0;
		double cpu_time = 0.0;
		if( do_raw_unpack(tf, wall_time, cpu_time) == false ) {
			fprintf(stderr, "do_raw_unpack failed from worker(id=%d,name=%s)\n", 
					worker_id, tf->worker->machine_name());
			return ABORT;
		}
	}else {
		// Detach this thread
		pthread_detach(thread);
	}

	return PENDING;
}


/// printresults
void MW_master::printresults() {

	assert(getPendingNum() == 0 );

}

/// Write app-specific master checkpoint info 
void MW_master::write_master_state( FILE *fp ) {
	// Nothing to be written 
}

/// Read app-specific master checkpoint info 
void MW_master::read_master_state( FILE *fp ) {
	// Nothing to be read 
}

/// Return a new application task object 
MWTask* MW_master::gimme_a_task() {
	return new MW_task;
}

/// Return a new driver object 
MWDriver* gimme_the_master() {
	return NULL;
}

void 
MW_master::AddTask( MWTask* t)
{
	addTask(t);
}

/// compute_next_starting_pair 
void MW_master::compute_next_starting_pair(int& seq_i, int& seq_j, int length, int total_seq_num) {

	int tmp_length = length;

	// Compute the next sequence pair 
	while( tmp_length > 0 ) {
		if( (seq_j + tmp_length) < total_seq_num ) {
			seq_j += tmp_length;
			break;
		} else {
			tmp_length -= total_seq_num - seq_j;
			seq_j = ++seq_i + 1;
		}
	}
}

/// get_proc_cpuinfo 
bool MW_master::get_proc_cpuinfo(char *arch) {
	char path[64];
	FILE *fp;
	char line_input[1024];

	memset(path, 0, sizeof(path));
	strcpy(path, "/proc/cpuinfo");

	if ((fp = fopen(path, "r")) != NULL ) {
		while (!feof(fp)) {
			fgets ( line_input, 1024, fp);
			if (strstr ( line_input, "flags" ) != NULL && strstr(line_input, " lm ") != NULL) {
				strcpy(arch,"X86_64");
				fclose ( fp );
				return true;
			}
		}
	}
	else {
		fclose (fp);
		return false;
	}

	strcpy(arch,"INTEL");
	fclose (fp);
	return true;

}



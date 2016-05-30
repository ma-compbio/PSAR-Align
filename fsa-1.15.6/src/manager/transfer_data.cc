/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include "fsa/fsa.h"
#include <errno.h>

#if defined(HAVE_CONDOR) 
#include "MW.h"
#include "MWSystem.h"
#endif

#define END_DATA    0
#define SPARSE_DATA 1
#define HEAP_DATA   2

using namespace fsa;

FSA *m_fsa = NULL;
DB_adapter *m_db_adapter = NULL;

Params m_params_seed;
Params m_pseudocounts;

MEM_Buffer m_mem_buffer;
DB_Buffer m_db_buffer;

static int num_probs = 0;
static int num_logs = 0;
static int num_edges = 0;

typedef struct MEM_Container {
	int what_data;

	int size1;
	char* buffer1;
	int size2;
	char* buffer2;

} MEM_Container;

static std::queue<MEM_Container*> buffer_queue;

bool final_buffer = false;
bool send_thread_waiting = false;

static pthread_mutex_t queue_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t queue_cond = PTHREAD_COND_INITIALIZER;

static double wall_time = 0.0;
static double cpu_time = 0.0;

void compute_data (SparseMatrix *ijMatrix, int i, int j); 

void insert_posterior_probability (int seq1, int pos1, int seq2, int pos2, float prob);
void insert_num_cells (int seq1, int seq2, int size);
void insert_edge (int seq1, int pos1, int seq2, int pos2, double weight, float delta);

void flush_on_the_fly (bool force);

void build_multiple_alignment();
void build_anchored_multiple_alignment();

static void enqueue_buffer(MEM_Container* buffer, bool final);
static MEM_Container* dequeue_buffer(void);

#if defined(HAVE_POSTGRES) 
void init_single_worker(const FSA *fsa, DB_adapter *db_adapter, Params &params_seed, Params &pseudocounts); 
#endif

#if defined(HAVE_CONDOR) 
static void* send_thread(void *ptr);
static bool do_raw_pack (MWRMComm* RMC, MEM_Container *mem_container);
static bool raw_pack(MWRMComm* _RMC, char* buf, int len);
void transfer_data(MW_task *t);
#endif


void insert_posterior_probability (int seq1, int pos1, int seq2, int pos2, float prob) {

	if (m_fsa->write_db) {
		char entry[MAX_DB_ENTRY_LENGTH];
		sprintf(entry,"%d\t%d\t%d\t%d\t%lf\n", seq1, pos1, seq2, pos2, prob);
		m_db_buffer.sparse_matrix.append (entry);
	}
	else {
		m_mem_buffer.sparse_matrix[num_probs].pos1 = pos1;
		m_mem_buffer.sparse_matrix[num_probs].pos2 = pos2;
		m_mem_buffer.sparse_matrix[num_probs].prob = prob;
	}

	num_probs++;
}

void insert_num_cells (int seq1, int seq2, int size) {


	if (m_fsa->write_db) {
		char entry[MAX_DB_ENTRY_LENGTH];
		sprintf(entry,"%d\t%d\t%d\n", seq1, seq2, size);
		m_db_buffer.num_cells.append (entry);
	}
	else {
		if (seq1 == 0 && seq2 == 0 )
			return;

		m_mem_buffer.num_cells[num_logs].seq1 = seq1;
		m_mem_buffer.num_cells[num_logs].seq2 = seq2;
		m_mem_buffer.num_cells[num_logs].size = size;
	}

	num_logs++;
}

void insert_edge (int seq1, int pos1, int seq2, int pos2, double weight, float delta) {


	if (m_fsa->write_db) {
		char entry[MAX_DB_ENTRY_LENGTH];
		sprintf(entry,"%d\t%d\t%d\t%d\t%lf\t%f\n", seq1, pos1, seq2, pos2, weight, delta);
		m_db_buffer.heap.append (entry);
	}
	else {

		m_mem_buffer.heap[num_edges].seq1 = seq1;
		m_mem_buffer.heap[num_edges].pos1 = pos1;
		m_mem_buffer.heap[num_edges].seq2 = seq2;
		m_mem_buffer.heap[num_edges].pos2 = pos2;
		m_mem_buffer.heap[num_edges].weight = weight;
		m_mem_buffer.heap[num_edges].delta = delta;
	}

	num_edges++;
}


void compute_data (SparseMatrix *ijMatrix, int i, int j) 
{
	if (ijMatrix == 0)
		return;

	int xlen = m_fsa->seq_db.get_seq (i).length();
	int ylen = m_fsa->seq_db.get_seq (j).length();

	int pre_num_probs = num_probs;

	for (int jj = 0; jj < ylen; jj++)  
		insert_posterior_probability (i, 0, j, jj+1, ijMatrix->get_gap_prob (1, jj+1));

	// for all entries in the SparseMatrix, calculate the corresponding weight
	// and store the edge
	for (int ii = 0; ii < xlen; ii++) { // note 0-based indexing for sequences
		float p_gap_ii = ijMatrix->get_gap_prob (0, ii + 1);

		insert_posterior_probability (i, ii+1, j, 0, p_gap_ii);

		for (std::vector<Matrix_entry>::iterator rowPtr = ijMatrix->GetRowPtr (ii + 1),
				rowEnd = rowPtr + ijMatrix->GetRowSize (ii + 1); rowPtr != rowEnd; rowPtr++) {
			int jj = rowPtr->first - 1; // convert from SparseMatrix's 1-based coords to the 0-based coords which we use here

			insert_posterior_probability (i, ii+1, j, jj+1, ijMatrix->get_match_prob(ii+1, jj+1));

			float p_match = rowPtr->second;
			if (!p_match)
				continue;

			float p_gap = p_gap_ii + ijMatrix->get_gap_prob (1, jj + 1);
			float weight = m_fsa->use_tgf ? (2 * p_match / p_gap) : (2 * p_match - m_fsa->gap_factor * p_gap);
			float delta = 2 * p_match - p_gap;

			float edge_weight_threshold = 0; //TODO it should be the same with in fsa

			// if the edge already scores too low to ever be accepted, then skip it
			if ((weight < edge_weight_threshold) || (m_fsa->use_tgf && weight < m_fsa->gap_factor))
				continue;

			insert_edge (i, ii, j, jj, weight, delta);
		}
	}

	insert_num_cells (i, j, (num_probs - pre_num_probs));
}

void build_anchored_multiple_alignment() {

  // read in weights for sequence pairs if information present
  if (m_fsa->tree_weights_file != "")
    m_fsa->tree_weights.from_file (m_fsa->seq_db, m_fsa->tree_weights_file);

    // get resolved anchors for the sequence pairs which we're going use for alignment
    Anchor_resolver anchor_resolver (m_fsa->seq_db, m_fsa->seq_db_internal, m_fsa->alignment_seq_pairs);

    // add Mercator constraints if present
    if (m_fsa->mercator_constraint_file != "")
      anchor_resolver.add_mercator_constraints (m_fsa->mercator_constraint_file);

    // now actually get resolved anchors
    std::vector<Anchors > resolved_anchors_list = anchor_resolver.get_resolved_anchors (m_fsa->tree_weights,
											m_fsa->anchor_minlen, m_fsa->anchor_max_join_length, m_fsa->use_translated,
										   m_fsa->exonerate, m_fsa->exonerate_minscore, m_fsa->softmasked,
										   m_fsa->hardmasked,
										   m_fsa->num_refinement_steps,
										   m_fsa->output_for_gui, m_fsa->gui_prefix);

    // loop through sequence database
    for (size_t cnt = 0; cnt < m_fsa->alignment_seq_pairs.size(); ++cnt) {

      // initialize all the sequence data
      int i = m_fsa->alignment_seq_pairs[cnt].first;
      int j = m_fsa->alignment_seq_pairs[cnt].second;

	  // initialize all the sequence data
      const Sequence& xseq = m_fsa->seq_db_internal.get_seq (i);
      const Sequence& yseq = m_fsa->seq_db_internal.get_seq (j);

      Post_probs post_probs = m_fsa->perform_anchored_pairwise_inference (m_params_seed, xseq, yseq,
										   m_pseudocounts,
										   resolved_anchors_list[cnt]);

      // we need the original sequences in case we've been hardmasking
      const Sequence& xseq_orig = (m_fsa->seq_db).get_seq (i);
      const Sequence& yseq_orig = (m_fsa->seq_db).get_seq (j);

      // if hardmasking, map coords back to original sequence
      if (m_fsa->hardmasked) {
	for (Post_probs::iterator p = post_probs.begin(); p != post_probs.end(); ++p) {
	  (*p).x = xseq_orig.map_stripped_to_orig ((*p).x);
	  (*p).y = yseq_orig.map_stripped_to_orig ((*p).y);
	}
      }

      SparseMatrix *ijMatrix = new SparseMatrix (i, j,
						 xseq_orig.length(), yseq_orig.length(),
						 post_probs);

	  compute_data (ijMatrix, i, j);

	  flush_on_the_fly (false);

	  delete ijMatrix;

    } // end loop over seq_db_internal

	insert_num_cells (0, 0, num_edges);

	flush_on_the_fly (true);

}


void build_multiple_alignment() {

  // read in weights for sequence pairs if information present
  if (m_fsa->tree_weights_file != "")
    m_fsa->tree_weights.from_file (m_fsa->seq_db, m_fsa->tree_weights_file);

	bool left_match = false;  // do not require homology beyond sequence boundaries
	bool right_match = false;

	// loop through sequence database
	for (size_t cnt = 0; cnt < m_fsa->alignment_seq_pairs.size(); ++cnt) {

		// initialize all the sequence data
		int i = m_fsa->alignment_seq_pairs[cnt].first;
		int j = m_fsa->alignment_seq_pairs[cnt].second;

		// initialize all the sequence data
		const Sequence& xseq = m_fsa->seq_db_internal.get_seq (i);
		const Sequence& yseq = m_fsa->seq_db_internal.get_seq (j);
		
		Params params_seed;
		params_seed.copy_all (m_params_seed);

		Post_probs post_probs = m_fsa->perform_pairwise_inference  (params_seed, xseq, yseq, 
										   left_match, right_match, m_fsa->ragged_ends, m_pseudocounts);

      // we need the original sequences in case we've been hardmasking
      const Sequence& xseq_orig = (m_fsa->seq_db).get_seq (i);
      const Sequence& yseq_orig = (m_fsa->seq_db).get_seq (j);

      // if hardmasking, map coords back to original sequence
      if (m_fsa->hardmasked) {
	for (Post_probs::iterator p = post_probs.begin(); p != post_probs.end(); ++p) {
	  (*p).x = xseq_orig.map_stripped_to_orig ((*p).x);
	  (*p).y = yseq_orig.map_stripped_to_orig ((*p).y);
	}
      }

      SparseMatrix *ijMatrix = new SparseMatrix (i, j,
						 xseq_orig.length(), yseq_orig.length(),
						 post_probs);

		compute_data(ijMatrix, i, j);

		delete ijMatrix;

		flush_on_the_fly (false);
	}

	insert_num_cells (0, 0, num_edges);

	flush_on_the_fly (true);
	
}

void flush_on_the_fly (bool force) {
	if (m_fsa->write_db) {
		if ( force || ( (m_db_buffer.sparse_matrix.size() + m_db_buffer.num_cells.size() + m_db_buffer.heap.size() ) > (BUFFER_RATIO * MAX_DB_BUFFER_SIZE))) {
#if defined(HAVE_POSTGRES) 
			(m_db_adapter->get_connection())->flush_sparse_matrix_buffer (m_db_buffer.sparse_matrix.c_str());
			m_db_buffer.sparse_matrix.clear();

			(m_db_adapter->get_connection())->flush_num_cells_buffer (m_db_buffer.num_cells.c_str());
			m_db_buffer.num_cells.clear();

			(m_db_adapter->get_connection())->flush_heap_buffer (m_db_buffer.heap.c_str());
			m_db_buffer.heap.clear();
#endif
		}

	} else  {
		
		if ( force || (num_probs * sizeof (Sparse_Matrix_Buffer) > (BUFFER_RATIO * MAX_MEM_BUFFER_SIZE))) {
			MEM_Container *mem_container = (MEM_Container*) malloc ( sizeof (MEM_Container) );
			assert (mem_container);

			mem_container->what_data = SPARSE_DATA;

			mem_container->size1 = num_logs * sizeof (Num_Cells_Buffer);
			mem_container->buffer1 = (char*) m_mem_buffer.num_cells;

			mem_container->size2 = num_probs * sizeof (Sparse_Matrix_Buffer);
			mem_container->buffer2 = (char*) m_mem_buffer.sparse_matrix;

			enqueue_buffer(mem_container, false);

			if (!force) {

				m_mem_buffer.num_cells = (Num_Cells_Buffer *) malloc ( 
						sizeof (Num_Cells_Buffer) * ((int) (MAX_DB_BUFFER_SIZE / sizeof (Num_Cells_Buffer))));
				assert(m_mem_buffer.num_cells);

				m_mem_buffer.sparse_matrix = (Sparse_Matrix_Buffer*) malloc (
						sizeof (Sparse_Matrix_Buffer) * ((int) (MAX_DB_BUFFER_SIZE / sizeof (Sparse_Matrix_Buffer))));
				assert(m_mem_buffer.sparse_matrix);
			} else {
				m_mem_buffer.num_cells = NULL;
				m_mem_buffer.sparse_matrix = NULL;
			}

			num_probs = 0;
			num_logs = 0;
		}

		if ( force || (num_edges * sizeof (Heap_Buffer) > (BUFFER_RATIO * MAX_MEM_BUFFER_SIZE))) {

			MEM_Container *mem_container = (MEM_Container*) malloc ( sizeof (MEM_Container) );
			assert (mem_container);

			mem_container->what_data = HEAP_DATA;

			mem_container->size1 = num_edges * sizeof (Heap_Buffer);
			mem_container->buffer1 = (char*) m_mem_buffer.heap;

			mem_container->size2 = 0;
			mem_container->buffer2 = NULL;

			enqueue_buffer(mem_container, false);

			if (!force) {
				m_mem_buffer.heap = (Heap_Buffer *) malloc ( 
						sizeof (Heap_Buffer) * ((int) (MAX_DB_BUFFER_SIZE / sizeof (Heap_Buffer))));
				assert(m_mem_buffer.heap);
			}
			else
				m_mem_buffer.heap = NULL;

			num_edges = 0;
		}

		if ( force ) {
			// send last signal
			MEM_Container *mem_container = (MEM_Container*) malloc ( sizeof (MEM_Container) );
			assert (mem_container);

			mem_container->what_data = END_DATA;
			mem_container->size1 = 0;
			mem_container->buffer1 = NULL;

			mem_container->size2 = 0;
			mem_container->buffer2 = NULL;
			
			enqueue_buffer(mem_container, true);
		}
	}
}



static void enqueue_buffer(MEM_Container* buffer, bool final)
{
	if( !buffer ) {
		return;
	}

	pthread_mutex_lock(&queue_lock);

	buffer_queue.push(buffer);
	if( final ) {
		final_buffer = true;
	}

	if( send_thread_waiting ) {
		pthread_cond_signal(&queue_cond);
	}

	pthread_mutex_unlock(&queue_lock);

}

static MEM_Container* dequeue_buffer (void) {

	MEM_Container *mem_container;

	if (buffer_queue.empty() )
		mem_container = NULL;
	else {
		mem_container = buffer_queue.front();
		buffer_queue.pop();
	}

	return  mem_container;
}

#if defined(HAVE_CONDOR) 

static void* send_thread(void *ptr) 
{
	
	MWTask* t = (MWTask*)ptr;

	while(1) {
		MEM_Container* mem_container = NULL;
		bool done = false;

		pthread_mutex_lock(&queue_lock);
		while( (mem_container = dequeue_buffer()) == NULL ) {

			if( final_buffer ) { 
				break;
			}else {
				send_thread_waiting = true;
				pthread_cond_wait(&queue_cond, &queue_lock);
				send_thread_waiting = false;
			}	
		}

		if(buffer_queue.empty()) {
			done = final_buffer;
		}

		pthread_mutex_unlock(&queue_lock);

		if( mem_container ) {
			if(!do_raw_pack (t->RMC, mem_container) ) {
				fprintf(stderr, "Failed to send data\n");
				exit(1);
			}

			if (mem_container->buffer1) {
				free (mem_container->buffer1);
			}

			if (mem_container->buffer2) {
				free (mem_container->buffer2);
			}
			

			free(mem_container);
		}

		if( done ) {
			// Now, we have sent all data
			// Here, we will send two extra double for stats
			wall_time += MWSystem::gettimeofday();
			cpu_time += MWSystem::getcputime();
			if (!raw_pack (t->RMC, (char*)&wall_time, sizeof(double))) {
				fprintf(stderr, "Failed to send stats data\n");
				exit(1);
			}
			if (!raw_pack (t->RMC, (char*)&cpu_time, sizeof(double))) {
				fprintf(stderr, "Failed to send stats data\n");
				exit(1);
			}
			break;
		}
	}

	return NULL;
}

static bool raw_pack(MWRMComm* _RMC, char* buf, int len)
{
	static int cnt = 0;

	cnt++;
	int ret = _RMC->raw_pack(buf, len);
	if( ret != len ) {
		fprintf(stderr, "raw_pack(%d) failed to write the size of data.\n", cnt);
		fprintf(stderr, "need to send (%d), but ret = %d\n", len, ret);
		if( ret < 0 ) {
			fprintf(stderr, "errno=%d, string=%s\n", errno, strerror(errno));
			fflush(stderr);
		}
		return false;
	}

	return true;
}

static bool do_raw_pack (MWRMComm* RMC, MEM_Container *mem_container) 
{
	int network_int = mem_container->what_data;
	if (!raw_pack (RMC, (char*) &network_int, sizeof (int)))
		return false;

	switch (mem_container->what_data) {
		case SPARSE_DATA : 
			{
				network_int = ((int) (mem_container->size1 / sizeof(Num_Cells_Buffer)));

				if (!raw_pack (RMC, (char*) &network_int, sizeof (int)))
					return false;

				if (!raw_pack (RMC, (char*) mem_container->buffer1, mem_container->size1))
					return false;

				network_int = ((int) mem_container->size2 / sizeof(Sparse_Matrix_Buffer));

				if (!raw_pack (RMC, (char*) &network_int, sizeof (int)))
					return false;

				if (!raw_pack (RMC, (char*) mem_container->buffer2, mem_container->size2))
					return false;

				break;
			}
		case HEAP_DATA : 
			{
				network_int = ((int) mem_container->size1 / sizeof (Heap_Buffer));
				if (!raw_pack (RMC, (char*) &network_int, sizeof(int))) {
					return false;
				}

				if (!raw_pack (RMC, (char*) mem_container->buffer1, mem_container->size1)) {
					return false;
				}

				break;
			}
		default :
			break;
	}

	return true;
}

void transfer_data(MW_task *t)
{

	// Record times for stats
	// Set our stopwatch
	wall_time -= MWSystem::gettimeofday();
	cpu_time -= MWSystem::getcputime();

	pthread_t thread;
	if( t->m_raw_flag ) {
		// we need to create a new thread
		if( pthread_create(&thread, NULL, send_thread, (void *)t ) != 0 ) {
			fprintf(stderr, "Failed to create a new thread\n");
			exit(1);
		}
	}

	if (m_fsa->anchored)
		build_anchored_multiple_alignment ();
	else
		build_multiple_alignment ();

	delete m_fsa;

	if( t->m_raw_flag ) {
		// Now, we have to wait for send_thread to be finished.
		void *ptr;
		int ret = pthread_join(thread, &ptr);
		if( ret != 0 ) {
			fprintf(stderr, "pthread_join error %d\n", ret);
			exit(1);
		}
		
	}
}

#endif

#if defined(HAVE_POSTGRES) 
void init_single_worker(const FSA *fsa, DB_adapter *db_adapter, Params &params_seed, Params &pseudocounts) {
	m_fsa = const_cast<FSA *> (fsa);
	m_db_adapter = db_adapter;

	m_params_seed.copy_all (params_seed);
	m_pseudocounts.copy_all (pseudocounts);

	int seqs_schema_id = 0;
	int params_table_id = 0;

	final_buffer = false;
	send_thread_waiting = false;

	bool is_db_running = false;

	m_db_buffer.sparse_matrix.clear();
	m_db_buffer.num_cells.clear();
	m_db_buffer.heap.clear();

	// set the sequence schema id and the parameter table id 
	seqs_schema_id  = (m_db_adapter->get_connection())->get_seqs_schema_id ();
	params_table_id = (m_db_adapter->get_connection())->get_params_table_id ();

	m_db_adapter->set_up_ids (seqs_schema_id, params_table_id, 0);
}
#endif


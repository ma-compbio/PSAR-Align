/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include "manager/manager.h"
#include "annealing/alignment_DAG.h"
#include "fsa/fsa.h"

using namespace fsa;


#ifdef HAVE_POSTGRES
extern void init_single_worker(const FSA *fsa, DB_adapter *db_adapter, Params &params_seed, Params &pseudocounts); 
extern void build_multiple_alignment();
extern void build_anchored_multiple_alignment();
#endif

Manager::Manager() {
	m_seq_db_internal = NULL;
	m_db_opts = NULL;

#ifdef HAVE_POSTGRES
	m_num_cells = NULL;
	m_available_sparse_matrices = NULL;
	m_sparse_matrix_scheduler = NULL;
#endif

}

Manager::Manager(const Sequence_database &seq_db_internal, DB_opts &db_opts) {

	m_seq_db_internal = &seq_db_internal;
	m_db_opts = &db_opts;

#ifdef HAVE_POSTGRES
	bool is_db_running, is_db_available;

	// make a connection to the database
	is_db_running = m_db_adapter.connect_db(m_db_opts->db_hostname.c_str(), m_db_opts->db_hostaddr.c_str(), 
			m_db_opts->db_name.c_str(), m_db_opts->db_port, m_db_opts->db_user.c_str(), m_db_opts->db_password.c_str());

	if (!is_db_running && m_db_opts->write_db)
		THROWEXPR ("ERROR: We can not access database server. exit");

	// try to find data of input sequences
	is_db_available = m_db_adapter.look_up_data (seq_db_internal, db_opts); 

	if (is_db_running && !m_db_opts->write_db && !is_db_available) 
		CL << "WARNING: there's no relevent data of the input sequences in the database." << endl;

	// initialize variables
	SPARSE_MATRIX_MAX_SIZE = 0;
	HEAP_WINDOW_SIZE = 0;

	m_db_merged_heap_buffer.clear();

	m_num_inserted_edges = 0;
	m_min_edge_weight = 0;
	m_merged_heap_offset = 0;
	m_num_cells = NULL;
	m_available_sparse_matrices = NULL;
	m_sparse_matrix_scheduler = NULL;

	m_last_edge_chunk =  (m_db_opts->db_max_ram > 0)? false : true;
#endif

}

Manager::~Manager() {
#ifdef HAVE_POSTGRES
	if (m_num_cells)
		delete m_num_cells;

	if (m_available_sparse_matrices)
		delete m_available_sparse_matrices;

	if (m_sparse_matrix_scheduler)
		delete m_sparse_matrix_scheduler;
#endif
}

bool Manager::check_available_sparse_matrices () { 

#ifdef HAVE_POSTGRES
	DB_postgres *db_connection = m_db_adapter.get_connection ();

	if (db_connection) {

		int num_jobs = db_connection->get_num_jobs ();
		int num_seqs = m_seq_db_internal->size();

		if (num_jobs > 0) {	
			int seq_pair_position = 1;
			int num_seq_pairs = num_seqs * (num_seqs - 1 ) / 2; 
			int num_of_pairs = num_seq_pairs / num_jobs; 
			int num_remains = num_seq_pairs % num_jobs;

			for (int i=0; i<num_jobs; i++) {
				m_sparse_matrix_ref.push_back (seq_pair_position);
				int length = (i < num_remains) ? num_of_pairs + 1 : num_of_pairs;
				seq_pair_position += length;
			}	
		}

		CTAG(4,MANAGER) << "Checking available sparse matrices." << endl;
		m_available_sparse_matrices = new std::vector<std::vector<int> > (num_seqs, std::vector<int> (num_seqs, -1)); 

		std::vector<int> added_workers;

		avg_sparse_matrix_size = 0.0;
		int num_pairs = 0;

		for (int i = 0; i < num_jobs; i++)
			added_workers.push_back (i);

		while (!added_workers.empty()) {
			int id = added_workers.front();
			if (db_connection->get_list_of_available_pairs (id, *m_available_sparse_matrices, m_seq_db_internal, avg_sparse_matrix_size, num_pairs, m_orig_edges_size))
				added_workers.erase (added_workers.begin());
		}


		if (m_db_opts && m_db_opts->db_max_ram > 0) {
		  m_sparse_matrix_scheduler = new std::queue <Sparse_matrix_entry> ();
			m_cnt_sparse_matrices = 0;
			
			// by the default, we assign 90% of total memory for the sparse matirces and the 10% of it to the priority queue
			SPARSE_MATRIX_MAX_SIZE = (int)(0.90 * m_db_opts->db_max_ram * 1024 * 1024) / (sizeof (SparseMatrix) + (int) avg_sparse_matrix_size);
			HEAP_WINDOW_SIZE = (int)(0.10 * m_db_opts->db_max_ram * 1024 * 1024 ) / sizeof (Edge);

			if (HEAP_WINDOW_SIZE > m_orig_edges_size) {
				// in this case, we can assign more memory for the sparse matrices
				HEAP_WINDOW_SIZE = m_orig_edges_size; 
				SPARSE_MATRIX_MAX_SIZE = (int)(m_db_opts->db_max_ram * 1024 * 1024 - sizeof(Edge) * m_orig_edges_size) / (sizeof (SparseMatrix) + (int) avg_sparse_matrix_size);
			} 
		} 
	}

	return true;
#endif

	return false;
}


#ifdef HAVE_POSTGRES
bool Manager::update_size (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges) {
	if (!m_db_opts)
		return false;

	SPARSE_MATRIX_MAX_SIZE = (int)(m_db_opts->db_max_ram * 1024 * 1024 - sizeof(Edge) * edges.size()) / (sizeof (SparseMatrix) + (int) avg_sparse_matrix_size);

	return true;
}
#endif

bool Manager::get_num_cells () { 

#ifdef HAVE_POSTGRES
	DB_postgres *db_connection = m_db_adapter.get_connection ();

	if (db_connection) {

		int num_jobs = db_connection->get_num_jobs ();
		int num_seqs = m_seq_db_internal->size();

		if (num_jobs > 0) {	
			int seq_pair_position = 1;
			int num_seq_pairs = num_seqs * (num_seqs - 1 ) / 2; 
			int num_of_pairs = num_seq_pairs / num_jobs; 
			int num_remains = num_seq_pairs % num_jobs;

			for (int i=0; i<num_jobs; i++) {
				m_sparse_matrix_ref.push_back (seq_pair_position);
				int length = (i < num_remains) ? num_of_pairs + 1 : num_of_pairs;
				seq_pair_position += length;
			}	
		}

		CTAG(4,MANAGER) << "Getting the number of cells of each sparse matrix." << endl;
		m_num_cells = new std::vector<std::vector<int> > (num_seqs, std::vector<int> (num_seqs, 0)); 

		std::vector<int> added_workers;

		for (int i = 0; i < num_jobs; i++)
			added_workers.push_back (i);

		while (!added_workers.empty()) {
			int id = added_workers.front();
			if (db_connection->get_num_cells (id, *m_num_cells) == DB_OK)
				added_workers.erase (added_workers.begin());
		}


		m_orig_edges_size = (*m_num_cells)[0][0]; /// [0][0] => heap size 
	}
	return true;
#endif

	return false;
}

bool Manager::is_edges_available () { 
#ifdef HAVE_CONDOR
	if (!(m_mem_buffers.heap_buffers).empty()) {
		return true;
	}
#endif

#ifdef HAVE_POSTGRES
	if (m_db_adapter.is_data_available()) {
		return true;
	}
#endif
	return false;
}

bool Manager::is_sparse_matrices_available () { 
#ifdef HAVE_CONDOR
	if (!(m_mem_buffers.sparse_matrix_buffers).empty()) {
		return true;
	}
#endif
#ifdef HAVE_POSTGRES
	if (m_db_adapter.is_data_available()) {
		return true;
	}
#endif

	return false;
}

bool Manager::is_sparse_matrix_available (const int i, const int j) {


#ifdef HAVE_POSTGRES
	if ( (m_num_cells && ((*m_num_cells)[i][j] > 0 || (*m_num_cells)[j][i]>0)) || 
			(m_available_sparse_matrices && ((*m_available_sparse_matrices)[i][j] > -1|| (*m_available_sparse_matrices)[j][i] > -1)) )
		return true;
#endif
	return false;

}

#ifdef HAVE_POSTGRES
int Manager::look_up_sparse_matrix_table_id (const int i, const int j) { 

	// binary search
	int num_seqs = m_seq_db_internal->size();
	int input_seq_position = (int) (i * (2 * num_seqs - i - 1) / 2) + (j - i);

	int first = 0;
	int last = m_sparse_matrix_ref.size() - 1;

	//binary search
	while (first <= last) {
		int mid = (first + last) / 2; //compute mid point

		if ( m_sparse_matrix_ref[mid] <= input_seq_position ) {
			if (((mid == m_sparse_matrix_ref.size()-1)) || (input_seq_position < m_sparse_matrix_ref[mid+1])) 
				return mid;
			else
				first = mid + 1;
		}

		if ( input_seq_position < m_sparse_matrix_ref[mid] ) {

			if (mid == 0 || (m_sparse_matrix_ref[mid-1] <= input_seq_position))  
				return mid-1;
			else
				last = mid - 1;
		}
	}

	return -1;
}
#endif

bool Manager::get_sparse_matrix (std::vector<std::vector<SparseMatrix*> >& sparse_matrices, const int i, const int j) { 

#ifdef HAVE_POSTGRES
	// get a sparse matrix from the database
	if ( !is_sparse_matrix_available (i, j) ) 
		return false;

	if (sparse_matrices[i][j] != NULL)
		return true;

	assert (sparse_matrices[i][j] == NULL);

	// see if transpose is possible
	if (sparse_matrices[j][i]) {
		sparse_matrices[i][j] = sparse_matrices[j][i]->ComputeTranspose();

		//check 
		if (m_cnt_sparse_matrices > SPARSE_MATRIX_MAX_SIZE) {
			SparseMatrix **sparse_matrix_ptr = m_sparse_matrix_scheduler->front();
			m_sparse_matrix_scheduler->pop();
			delete (*sparse_matrix_ptr);
			(*sparse_matrix_ptr) = NULL;
			m_cnt_sparse_matrices--;
		}

		if (m_sparse_matrix_scheduler->size() < MAX_SPARSE_MATRIX_COUNT) {
			m_sparse_matrix_scheduler->push (&sparse_matrices[i][j]);
			m_cnt_sparse_matrices++;
		} else
			m_cnt_sparse_matrices++;

		return true;
	}

	DB_postgres *db_connection = m_db_adapter.get_connection ();

	if (!db_connection)
		return false;

	int seq1 = (i<j)? i : j;
	int seq2 = (i<j)? j : i;

	int seq1Length = m_seq_db_internal->get_seq (seq1).length();
	int seq2Length = m_seq_db_internal->get_seq (seq2).length();

	/*get the table for (seq1, seq2) from the database */
	int id = (*m_available_sparse_matrices)[i][j];
	if (id == -1) 
		return false;

	bool res = db_connection->get_sparse_matrix (id, seq1, seq2, seq1Length, seq2Length, sparse_matrices);

	assert (sparse_matrices[i][j]);

	//check 
	if (m_cnt_sparse_matrices > SPARSE_MATRIX_MAX_SIZE) {
		SparseMatrix **sparse_matrix_ptr = m_sparse_matrix_scheduler->front();
		m_sparse_matrix_scheduler->pop();
		delete (*sparse_matrix_ptr);
		(*sparse_matrix_ptr) = NULL;
		m_cnt_sparse_matrices--;
	}

	if (m_sparse_matrix_scheduler->size() < MAX_SPARSE_MATRIX_COUNT) {
		m_sparse_matrix_scheduler->push (&sparse_matrices[i][j]);
		m_cnt_sparse_matrices++;
	} else
		m_cnt_sparse_matrices++;

	if (res && sparse_matrices[i][j] == NULL) {
		sparse_matrices[i][j] = sparse_matrices[j][i]->ComputeTranspose();

		//check 
		if (m_cnt_sparse_matrices > SPARSE_MATRIX_MAX_SIZE) {
			SparseMatrix **sparse_matrix_ptr = m_sparse_matrix_scheduler->front();
			m_sparse_matrix_scheduler->pop();
			delete (*sparse_matrix_ptr);
			(*sparse_matrix_ptr) = NULL;
			m_cnt_sparse_matrices--;
		}

		if (m_sparse_matrix_scheduler->size() < MAX_SPARSE_MATRIX_COUNT) {
			m_sparse_matrix_scheduler->push (&sparse_matrices[i][j]);
			m_cnt_sparse_matrices++;
		} else
			m_cnt_sparse_matrices++;
	}

	assert (sparse_matrices[i][j]);

	return res;
#else
	return false;
#endif
}

bool Manager::get_all_sparse_matrices (std::vector<std::vector<SparseMatrix*> >& sparse_matrices) {

#ifdef HAVE_CONDOR
	// get pairwise posterior probabilties from the database
	Cells_vector &num_cells_buffers = m_mem_buffers.num_cells_buffers;
	Sparse_vector &sparse_matrix_buffers = m_mem_buffers.sparse_matrix_buffers;

	if (!num_cells_buffers.empty()) {

		Num_Cells_Buffer *num_cells_buffer = NULL;
		Sparse_Matrix_Buffer *sparse_matrix_buffer = NULL;
		int size = 0;

		while (!num_cells_buffers.empty()) {
			Cells_pair &cells_pair = num_cells_buffers.front();
			Sparse_pair &sparse_pair = sparse_matrix_buffers.front();

			size = cells_pair.first;

			Num_Cells_Buffer *num_cells_buffer = cells_pair.second;
			Sparse_Matrix_Buffer *sparse_matrix_buffer = sparse_pair.second;

			int pre_len = 0;
			int len = 0;

			// create a sparse matrix and the trasposed one
			for (int cnt = 0; cnt < size; cnt++) {
				int i = num_cells_buffer[cnt].seq1;
				int j = num_cells_buffer[cnt].seq2;

				len = num_cells_buffer[cnt].size;

				int seq1Length = m_seq_db_internal->get_seq (i).length();
				int seq2Length = m_seq_db_internal->get_seq (j).length();

				sparse_matrices[i][j] = new SparseMatrix (i, j,
									  seq1Length, seq2Length,
									  pre_len, len, sparse_matrix_buffer);
				sparse_matrices[j][i] = sparse_matrices[i][j]->ComputeTranspose();

				pre_len += len;
			}

			assert (pre_len == sparse_pair.first);
			if (num_cells_buffer)
				free (num_cells_buffer); // prevent memory leak

			if (sparse_matrix_buffer)
				free (sparse_matrix_buffer);
			num_cells_buffers.erase (num_cells_buffers.begin()); //erase
			sparse_matrix_buffers.erase (sparse_matrix_buffers.begin());
		}
		assert (m_mem_buffers.num_cells_buffers.size() == 0);
		assert (m_mem_buffers.sparse_matrix_buffers.size() == 0);

		return true;
	}

#endif

#ifdef HAVE_POSTGRES
	// get pairwise posterior probabilities from the workers
	DB_postgres *db_connection = m_db_adapter.get_connection ();

	if (!db_connection)
		return false;

	int num_seqs = m_seq_db_internal->size();

	// get the nubmer of jobs
	int num_jobs = db_connection->get_num_jobs ();

	if ( num_jobs == 0 )
		return false;

	std::vector<int> added_workers;

	for (int i = 0; i < num_jobs; i++)
		added_workers.push_back (i);

	while (!added_workers.empty()) {
		int id = added_workers.front();
		if (db_connection->get_sparse_matrices (id, *m_seq_db_internal, *m_num_cells, sparse_matrices) == DB_OK) 
			added_workers.erase (added_workers.begin());
	}

	return true;
#endif
	return false;
}


bool Manager::get_all_edges (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
		std::vector<Seq_pos_col_map> &seq_pos_col_maps) {  

#ifdef HAVE_CONDOR
	// get all edges from the workers
	Heap_vector &heap_buffers = m_mem_buffers.heap_buffers;

	if (!heap_buffers.empty()) {

		while (!heap_buffers.empty()) {
			Heap_pair &heap_pair = heap_buffers.front();

			int size = heap_pair.first;
			Heap_Buffer *heap_buffer = heap_pair.second;
			
			// create a new edge and push it into the priority queue
			for (int i = 0; i < size; i++) {
				Edge *edge = new Edge(seq_pos_col_maps[heap_buffer[i].seq1][heap_buffer[i].pos1], 
						seq_pos_col_maps[heap_buffer[i].seq2][heap_buffer[i].pos2], 
						      std::pair<int, int> (heap_buffer[i].seq1, heap_buffer[i].pos1), 
						      std::pair<int, int> (heap_buffer[i].seq2, heap_buffer[i].pos2), 
						heap_buffer[i].weight, heap_buffer[i].delta, 2);

				edges.push(edge);

			}
			if (heap_buffer) 
				free (heap_buffer); // prevent memory leak
			heap_buffers.erase (heap_buffers.begin()); //erase

		}
		return true;
	}
#endif

#ifdef HAVE_POSTGRES
	// get all edges from the database
	DB_postgres *db_connection = m_db_adapter.get_connection ();

	if (!db_connection)
		return false;

	int num_seqs = m_seq_db_internal->size();

	// get the number of jobs ( = workers)
	int num_jobs = db_connection->get_num_jobs ();

	if ( num_jobs == 0 )
		return false;

	std::vector<int> added_workers;

	for (int i = 0; i < num_jobs; i++)
		added_workers.push_back (i);

	while (!added_workers.empty()) {
		int id = added_workers.front();
		if (db_connection->get_heaps (id, m_min_edge_weight, seq_pos_col_maps, edges) == DB_OK)
			added_workers.erase (added_workers.begin());
	}

	return true;
#endif

	return false;
}


bool Manager::mw_master_run (int argc, char** argv, const Params &params_seed, const Params &pseudocounts) {

#ifdef HAVE_CONDOR
	if (m_db_opts && m_db_opts->write_db)  {

	int seqs_schema_id = 0;
	int params_table_id = 0;

#ifdef HAVE_POSTGRES
		if ( !m_db_adapter.init_database (*m_seq_db_internal, *m_db_opts) ) 
			THROWEXPR ("ERROR: We can not initialize database server. exit");

		// set the sequence schema id and the parameter table id 
		seqs_schema_id  = (m_db_adapter.get_connection())->get_seqs_schema_id ();
		params_table_id = (m_db_adapter.get_connection())->get_params_table_id ();
#endif 
	
	}
	// if it is not the database mode , just send dummy values for seqs_schema_id and params_table_id
	m_mw_adapter.master_run (argc, argv, params_seed, pseudocounts, m_mem_buffers, 
			m_seq_db_internal->size(), m_db_opts->num_alignment_pairs, m_db_opts->num_parallelized_jobs,
			seqs_schema_id, params_table_id);

#ifdef HAVE_POSTGRES
	m_db_adapter.look_up_data(*m_seq_db_internal, *m_db_opts);
#endif //endif for HAVE_POSTGRES
	return true;

#else
	// else for HAVE_CONDOR
	return false;
#endif
}

bool Manager::mw_worker_run(int argc, char** argv) {
	// run the worker instance

#ifdef HAVE_CONDOR
	return m_mw_adapter.worker_run (argc, argv, m_db_adapter);
#else

	return false;
#endif

}

bool Manager::mw_single_worker_run (Params &params_seed, Params &pseudocounts) {

#ifdef HAVE_POSTGRES
	
	assert(m_db_opts);

	m_db_opts->num_parallelized_jobs = 1;

	if ( !m_db_adapter.init_database (*m_seq_db_internal, *m_db_opts) ) 
		THROWEXPR ("ERROR: We can not initialize database server. exit");

	init_single_worker(m_db_opts->fsa, &m_db_adapter, params_seed, pseudocounts);
	if (m_db_opts->anchored)
		build_anchored_multiple_alignment ();
	else
		build_multiple_alignment ();

	m_db_adapter.look_up_data(*m_seq_db_internal, *m_db_opts);

	return true;
#else
	return false;
#endif

}

void Manager::push_edge(std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges, Edge *edge) {

#ifdef HAVE_POSTGRES
	edges.push(edge);
	return;

	if ( edge->weight >= m_min_edge_weight) 
		edges.push (edge);
	else {
	
		DB_postgres *db_connection = m_db_adapter.get_connection ();

		Seq_pos_map::const_iterator seq_pos_source= edge->source->get_seq_pos_map().begin();
		Seq_pos_map::const_iterator seq_pos_dest = edge->dest->get_seq_pos_map().begin();

		// put the edge into the buffer
		char entry[MAX_DB_ENTRY_LENGTH];
		sprintf(entry,"%d\t%d\t%d\t%d\t%lf\t%f\n", 
				seq_pos_source->first, seq_pos_source->second, seq_pos_dest->first, seq_pos_dest->second, edge->weight, edge->delta);

		m_db_merged_heap_buffer.append (entry);

		// increase the number of inserted re-weighted edges by one
		++m_num_inserted_edges;

		delete edge;
	}
#else

	edges.push (edge);
#endif

}


bool Manager::get_next_edges (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
		std::vector<Seq_pos_col_map> &seq_pos_col_maps) { 

#ifdef HAVE_POSTGRES
	// that is, this method is called only when the --db-maxram option is used
	if (m_db_opts && m_db_opts->db_max_ram == 0)
		return false;

	DB_postgres *db_connection = m_db_adapter.get_connection ();

	if (!db_connection)
		return false;

	// before getting the next edges, put re-weighted edges in the merged_heap_buffer to 
	// the merged heap table
	db_connection->flush_merged_heap_buffer ( m_db_merged_heap_buffer.c_str());
	m_db_merged_heap_buffer.clear();

	// get the next edges
	if ( (db_connection->get_merged_heap (HEAP_WINDOW_SIZE, m_merged_heap_offset, m_min_edge_weight, 
					seq_pos_col_maps, edges, m_last_edge_chunk)) == DB_BAD)
		return false;
	++m_merged_heap_offset;

	return true;
#else
	
	return false;
#endif

}


int Manager::get_edges_size(std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges) {

#ifdef HAVE_POSTGRES
	if ( is_edges_available ())
		// the merged heap table is needed to be considered
		return (edges.size() + m_orig_edges_size - (m_merged_heap_offset * HEAP_WINDOW_SIZE) + m_num_inserted_edges); 
	else 
		return edges.size();
#else

	return edges.size();
#endif

}


void Manager::get_sparse_matrices (std::vector<std::vector<SparseMatrix*> >& sparse_matrices) { 

#ifdef HAVE_POSTGRES
	if (m_db_opts && m_db_opts->db_max_ram > 0) {

		DB_postgres *db_connection = m_db_adapter.get_connection ();
		int num_jobs = db_connection->get_num_jobs ();

		// create an index on each table.
		for (int i=0;i<num_jobs;i++) 
			db_connection->create_sparse_matrix_table_index(i);

		// drop merged heap table index
		//   it is much faster to insert tuples into a table that does not have an index 
		//   than the one having an index
		db_connection->drop_merged_heap_table_index();
		db_connection->delete_merged_heap_table();

		// copy to merged heap
		for (int i=0; i<num_jobs;i++) {
			db_connection->copy_to_merged_heap_table (i);
		}

		// create merged_heap index
		db_connection->create_merged_heap_table_index();

		//m_orig_edges_size = db_connection->get_merged_heap_size();
		CTAG(1,MANAGER) << "Original edge size : " << m_orig_edges_size << endl;

		check_available_sparse_matrices();
	} else 
#endif
	{

		// in the case of database || parallelization
		get_num_cells ();
		get_all_sparse_matrices (sparse_matrices);
	}

}

void Manager::get_edges (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
		std::vector<Seq_pos_col_map> &seq_pos_col_maps) {  

#ifdef HAVE_POSTGRES
	if (m_db_opts && m_db_opts->db_max_ram > 0) {
		get_next_edges (edges, seq_pos_col_maps);
	} else 
#endif
	{
		// in the case of database || parallelization
		get_all_edges (edges, seq_pos_col_maps);
	}

}

Edge* Manager::get_next_top_edge (std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges,
		std::vector<Seq_pos_col_map> &seq_pos_col_maps) { 

	Edge *edge = edges.top();

#ifdef HAVE_POSTGRES
	if (m_db_opts && m_db_opts->db_max_ram > 0) {
		// if --db-maxram option is used
		if (m_last_edge_chunk== true )  {
			update_size(edges);
		}

		if (m_last_edge_chunk == false && edge->weight < m_min_edge_weight || edges.size() == 1) {
			// in the case of getting next edges
			get_next_edges (edges, seq_pos_col_maps);
			edge = edges.top();
		}
	}
#endif

    edges.pop();
	
	return edge;

}

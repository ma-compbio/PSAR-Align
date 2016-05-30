/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include <string.h>

#include "manager/db_postgres.h"
#include "annealing/alignment_DAG.h"


using namespace fsa;

/// connect_db
bool DB_postgres::connect_db (const char *hostname, const char *hostaddr, const char* dbname, 
		const int port, const char *user, const char *password) {

	char conninfo[DB_MAX_QUERY_LENGTH];
	memset( conninfo, 0, sizeof(conninfo) );

	// set conninfo string 
	if ( strlen(hostname) > 0 )
		sprintf( conninfo, "host = '%s'", hostname );
	if ( strlen(hostaddr) > 0 )
		sprintf( conninfo, "%s hostaddr = '%s'", conninfo, hostaddr );
	if ( strlen(dbname) > 0 )
		sprintf( conninfo, "%s dbname = '%s'", conninfo, dbname);
	if ( port > 0 )
		sprintf( conninfo, "%s port = %d", conninfo, port);
	if ( strlen(user) > 0 )
		sprintf( conninfo, "%s user = '%s'", conninfo, user);
	if ( strlen(password) > 0 )
		sprintf( conninfo, "%s password = '%s'", conninfo, password);

	// make a connection to the database	
	m_conn = PQconnectdb(conninfo);

	// check to see that the backend connection was successfully made 
	if ( PQstatus(m_conn) != CONNECTION_OK ) {
		CTAG(5,DB) << "Connection to database '" << PQdb(m_conn) << "' failed." << endl;
		CTAG(3,DB) << PQerrorMessage(m_conn) << endl;
		PQfinish(m_conn);

		return DB_BAD;
	}

	CTAG(9,DB) << "Connection to database '" << PQdb(m_conn) << "' succeeded." << endl;

	return DB_OK;
}


int DB_postgres::get_seqs_schema_id () {
	return m_seqs_schema_id;
}

int DB_postgres::get_params_table_id () {
	return m_params_table_id;
}

int DB_postgres::get_seqs_schema_id (const uint32_t hash_key, const int num_seqs, const float avg_length) {
	//TODO: what happen if nTuples>1

	int num_tuples, seqs_schema_id;

	// make a query
	sprintf(m_query, "SELECT id FROM %s WHERE hash = %u AND num_seqs = %d AND avg_length = %f",
			m_fsa_table.c_str(), hash_key, num_seqs, avg_length);

	// execute the query
	if (execute_query() == DB_BAD)
		return DB_NOT_AVAILABLE;

	// get the number of tuples
	num_tuples = PQntuples(m_res);

	// error check
	if (num_tuples == 0) {
		PQclear (m_res);
		return DB_NOT_AVAILABLE;
	}

	// get the sequence schema id
	seqs_schema_id = atoi (PQgetvalue (m_res, 0, 0));

	PQclear(m_res);

	return seqs_schema_id;
}


int DB_postgres::get_params_table_id (const DB_opts &db_opts) {

	int num_tuples, params_table_id;

	// make a query
	sprintf(m_query, "SELECT id, num_jobs FROM %s WHERE learn_gap = %s AND learn_emit_all = %s AND regularize = %s AND anchored = %s AND num_refinement_steps = %d",
			m_params_table.c_str(), (db_opts.learn_gap)? "true" : "false", (db_opts.learn_emit_all)? "true" : "false",
			(db_opts.regularize)? "true" : "false", (db_opts.anchored)? "true" : "false", db_opts.num_refinement_steps);

	// execute the query
	if (execute_query() == DB_BAD)
		return DB_NOT_AVAILABLE;

	// get the number of tuples	
	num_tuples = PQntuples(m_res);

	// error check
	if ( num_tuples == 0 ) {
		PQclear(m_res);
		return DB_NOT_AVAILABLE;
	}

	// get the parameter table id and the number of workers that have been used when generating the data
	params_table_id = atoi (PQgetvalue (m_res, 0, 0));
	m_num_jobs = atoi (PQgetvalue (m_res, 0, 1));

	PQclear(m_res);

	return params_table_id;
}

int DB_postgres::get_num_jobs () {
	return m_num_jobs;
}

bool DB_postgres::update_fsa_table() {

	bool res;

	// make a query
	sprintf(m_query, "UPDATE %s SET last_modified = CURRENT_TIMESTAMP WHERE id = %d", 
			m_fsa_table.c_str(), m_seqs_schema_id);

	// execute the query
	if (res = execute_query())
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_fsa_table () {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s (id SERIAL, hash BIGINT, num_seqs INT, avg_length INT, first_generated TIMESTAMP with time zone, last_modified TIMESTAMP with time zone)",
			m_fsa_table.c_str());

	// execute the query
	if (res = execute_query())
		PQclear(m_res);

	return res;
}

bool DB_postgres::copy_to_merged_heap_table (const int id) {

	bool res;

	// make a query
	sprintf(m_query, "INSERT INTO %s (seq1, pos1, seq2, pos2, weight, delta ) SELECT * FROM %s_%d",
			m_merged_heap_table.c_str(), m_heap_table.c_str(), id);

	// exeucte the query
	if ( res = execute_query () )
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_merged_heap_table_index () {

	bool res;

	// make a query
	sprintf(m_query, "CREATE INDEX %s_%s_idx ON %s (weight DESC, seq1, seq2, pos1, pos2)",
			m_seqs_schema.c_str(), m_merged_heap_prefix.c_str(), m_merged_heap_table.c_str());

	// execute the query
	res = execute_query();

	if (res)
		PQclear(m_res);

	return res;
}

bool DB_postgres::insert_fsa_schema (const uint32_t &hash_key, const int &num_seqs, const float &avg_length) {

	bool res;

	// make a query
	sprintf(m_query, "INSERT INTO %s (hash, num_seqs, avg_length, first_generated, last_modified) VALUES (%u, %d, %f, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)",
			m_fsa_table.c_str(), hash_key, num_seqs, avg_length);

	// execute the query
	if (res = execute_query())
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_seqs_schema() {

	bool res;

	// make a query
	sprintf(m_query, "CREATE SCHEMA %s", m_seqs_schema.c_str());	

	// execute the query
	if (res = execute_query())
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_params_table() {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s (id SERIAL, learn_gap BOOLEAN, learn_emit_all BOOLEAN, regularize BOOLEAN, anchored BOOLEAN, num_refinement_steps INT, first_generated TIMESTAMP with time zone, last_modified TIMESTAMP with time zone, num_jobs INT)", 
			m_params_table.c_str());

	// execute the query
	if (res = execute_query())
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_seqs_table() {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s (id SERIAL, sequence TEXT, length INT, hash BIGINT)",
			m_seqs_table.c_str());

	// execute the query
	if (res = execute_query())
		PQclear(m_res);

	return res;
}

bool DB_postgres::insert_seqs_table(const Sequence_database &seq_db_internal) {

	sstring seqs_buffer;
	sstring	table_name;

	seqs_buffer.clear();
	table_name.clear();

	table_name = m_seqs_table.c_str();
	table_name += " (sequence, length, hash)";

	// make a buffer to contain input sequence strings
	for (int i = 0; i < seq_db_internal.size(); i++) {
	  seqs_buffer.append(seq_db_internal.get_seq (i).seq);
	  sprintf(m_exec, "\t%d\t%lu\n", seq_db_internal.get_seq (i).length(), Hash_functions::hsieh_hash (seq_db_internal.get_seq (i).seq.c_str()));
		seqs_buffer.append(m_exec);
	}

	// put the input sequence strings into the database
	return copy_stdin (table_name.c_str(), seqs_buffer.c_str());
}

bool DB_postgres::get_merged_heap (const int size, const int offset, double &min_edge_weight, 
		std::vector<Seq_pos_col_map> &seq_pos_col_maps, std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges, bool &last) {

	// make a query
	sprintf(m_query,"SELECT seq1, pos1, seq2, pos2, weight, delta FROM %s ORDER BY WEIGHT DESC LIMIT %d OFFSET %d",
			m_merged_heap_table.c_str(), size, offset * size);

	// execute the query
	if (!execute_query())
		return DB_BAD;

	int num_tuples, seq1, pos1, seq2, pos2;
	double weight;
	float delta;
	Edge *edge = NULL;

	// get the number of tuples
	num_tuples = PQntuples(m_res);

	// if the number of result tuples is less than size, then 
	// it means that the reult is the last set of tuples 
	if (num_tuples < size)
		last = true;

	// error check
	if (num_tuples == 0) {
		PQclear (m_res);
		return DB_BAD;
	}

	// put edges into the priority queue
	for (int i = 0; i < num_tuples; i++) {
		seq1   = atoi (PQgetvalue (m_res, i, 0));
		pos1   = atoi (PQgetvalue (m_res, i, 1));
		seq2   = atoi (PQgetvalue (m_res, i, 2));
		pos2   = atoi (PQgetvalue (m_res, i, 3));
		weight = atof (PQgetvalue (m_res, i, 4));
		delta  = atof (PQgetvalue (m_res, i, 5));

		// min_edge_weight == 0.0: initial condition
		// find the smallest weight
		if (min_edge_weight == 0.0 || min_edge_weight > weight) 
			min_edge_weight = weight;

		// create a new edge and push it into the queue
		edge = new Edge(seq_pos_col_maps[seq1][pos1], seq_pos_col_maps[seq2][pos2], 
				std::pair<int, int> (seq1, pos1), std::pair<int, int> (seq2, pos2), weight, delta, 2);
		edges.push(edge);
	}

	PQclear(m_res);

	return DB_OK;
}

int DB_postgres::get_heaps (const int worker_id, double &min_edge_weight, std::vector<Seq_pos_col_map> &seq_pos_col_maps, 
		std::priority_queue<Edge*, std::vector<Edge*>, smaller_weight> &edges) {

	sstring	table_name;
	std::stringstream sstream;

	table_name.clear();
	sstream << worker_id;
	
	// set the table name
	table_name = m_heap_table.c_str();
	table_name += "_";
	table_name += sstream.str();

	// move tuples in database tables to the standard output
	if ( copy_stdout (table_name.c_str()) != DB_OK)
		return DB_BAD;

	char *buf;
	int seq1, seq2, pos1, pos2;
	double weight;
	float delta;
	Edge *edge = NULL;

	// put edges into the priority queue
	if (PQgetCopyData(m_conn, &buf, 0) != -1) {
		do {
			char *pch;
			pch = strtok(buf, "\t");  seq1   = atoi(pch);
			pch = strtok(NULL, "\t"); pos1   = atoi(pch);
			pch = strtok(NULL, "\t"); seq2   = atoi(pch);
			pch = strtok(NULL, "\t"); pos2   = atoi(pch);
			pch = strtok(NULL, "\t"); weight = atof(pch);
			pch = strtok(NULL, "\t"); delta  = atof(pch);

			// deallocated the memory
			PQfreemem(buf);

			// find the smallest weight
			if ( min_edge_weight > weight) 
				min_edge_weight = weight;

			// create a new edge and push it into the queue
			edge = new Edge(seq_pos_col_maps[seq1][pos1], seq_pos_col_maps[seq2][pos2], 
					std::pair<int, int> (seq1, pos1), std::pair<int, int> (seq2, pos2), weight, delta, 2);
			edges.push(edge);

		} while ( (PQgetCopyData(m_conn, &buf, 0)) != -1);
	} else
		return DB_NOT_AVAILABLE;

	return DB_OK;
}

int DB_postgres::get_sparse_matrices (const int worker_id, const Sequence_database &seq_db_internal, 
		const std::vector<std::vector<int> > &num_cells, std::vector<std::vector<SparseMatrix*> > &sparse_matrices) {

	sstring	table_name;
	std::stringstream sstream;

	table_name.clear();
	sstream << worker_id;

	// set the table name
	table_name = m_sparse_matrix_table.c_str();
	table_name += "_";
	table_name += sstream.str();

	// move tuples in database tables to the standard output
	if ( copy_stdout (table_name.c_str()) != DB_OK)
		return DB_BAD;

	char *buf;
	int seq1, seq2, pos1, pos2, offset = 0;
	int ex_seq1 = -1, ex_seq2 = -1;
	float prob;

	int cnt = 0;

	// put edges into the priority queue
	if ( (PQgetCopyData (m_conn, &buf, 0 )) != -1) {

		do {
			++cnt;
			char *pch;
			pch = strtok(buf, "\t");  seq1 = atoi(pch);
			pch = strtok(NULL, "\t"); pos1 = atoi(pch);
			pch = strtok(NULL, "\t"); seq2 = atoi(pch);
			pch = strtok(NULL, "\t"); pos2 = atoi(pch);
			pch = strtok(NULL, "\t"); prob = atof(pch);

			// deallocated the memory
			PQfreemem(buf);

			if ( sparse_matrices[seq1][seq2] == NULL ) {

				// once the sparse matrix of the sequence pair (ex_seq1, ex_seq2) has been constructed, 
				// create a new sparse matrix of sequence pair (ex_seq2, ex_seq1)

				if (ex_seq1 != -1 && ex_seq2 != -1)  
					sparse_matrices[ex_seq2][ex_seq1] = sparse_matrices[ex_seq1][ex_seq2]->ComputeTranspose();

				offset = 0;

				int seq1Length = seq_db_internal.get_seq (seq1).length();
				int seq2Length = seq_db_internal.get_seq (seq2).length();

				ex_seq1 = seq1;	
				ex_seq2 = seq2;

				// create the sparse matrix for the sequence pair (seq1, seq2)
				sparse_matrices[seq1][seq2] = new SparseMatrix (seq1, seq2,
										seq1Length, seq2Length,
										num_cells[seq1][seq2] - seq1Length - seq2Length);
			}

			// insert a pairwise posterior probabitlity into the sparse matrix
			sparse_matrices[seq1][seq2]->add_probability (pos1, pos2, prob, offset);

			if (pos1 != 0 && pos2 != 0 ) 
				offset++;

		} while ( (PQgetCopyData(m_conn, &buf, 0)) != -1);

	}
	else 
		return DB_NOT_AVAILABLE;

	// compute the transpose for the last sequence pair
	sparse_matrices[ex_seq2][ex_seq1] = sparse_matrices[ex_seq1][ex_seq2]->ComputeTranspose();

	return DB_OK;
}

bool DB_postgres::get_sparse_matrix (const int worker_id, const int seq1, const int seq2, 
		const int seq1Length, const int seq2Length, 
		std::vector<std::vector<SparseMatrix*> > &sparse_matrices) {

	// make a query
	sprintf(m_query,"SELECT pos1, pos2, prob FROM %s_%d where seq1 = %d and seq2 = %d order by pos1, pos2", 
			m_sparse_matrix_table.c_str(), worker_id, seq1, seq2);

	// execute the query
	if (!execute_query())
		return DB_BAD;

	int num_tuples, pos1, pos2, offset = 0;
	float prob;

	// get the number of tuples
	num_tuples = PQntuples(m_res);

	if (num_tuples == 0) {
		PQclear (m_res);
		return DB_BAD;
	}

	// create the sparse matrix for the sequence pair (seq1, seq2)
	sparse_matrices[seq1][seq2] = new SparseMatrix (seq1, seq2,
							seq1Length, seq2Length,
							num_tuples - seq1Length - seq2Length);

	for (int i = 0; i < num_tuples; i++) {
		pos1 = atoi (PQgetvalue (m_res, i, 0));
		pos2 = atoi (PQgetvalue (m_res, i, 1));
		prob = atof (PQgetvalue (m_res, i, 2));

		// insert a pairwise posterior probability to the sparse matrix
		sparse_matrices[seq1][seq2]->add_probability (pos1, pos2, prob, offset);

		if (pos1 != 0 && pos2 != 0)
			offset++;
	}	

	PQclear (m_res);

	return DB_OK;
}

int DB_postgres::get_list_of_available_pairs (const int worker_id, std::vector<std::vector<int> > &available_sparse_matrices, 
		const Sequence_database *seq_db_internal, double &avg_sparse_matrix_size, int &num_pairs, int &orig_edges_size) {

	sstring	table_name;
	std::stringstream sstream;

	table_name.clear();
	sstream << worker_id;

	// set the table name
	table_name = m_num_cells_table.c_str();
	table_name += "_";
	table_name += sstream.str();

	// move tuples in the database tables to the standard output
	if ( copy_stdout (table_name.c_str()) != DB_OK)
		return DB_BAD;

	char *buf;
	int seq1, seq2, nCells;

	double avg_size = 0.0;
	int cnt_pairs = 0;

	if (PQgetCopyData(m_conn, &buf, 0)!= -1) {
		do {	
			char *pch;
			pch = strtok(buf, "\t");  seq1 = atoi(pch);
			pch = strtok(NULL, "\t"); seq2 = atoi(pch);
			pch = strtok(NULL, "\t"); nCells = atoi(pch);

			// deallocate the memory
			PQfreemem(buf);

			int seq1Length = seq_db_internal->get_seq (seq1).length();
			int seq2Length = seq_db_internal->get_seq (seq2).length();

			if (seq1 == 0 && seq2 == 0 )
				orig_edges_size += nCells;
			else {
				cnt_pairs++;

				// flag the availability of the sparse matrix (seq1, seq2) on by setting the worker id
				// which has generated the sparse matrix
				available_sparse_matrices[seq1][seq2] = worker_id;
				avg_size += sizeof (int) * (seq1Length + 1) + sizeof (int) * (seq1Length + 1) + 
					sizeof(float) * (seq1Length + seq2Length + 2) + (sizeof (int) + sizeof(float)) * nCells;
			}
		} while ( (PQgetCopyData(m_conn, &buf, 0)) != -1);
	} else 
		return DB_NOT_AVAILABLE;

	// compute the average sparse matrix size
	avg_sparse_matrix_size = ((avg_sparse_matrix_size * (double) num_pairs ) + avg_size) / ((double) (num_pairs + cnt_pairs));

	// compute the sum of available sequence pairs	
	num_pairs += cnt_pairs;

	return DB_OK;
}

int DB_postgres::get_num_cells (const int worker_id, std::vector<std::vector<int> > &num_cells) {

	sstring	table_name;
	std::stringstream sstream;

	table_name.clear();
	sstream << worker_id;

	// set the table name
	table_name = m_num_cells_table.c_str();
	table_name += "_";
	table_name += sstream.str();

	// move the tuples in the database tables to the standard output
	if ( copy_stdout (table_name.c_str()) != DB_OK)
		return DB_BAD;

	char *buf;
	int seq1, seq2, nCells;

	if (PQgetCopyData(m_conn, &buf, 0)!= -1) {
		do {	
			char *pch;
			pch = strtok(buf, "\t");  seq1 = atoi(pch);
			pch = strtok(NULL, "\t"); seq2 = atoi(pch);
			pch = strtok(NULL, "\t"); nCells = atoi(pch);

			// deallocate the memory
			PQfreemem(buf);

			// set the number of data cells that are used for the sparse matrix (seq1, seq2)
			num_cells[seq1][seq2] += nCells;

		} while ( (PQgetCopyData(m_conn, &buf, 0)) != -1);
	} else 
		return DB_NOT_AVAILABLE;

	return DB_OK;
}

bool DB_postgres::copy_stdout (const char *table_name) {

	// make a query
	sprintf(m_query, "COPY %s TO STDOUT", table_name);

	// execute the query
	execute_query();

	PQclear(m_res);

	return DB_OK;
}

bool DB_postgres::copy_stdin (const char *table_name, const char* copy_string) {

	// make a query
	sprintf(m_query, "COPY %s FROM STDIN", table_name);

	// execute the query
	execute_query();

	PQclear(m_res);

	PQputCopyData(m_conn, copy_string, strlen(copy_string));

	PQputCopyEnd(m_conn, NULL);

	m_res = PQgetResult(m_conn);

	// error check
	if (PQresultStatus(m_res) != PGRES_COPY_IN && PQresultStatus(m_res) != PGRES_COMMAND_OK) {
		PQclear(m_res);
		return DB_BAD;
	}

	PQclear(m_res);
	return DB_OK;
}

bool DB_postgres::insert_params_table (const DB_opts &db_opts) {

	bool res;

	// make a query
	sprintf(m_query, "INSERT INTO %s (learn_gap, learn_emit_all, regularize, anchored, num_refinement_steps, first_generated, last_modified, num_jobs) VALUES (%s, %s, %s, %s, %d, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP, %d)", 
			m_params_table.c_str(), (db_opts.learn_gap)? "true" : "false", 
			(db_opts.learn_emit_all)? "true" : "false",
			(db_opts.regularize)? "true" : "false", 
			(db_opts.anchored)? "true" : "false", db_opts.num_refinement_steps, 
			db_opts.num_parallelized_jobs);

	// execute the query
	if (res = execute_query ())
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_num_cells_table (const int id) {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s_%d (seq1 INT, seq2 INT, num_cells INT)",
			m_num_cells_table.c_str(), id);

	// execute the query
	if (res = execute_query ())
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_sparse_matrix_table (const int id) {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s_%d (seq1 INT, pos1 INT, seq2 INT, pos2 INT, prob FLOAT)",
			m_sparse_matrix_table.c_str(), id);

	// execute the query
	if (res = execute_query ())
		PQclear(m_res);

	return res;
}

bool DB_postgres::delete_num_cells_table(const int id) {

	bool res;

	// make a query
	sprintf(m_query, "TRUNCATE %s_%d", m_num_cells_table.c_str(), id);

	// execute the query
	if (res = execute_query ())
		PQclear(m_res);

	return res;
}

bool DB_postgres::delete_sparse_matrix_table(const int id) {

	bool res;

	// make a query
	sprintf(m_query, "TRUNCATE %s_%d", m_sparse_matrix_table.c_str(), id);

	// execute the query
	if (res = execute_query ())
		PQclear(m_res);

	return res;
}

bool DB_postgres::update_params_table (const int num_parallelized_jobs) {

	bool res;

	// make a query
	sprintf(m_query, "UPDATE %s SET last_modified = CURRENT_TIMESTAMP, num_jobs = %d WHERE id = %d",
			m_params_table.c_str(), num_parallelized_jobs, m_params_table_id);

	// execute the query
	if (res = execute_query ())
		PQclear(m_res);

	// set the number of jobs (Note that this is equivalent to the nubmer of workers)
	m_num_jobs = num_parallelized_jobs;

	return res;
}

bool DB_postgres::drop_sparse_matrix_table_index (const int id) {

	bool res;

	// make a query
	sprintf(m_query, "DROP INDEX %s.%s_%s_%d_idx", 
			m_seqs_schema.c_str(), m_seqs_schema.c_str(), m_sparse_matrix_prefix.c_str(), id);

	// execute the query
	if (res = execute_query ())
		PQclear (m_res);

	return res;
}

bool DB_postgres::drop_merged_heap_table_index () {

	bool res;

	// make a query
	sprintf(m_query, "DROP INDEX %s.%s_%s_idx", 
			m_seqs_schema.c_str(), m_seqs_schema.c_str(), m_merged_heap_prefix.c_str());

	// execute the query
	if (res = execute_query ())
		PQclear (m_res);

	return res;
}

bool DB_postgres::create_heap_table (const int id) {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s_%d (seq1 INT, pos1 INT, seq2 INT, pos2 INT, weight DOUBLE PRECISION, delta FLOAT )",
			m_heap_table.c_str(), id);

	// execute the query
	if ( res = execute_query () )
		PQclear(m_res);

	return res;
}

bool DB_postgres::create_merged_heap_table () {

	bool res;

	// make a query
	sprintf(m_query, "CREATE TABLE %s (seq1 INT, pos1 INT, seq2 INT, pos2 INT, weight DOUBLE PRECISION, delta FLOAT ) ",
			m_merged_heap_table.c_str());

	// execute the query
	if ( res = execute_query () )
		PQclear(m_res);

	return res;
}

bool DB_postgres::delete_merged_heap_table () {

	bool res;	

	// make a query
	sprintf(m_query, "TRUNCATE %s", 
			m_merged_heap_table.c_str());

	// execute the query
	if (res = execute_query() )
		PQclear(m_res);

	return res;
}

bool DB_postgres::delete_heap_table (const int id) {

	bool res;	

	// make a query
	sprintf(m_query, "TRUNCATE %s_%d", 
			m_heap_table.c_str(),id);

	// execute the query
	if (res = execute_query() )
		PQclear(m_res);

	return res;
}

void DB_postgres::flush_merged_heap_buffer(const char* string) 
{
	// move the contents in the buffer (i.e., char* string) to the database table
	if (strlen (string) > 0 )
		copy_stdin (m_merged_heap_table.c_str(), string);
}

void DB_postgres::flush_sparse_matrix_buffer(const char* string) 
{
	// move the contents in the buffer (i.e., char* string) to the database table
	if (strlen (string) > 0) 
		copy_stdin (m_sparse_matrix_table.c_str(), string);
}

void DB_postgres::flush_num_cells_buffer(const char* string) 
{
	// move the contents in the buffer (i.e., char* string) to the database table
	if (strlen(string) > 0) 
		copy_stdin (m_num_cells_table.c_str(), string);
}

void DB_postgres::flush_heap_buffer(const char *string) 
{
	// move the contents in the buffer (i.e., char* string) to the database table
	if (strlen (string) > 0 ) {
		copy_stdin (m_heap_table.c_str(), string);
	}
}

void DB_postgres::set_worker_id (const int worker_id) {

	// append the worker id to the heap table name, the sparse matrix table name, 
	// and the num cells table name
	std::stringstream sstream;
	sstream << worker_id;

	m_heap_table += "_" + sstream.str();
	m_sparse_matrix_table += "_" + sstream.str();
	m_num_cells_table += "_" + sstream.str();
}

void DB_postgres::set_seqs_schema_id (const int seqs_schema_id) {

	// make the sequence schema name, the sequence table name, the sparse matrix info table name,
	// the heap info table name, and the parameter table name
	m_seqs_schema_id = seqs_schema_id;
	std::stringstream sstream;
	sstream << m_seqs_schema_id;

	m_seqs_schema.clear();
	m_seqs_schema = DB_FSA_SUFIX;
	m_seqs_schema += "_";
	m_seqs_schema += sstream.str();

	m_seqs_table.clear();
	m_seqs_table = m_seqs_schema;
	m_seqs_table += ".";
	m_seqs_table += DB_SEQUENCE_INFO;

	m_sparse_matrix_info.clear();
	m_sparse_matrix_info = m_seqs_schema;
	m_sparse_matrix_info += ".";
	m_sparse_matrix_info += DB_SPARSE_MATRIX_INFO;

	m_heap_info.clear();
	m_heap_info = m_seqs_schema;
	m_heap_info += ".";
	m_heap_info += DB_HEAP_INFO;

	m_params_table.clear();
	m_params_table = m_seqs_schema;
	m_params_table += ".";
	m_params_table += DB_PARAMETER_INFO;
}

void DB_postgres::set_params_table_id (const int params_table_id) {

	// make the heap_table name, the merged_heap_table name, the sparse_matrix_table name, 
	// and the num cells table name
	m_params_table_id = params_table_id;
	std::stringstream sstream;
	sstream << m_params_table_id;

	m_heap_prefix.clear();
	m_heap_prefix = DB_HEAP_SUFIX;
	m_heap_prefix += sstream.str();

	m_heap_table.clear();
	m_heap_table = m_seqs_schema;
	m_heap_table += ".";
	m_heap_table += m_heap_prefix;

	m_merged_heap_prefix.clear();
	m_merged_heap_prefix = DB_MERGED_HEAP_SUFIX;
	m_merged_heap_prefix += sstream.str();

	m_merged_heap_table.clear();
	m_merged_heap_table = m_seqs_schema;
	m_merged_heap_table += ".";
	m_merged_heap_table += m_merged_heap_prefix;

	m_sparse_matrix_prefix.clear();
	m_sparse_matrix_prefix = DB_SPARSE_MATRIX_SUFIX;
	m_sparse_matrix_prefix += sstream.str();

	m_sparse_matrix_table.clear();
	m_sparse_matrix_table = m_seqs_schema;
	m_sparse_matrix_table += ".";
	m_sparse_matrix_table += m_sparse_matrix_prefix;

	m_num_cells_table.clear();
	m_num_cells_table = m_seqs_schema;
	m_num_cells_table += ".";
	m_num_cells_table += DB_NUM_CELLS_SUFIX;
	m_num_cells_table += sstream.str();
}

bool DB_postgres::execute_query() {

	// execute the query
	m_res = PQexec(m_conn, m_query);

	// error check
	if (PQresultStatus(m_res) != PGRES_COMMAND_OK && PQresultStatus(m_res) != PGRES_TUPLES_OK
			&& PQresultStatus(m_res) != PGRES_COPY_IN && PQresultStatus(m_res) != PGRES_COPY_OUT) {

		PQclear(m_res);
		m_res = NULL;
		return DB_BAD;
	}

	return DB_OK;
}

DB_postgres::DB_postgres() {
	m_conn = NULL;

	// make the fsa_table name, which is the main table name
	m_fsa_table.clear();
	m_fsa_table = DB_FSA_SUFIX;
	m_fsa_table += "_";
	m_fsa_table += DB_FSA_TABLE_PREFIX;
}

DB_postgres::~DB_postgres() {}

void DB_postgres::disconnect_db() {
	PQfinish(m_conn);
}

bool DB_postgres::create_sparse_matrix_table_index (const int id) {

	bool res;

	// make a query
	sprintf(m_query, "CREATE INDEX %s_%s_%d_idx ON %s_%d (seq1, seq2, pos1, pos2)",
			m_seqs_schema.c_str(), m_sparse_matrix_prefix.c_str(), id, m_sparse_matrix_table.c_str(), id);

	// execute the query
	res = execute_query();

	if (res)
		PQclear(m_res);

	return res;
}

int DB_postgres::get_merged_heap_size() {
	
	bool res;

	// make a query
	sprintf(m_query, "SELECT COUNT(weight) FROM %s",
			m_merged_heap_table.c_str());

	// execute the query
	res = execute_query();

	if (res)
		PQclear(m_res);

	return res;
}

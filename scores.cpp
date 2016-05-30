#include <iostream>
#include <fstream>
#include <cassert>
#include <sys/types.h>
#include <dirent.h>
#include <cstdlib>
#include <iomanip>
#include "scores.h"

Scores::Scores(vector<Sequence>& sequences)
{
	for (int i = 0; i < sequences.size()-1; i++) {
		Sequence seq1 = sequences.at(i);
		int seq1len = seq1.length();
		for (int j = i+1; j < sequences.size(); j++) {
			Sequence seq2 = sequences.at(j);
			int seq2len = seq2.length();
			
			assert(seq1len == seq2len);

			int pos1 = 1;
			int pos2 = 1; 
			for (int col = 0; col < seq1len; col++) {
				int ch1 = seq1.at(col);
				int ch2 = seq2.at(col);
				if (ch1 != GAP && ch2 != GAP) {
					ResiduePair rp(i+1, j+1, pos1, pos2);
					_pairscores[rp] = 0.0;

					_coliter = _columns.find(col+1);
					if (_coliter == _columns.end()) {
						vector<ResiduePair> vec;	
						vec.push_back(rp);
						_columns[col+1] = vec;
					} else {
						_coliter->second.push_back(rp);
					} // end of if else
				} // end of if

				if (ch1 != GAP) pos1++;	
				if (ch2 != GAP) pos2++;	
			} // end of for s		
		} // end of for j
	} // end of for i
}

void Scores::compute_scores (string dir)
{
	// read sample filenames
	vector<string> filenames;
	DIR *dp;
	struct dirent *dirp;
	if ((dp = opendir(dir.c_str())) == NULL) {
		cerr << "\n[ERROR] Unable to read directory " << dir << endl;
		exit(1); 
	} // end of if

	while ((dirp = readdir(dp)) != NULL) {
		string fname = string(dirp->d_name);
		if (fname.at(0) != '.') filenames.push_back(dir + "/" + fname);
	}

	closedir(dp);	

	// for each alignment samples, check consistency with the input
	for (int f = 0; f < filenames.size(); f++) {
		string fname = filenames.at(f);
		vector<Sequence> sequences;
		Sequence::read_alignments(fname.c_str(), sequences);
	
		for (int i = 0; i < sequences.size()-1; i++) {
			Sequence seq1 = sequences.at(i);
			int seq1len = seq1.length();
			for (int j = i+1; j < sequences.size(); j++) {
				Sequence seq2 = sequences.at(j);
				int seq2len = seq2.length();
			
				assert(seq1len == seq2len);

				int pos1 = 1;
				int pos2 = 1; 
				for (int col = 0; col < seq1len; col++) {
					int ch1 = seq1.at(col);
					int ch2 = seq2.at(col);
					if (ch1 != GAP && ch2 != GAP) {
						ResiduePair rp(i+1, j+1, pos1, pos2);
						_pairiter = _pairscores.find(rp);
						if (_pairiter != _pairscores.end()) _pairiter->second++;
					} // end of if

					if (ch1 != GAP) pos1++;	
					if (ch2 != GAP) pos2++;	
				} // end of for s		
			} // end of for j
		} // end of for i
	} // end of for f 
	
	// compute average pair scores
	int total = filenames.size();
	for (_pairiter = _pairscores.begin(); _pairiter != _pairscores.end(); _pairiter++) 
		_pairiter->second /= (double)total;
}

void Scores::writetofile (string fnamepair, string fnamecol) 
{
	// pair scores
	fstream pairfile(fnamepair.c_str(), ios::out);
	pairfile << "#seq1\tseq2\tpos1\tpos2\tscore" << endl;
	for (_pairiter = _pairscores.begin(); _pairiter != _pairscores.end(); _pairiter++) {
		ResiduePair rp = _pairiter->first;
		double score = _pairiter->second;
		pairfile << rp.seq1 << "\t" << rp.seq2 << "\t" << rp.pos1 << "\t" << rp.pos2 << "\t" << score << endl;;
	}
	pairfile.close();
	
	// compute average column scores
	fstream colfile(fnamecol.c_str(), ios::out);
	colfile << "#col\tscore" << endl;
	for (_coliter = _columns.begin(); _coliter != _columns.end(); _coliter++) {	
		int col = _coliter->first;
		vector<ResiduePair> vec = _coliter->second;

		double sum = 0.0;
		for (int i = 0; i < vec.size(); i++) {
			ResiduePair rp = vec.at(i);
			_pairiter = _pairscores.find(rp);
			assert(_pairiter != _pairscores.end());
			sum += _pairiter->second;
		} // end of for i

		int count = vec.size();
		colfile << col << "\t" << setprecision(6) << sum/(double)count << endl; 
	} // end of for	
	colfile.close();	
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include "sequence.h"
#include "model.h"
#include "scores.h"
#include "params.h"

using namespace std;

const char *const program_name = "PSAR";
const char *const program_version = "Version 1.0";

void usage(int exit_value = 0)
{
	cerr << "[USAGE] " << program_name << " <alignment fasta file> <parameter file> <output directory> [-s <sampling threshold in (0,1)>] [-p]" << endl;
	cerr << "\t-s: Sampling threshold in (0,1), default=0.5" << endl;
	cerr << "\t-p: Print estimated parameters" << endl << endl; 
	exit(exit_value);
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	
	// check parameters
	if (argc < 4 || argc > 7) {
		cerr << "\n[ERROR] Parameter errors." << endl;
		usage(1);
	}

	char* filename = argv[1];
	char* parname = argv[2];
	char* out_dir = argv[3];	

	bool bPrintParams = false;	
	double threshold = 0.5;
	if (argc == 4) {
		;	
	} else if (argc == 5 && strcmp(argv[4], "-p") == 0 ) { 
		bPrintParams = true;
	} else if (argc == 6 && strcmp(argv[4], "-s") == 0) {
		threshold = atof(argv[5]);
		if (threshold < 0.0 || threshold > 1.0) {
			cerr << "\n[ERROR] The sample threshold should be in range [0,1]." << endl;
			usage(1);
		}
		cerr << "Sampling threshold : " << threshold << endl;
	} else if (argc == 7) {
		if (strcmp(argv[4], "-p") == 0 && strcmp(argv[5], "-s") == 0) {
			threshold = atof(argv[6]);	
		} else if (strcmp(argv[4], "-s") == 0 && strcmp(argv[6], "-p") == 0) {
			threshold = atof(argv[5]);	
		} else {
			cerr << "\n[ERROR] Parameter errors." << endl;
			usage(1);
		}
		if (threshold < 0.0 || threshold > 1.0) {
			cerr << "\n[ERROR] The sample threshold should be in range [0,1]." << endl;
			usage(1);
		}
		cerr << "Sampling threshold : " << threshold << endl;
		bPrintParams = true; 
	} else {
		cerr << "\n[ERROR] Parameter errors." << endl;
		usage(1);
	} 

	cerr << endl;
	cerr << "::::: " << program_name << " " << program_version << " :::::" << endl << endl;	
	cerr << "Alignment file : " << filename << endl;
	cerr << "Parameter file : " << parname << endl;
	cerr << "Output directory : " << out_dir << endl;
	if (threshold < 2.0) cerr << "Sampling threshold : " << threshold << endl;
	cerr << endl;
	
	// create output directory
	mkdir (out_dir, 0755);

	// read alignment
	vector<Sequence> sequences;
	Sequence::read_alignments(filename, sequences);

	// read parameters
	Params params(parname);
		
	int gapcolumncnt = Sequence::getNumofGappedColumns(sequences);
	int numSampling = 0;
	int totalgapsize = 0;
	int totalgapcnt = 0;
	map<string,string> samples;
	vector<double> vecParams(6,0.0);

	// iterate for each sequence
	for (int i = 0; i < sequences.size(); i++) {
		Sequence seqTarget = sequences.at(i);
		vector<int> gapstats = seqTarget.removeGaps();
	
		// skip if there is no gap
		if (gapstats.at(1) == 0) {
			cerr << "Skip sequence" << i+1 << " : no gaps" << endl;
			continue;
		}
		cerr << "Sampling : sequence " << i+1 << " vs others ..." << endl;
	
		totalgapsize += gapstats.at(0);	// total gap size
		totalgapcnt += gapstats.at(1);	// the number of contiguous gaps
		numSampling += gapstats.at(1);
		
		vector<Sequence> seqFamily;
		for (int j = 0; j < sequences.size(); j++) {
			if (i == j) continue;
			seqFamily.push_back(sequences.at(j));
		} // end of for

		seqFamily = Sequence::removeNullColumns(seqFamily);

		Model m(&seqTarget, &seqFamily, i, &params);
		if (params.hasTransProb()) m.estimateParameters(false, vecParams);
		else m.estimateParameters(true, vecParams);
		
		for (int sp = 0; sp < gapstats.at(0); sp++) {
			string straln;
			double lls = m.sample_alignments(straln, sequences, threshold);
			ostringstream os;
			os << exp(lls);	
			if (samples.find(straln) == samples.end()) {
				samples.insert(make_pair(straln, os.str()));
			}
		} // end of for sp
	} // end of for i

	int mapsize = samples.size();
	double avgnumsamples = (double)mapsize;

	cerr << endl << "Number of distinct samples :" << avgnumsamples << endl;

	if (bPrintParams) {
		cerr << endl << "Estimated parameters: average transition probabilities among states M, IS, and IA" << endl;
		for (int vi = 0; vi < 6; vi++) {
			vecParams[vi] /= sequences.size();
		}	
		cerr << "MtoIS = " << vecParams[0] << endl; 	
		cerr << "MtoIA = " << vecParams[1] << endl; 	
		cerr << "IStoIS = " << vecParams[2] << endl; 	
		cerr << "IAtoIA = " << vecParams[3] << endl; 	
		cerr << "IStoIA = " << vecParams[4] << endl; 	
		cerr << "IAtoIS = " << vecParams[5] << endl; 	
	}


	if (mapsize > 0) {
		int filecnt = 1;
		map<string, string>::iterator iter;
		for (iter = samples.begin(); iter != samples.end(); iter++) {
			string aln = iter->first;
			string strls = iter->second;

			ostringstream os;
			os << out_dir << "/sample" << filecnt << ".fa";
			string out_f = os.str();

			fstream outfile(out_f.c_str(), ios::out);
			outfile << aln;
			outfile.close();	

			filecnt++;	
		} // end of for iter	
	}
	if (mapsize == 0) {
		cerr << endl << "No scoring because there is no sampled alignment..." << endl;
		return 0;
	} 

	cerr << endl << "Scoring..." << endl;
	// scoring
	Scores sc(sequences);	
	sc.compute_scores(out_dir);
	sc.writetofile("PSAR_pair.txt", "PSAR_column.txt");

	cerr << endl << "Done." << endl;

	return 0;
}


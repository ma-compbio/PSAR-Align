#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <map>
#include <cmath>
#include "sequence.h"
#include "params.h"
  
using namespace std;

const double DEFAULT_DELTA = 0.02;
const double DEFAULT_EPSILON = 0.8;
const double DEFAULT_TAU = 0.001;
const double DEFAULT_RHO = 0.00001;

// HMM states
enum {M,I,D,S,E};
const int NUM_STATES = 3; // only for M,I,D
const double NEGINF = log(0);
const int NUM_ITER = 2;
const double DOUBLE_MIN = 0.000001;

class Model
{
private:
	// Sequences
	Sequence* _seqTarget;
	vector<Sequence>* _seqFamily;
	Params* _params;
	int _numrows;
	int _numcols;
	int _familysize;
	int _targetIndex;

	// For emission probabilities (F81 model)
	double _time;
	double _bkgProbs[4];
	double _logBkgProbs[4];
	double _substMatrix[4][4];	
	double _logSubstMatrix[4][4];	

	map<string, double> _emissionProbs;

	// For transition probabilities (3 states pair-HMM)
	double _tpMI, _tpMD;
	double _tpII, _tpDD;
	double _tpID, _tpDI;
	double _tpSD, _tpSI;
	double _tp2E;
	double _logTransMatrix[NUM_STATES+2][NUM_STATES+2];

	// Dynamic programming tables
	double** _forwardTable[NUM_STATES]; // only for M,I,D
	double** _backwardTable[NUM_STATES]; // only for M,I,D
	double _logfullprob;

	// Private functions
	void updateTransMatrix();
	double logsum(double logx, double logy);
	double getLogJointProb(int r, int c);
	double getLogJointProbSingle(int c);
	double getNextTableIndex(double mprob, double iprob, double dprob, double threshold);
	
public:
	//Model() {}
	Model(Sequence* seqTarget, vector<Sequence>* seqFamily, int targetIndex, Params* params);
	~Model();

	void estimateParameters(bool est_param, vector<double>& vecParams);
	void forward();
	void backward();
	void dump();
	double sample_alignments(string& straln, vector<Sequence>& sequences, double threshold);
};

#endif

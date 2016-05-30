#include "model.h"
#include <iostream>
#include <sstream>
#include <cstdlib>

Model::Model(Sequence* seqTarget, vector<Sequence>* seqFamily, int targetIndex, Params* params)
{
	_seqTarget = seqTarget;
	_seqFamily = seqFamily;
	_params = params;
	_numrows = _seqTarget->length();
	_numcols = _seqFamily->at(0).length();
	_familysize = _seqFamily->size();
	_targetIndex = targetIndex;

	_time = _params->getNumSubst();
	_bkgProbs[0] = _params->getProbA();	// A
	_bkgProbs[1] = _params->getProbC(); // C
	_bkgProbs[2] = _params->getProbG(); // G
	_bkgProbs[3] = _params->getProbT(); // T
	for (int i = 0; i < 4; i++) _logBkgProbs[i] = log(_bkgProbs[i]);

	// F81 substituion model
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) _substMatrix[i][j] = exp(-_time) + (1.0 - exp(-_time)) * _bkgProbs[j];
			else _substMatrix[i][j] = (1.0 - exp(-_time)) * _bkgProbs[j];
		} // end of for j
	} // end of for i
	for (int i = 0; i < 4; i++) 
		for (int j = 0; j < 4; j++) 
			_logSubstMatrix[i][j] = log(_substMatrix[i][j]);

	double val = _params->getMtoI();
	_tpMI = (val != -1) ? val : DEFAULT_DELTA; 
	val = _params->getMtoD();
	_tpMD = (val != -1) ? val : DEFAULT_DELTA;
	val = _params->getItoI();
	_tpII = (val != -1) ? val : DEFAULT_EPSILON;
	val = _params->getDtoD();
	_tpDD = (val != -1) ? val : DEFAULT_EPSILON;
	val = _params->getItoD();
	_tpID = (val != -1) ? val : DEFAULT_RHO;
	val = _params->getDtoI();
	_tpDI = (val != -1) ? val : DEFAULT_RHO;
	_tpSI = _tpSD = DEFAULT_DELTA;
	_tp2E = DEFAULT_TAU;

	updateTransMatrix();

	// Allocate memory for DP tables
	for (int i = 0; i < NUM_STATES; i++) {
		_forwardTable[i] = new double* [_numrows+1];
		_backwardTable[i] = new double* [_numrows+1];
		for (int j = 0; j < _numrows+1; j++) {
			_forwardTable[i][j] = new double [_numcols+1];
			_backwardTable[i][j] = new double [_numcols+1];
		} // end of for j
	} // end of for i
}

Model::~Model()
{
	for (int i = 0; i < NUM_STATES; i++) {
		for (int j = 0; j < _numrows+1; j++) {
			delete [] _forwardTable[i][j];
			delete [] _backwardTable[i][j];
		} // end of for j
		delete [] _forwardTable[i];
		delete [] _backwardTable[i];
	} // end of for i
}

void Model::dump()
{
	cerr << "t=" << _time << endl;
	cerr << "bkg prob=(" << _bkgProbs[0] << "," << _bkgProbs[1] << "," << _bkgProbs[2] << "," << _bkgProbs[3] << ")" << endl;
	cerr << "subst prob:" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) cerr << _substMatrix[i][j] << "\t";
		cerr << endl;
	}
	cerr << "trans prob:" << endl;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) cerr << exp(_logTransMatrix[i][j]) << "\t";
		cerr << endl;
	}
	cerr << "forward table:M" << endl;
	for (int i = 0; i <= _numrows; i++) {
		for (int j = 0; j <= _numcols; j++) {
			cerr << exp(_forwardTable[M][i][j]) << "\t";
		}
		cerr << endl;	
	}	
	cerr << "forward table:I" << endl;
	for (int i = 0; i <= _numrows; i++) {
		for (int j = 0; j <= _numcols; j++) {
			cerr << exp(_forwardTable[I][i][j]) << "\t";
		}
		cerr << endl;	
	}	
	cerr << "forward table:D" << endl;
	for (int i = 0; i <= _numrows; i++) {
		for (int j = 0; j <= _numcols; j++) {
			cerr << exp(_forwardTable[D][i][j]) << "\t";
		}
		cerr << endl;	
	}	
}
	

void Model::updateTransMatrix()
{
	double logtpMI = log(_tpMI);
	double logtpMD = log(_tpMD);
	double logtpII = log(_tpII);
	double logtpDD = log(_tpDD);
	double logtpID = log(_tpID);
	double logtpDI = log(_tpDI);
	double logtpSI = log(_tpSI);
	double logtpSD = log(_tpSD);
	double logtp2E = log(_tp2E);

	for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            _logTransMatrix[i][j] = NEGINF;

	_logTransMatrix[M][M] = log(1.0-_tpMI-_tpMD-_tp2E);
	_logTransMatrix[M][I] = logtpMI; 
	_logTransMatrix[M][D] = logtpMD;	
	_logTransMatrix[M][E] = logtp2E;

	_logTransMatrix[I][M] = log(1.0-_tpII-_tpID-_tp2E);	
	_logTransMatrix[I][I] = logtpII;	
	_logTransMatrix[I][D] = logtpID;	
	_logTransMatrix[I][E] = logtp2E;
	
	_logTransMatrix[D][M] = log(1.0-_tpDD-_tpDI-_tp2E);	
	_logTransMatrix[D][D] = logtpDD;	
	_logTransMatrix[D][I] = logtpDI;	
	_logTransMatrix[D][E] = logtp2E;

	_logTransMatrix[S][M] = _logTransMatrix[M][M];
    _logTransMatrix[S][I] = _logTransMatrix[M][I];
    _logTransMatrix[S][D] = _logTransMatrix[M][D];
    _logTransMatrix[S][E] = logtp2E;
}

void Model::estimateParameters(bool est_param, vector<double>& vecParams)
{
	if (est_param == false) {
		forward();
		backward();	
		return;
	}

	double prevscore = NEGINF;
	for (int it = 0; it < NUM_ITER; it++) {
		forward();
		backward();	

		// expected transition counts
		double logTrCnts[NUM_STATES][NUM_STATES];
		for (int i = 0; i < NUM_STATES; i++) 
			for (int j = 0; j < NUM_STATES; j++) logTrCnts[i][j] = NEGINF;
 
		for (int i = 0; i < _numrows; i++) {
			for (int j = 0; j < _numcols; j++) {
				//MtoM
				double logval = _forwardTable[M][i][j] + _logTransMatrix[M][M] + getLogJointProb(i+1, j+1) + _backwardTable[M][i+1][j+1]; 
				logTrCnts[M][M] = logsum(logTrCnts[M][M], logval);
			
				//MtoI	
				logval = _forwardTable[M][i][j] + _logTransMatrix[M][I] + _logBkgProbs[_seqTarget->at(i)] + _backwardTable[I][i+1][j]; 
				logTrCnts[M][I] = logsum(logTrCnts[M][I], logval);	 
			
				//MtoD	
				logval = _forwardTable[M][i][j] + _logTransMatrix[M][D] + getLogJointProbSingle(j+1) + _backwardTable[D][i][j+1]; 
				logTrCnts[M][D] = logsum(logTrCnts[M][D], logval);	 

				//ItoM
				logval = _forwardTable[I][i][j] + _logTransMatrix[I][M] + getLogJointProb(i+1, j+1) + _backwardTable[M][i+1][j+1];
				logTrCnts[I][M] = logsum(logTrCnts[I][M], logval);	
				
				//ItoI
				logval = _forwardTable[I][i][j] + _logTransMatrix[I][I] + _logBkgProbs[_seqTarget->at(i)] + _backwardTable[I][i+1][j];
				logTrCnts[I][I] = logsum(logTrCnts[I][I], logval);	
				
				//ItoD
				logval = _forwardTable[I][i][j] + _logTransMatrix[I][D] + getLogJointProbSingle(j+1) + _backwardTable[D][i][j+1];
				logTrCnts[I][D] = logsum(logTrCnts[I][D], logval);	
				
				//DtoM
				logval = _forwardTable[D][i][j] + _logTransMatrix[D][M] + getLogJointProb(i+1, j+1) + _backwardTable[M][i+1][j+1];
				logTrCnts[D][M] = logsum(logTrCnts[D][M], logval);	
				
				//DtoD	
				logval = _forwardTable[D][i][j] + _logTransMatrix[D][D] + getLogJointProbSingle(j+1) + _backwardTable[D][i][j+1]; 
				logTrCnts[D][D] = logsum(logTrCnts[D][D], logval);	 
				
				//DtoI	
				logval = _forwardTable[D][i][j] + _logTransMatrix[D][I] + _logBkgProbs[_seqTarget->at(i)] + _backwardTable[I][i+1][j]; 
				logTrCnts[D][I] = logsum(logTrCnts[D][I], logval);	 
			} // end of for j
		} // end of for i

		for (int i = 0; i < NUM_STATES; i++)  
			for (int j = 0; j < NUM_STATES; j++) { 
				logTrCnts[i][j] = exp(logTrCnts[i][j] - _logfullprob);
			}

			
		_tpMI = logTrCnts[M][I] / (logTrCnts[M][M] + logTrCnts[M][I] + logTrCnts[M][D]);
		if (_tpMI < DOUBLE_MIN) _tpMI = DOUBLE_MIN;
		_tpMD = logTrCnts[M][D] / (logTrCnts[M][M] + logTrCnts[M][I] + logTrCnts[M][D]);
		if (_tpMD < DOUBLE_MIN) _tpMD = DOUBLE_MIN;
		
		_tpII = logTrCnts[I][I] / (logTrCnts[I][I] + logTrCnts[I][M] + logTrCnts[I][D]);
		if (_tpII < DOUBLE_MIN) _tpII = DOUBLE_MIN;
		_tpDD = logTrCnts[D][D] / (logTrCnts[D][D] + logTrCnts[D][M] + logTrCnts[D][I]);
		if (_tpDD < DOUBLE_MIN) _tpDD = DOUBLE_MIN;
	
		_tpID = logTrCnts[I][D] / (logTrCnts[I][I] + logTrCnts[I][M] + logTrCnts[I][D]);
		if (_tpID < DOUBLE_MIN) _tpID = DOUBLE_MIN;
		_tpDI = logTrCnts[D][I] / (logTrCnts[D][D] + logTrCnts[D][M] + logTrCnts[D][I]);
		if (_tpDI < DOUBLE_MIN) _tpDI = DOUBLE_MIN;

		updateTransMatrix();
	} // end of for
	
	vecParams.at(0) += _tpMI;
	vecParams.at(1) += _tpMD;
	vecParams.at(2) += _tpII;
	vecParams.at(3) += _tpDD;
	vecParams.at(4) += _tpID;
	vecParams.at(5) += _tpDI;
}

void Model::forward()
{
	// initialize a table
	_forwardTable[M][0][0] = log(1);
	_forwardTable[I][0][0] = _forwardTable[D][0][0] = NEGINF; 	

	for (int r = 0; r <= _numrows; r++) {
		for (int c = 0; c <= _numcols; c++) {
			if (r == 0 && c == 0) continue;
			double Mr_1c_1 = (r-1 > -1 && c-1 > -1) ? _forwardTable[M][r-1][c-1] : NEGINF;	
			double Mr_1c = (r-1 > -1) ? _forwardTable[M][r-1][c] : NEGINF; 			
			double Mrc_1 = (c-1 > -1) ? _forwardTable[M][r][c-1] : NEGINF;	
			
			double Ir_1c_1 = (r-1 > -1 && c-1 > -1) ? _forwardTable[I][r-1][c-1] : NEGINF;	
			double Ir_1c = (r-1 > -1) ? _forwardTable[I][r-1][c] : NEGINF; 			
			double Irc_1 = (c-1 > -1) ? _forwardTable[I][r][c-1] : NEGINF; 			
			
			double Dr_1c_1 = (r-1 > -1 && c-1 > -1) ? _forwardTable[D][r-1][c-1] : NEGINF;	
			double Dr_1c = (r-1 > -1) ? _forwardTable[D][r-1][c] : NEGINF;
			double Drc_1 = (c-1 > -1) ? _forwardTable[D][r][c-1] : NEGINF;
	
			double logdouble = getLogJointProb(r, c);
			
			double logval = logsum(_logTransMatrix[M][M]+Mr_1c_1, _logTransMatrix[I][M]+Ir_1c_1);
			logval = logsum(logval, _logTransMatrix[D][M]+Dr_1c_1);
			_forwardTable[M][r][c] = logdouble + logval;

			double logsingle = (r-1 > -1) ? _logBkgProbs[_seqTarget->at(r-1)] : NEGINF;
			logval = logsum(_logTransMatrix[M][I]+Mr_1c, _logTransMatrix[I][I]+Ir_1c);
			logval = logsum(logval, _logTransMatrix[D][I]+Dr_1c);
			_forwardTable[I][r][c] = logsingle + logval;

			logsingle = getLogJointProbSingle(c);	
			logval = logsum(_logTransMatrix[M][D]+Mrc_1, _logTransMatrix[D][D]+Drc_1);
			logval = logsum(logval, _logTransMatrix[I][D]+Irc_1);
			_forwardTable[D][r][c] = logsingle + logval;
		} // end of for c
	} // end of for r

	// compute full probability
	double logfullprob = _logTransMatrix[M][E]+_forwardTable[M][_numrows][_numcols];
	logfullprob = logsum(logfullprob, _logTransMatrix[M][E]+_forwardTable[I][_numrows][_numcols]);
	logfullprob = logsum(logfullprob, _logTransMatrix[M][E]+_forwardTable[D][_numrows][_numcols]);
	
	_logfullprob = logfullprob;	
}

void Model::backward()
{
	// initialize a table
	_backwardTable[M][_numrows][_numcols] = log(_tp2E);
	_backwardTable[I][_numrows][_numcols] = log(_tp2E);
	_backwardTable[D][_numrows][_numcols] = log(_tp2E);

	for (int r = _numrows; r >= 0; r--) {
		for (int c = _numcols; c >= 0; c--) {
			if (r == _numrows && c == _numcols) continue;
			double Mr__1c__1 = (r+1 <= _numrows && c+1 <= _numcols) ? _backwardTable[M][r+1][c+1] : NEGINF;
            double Mr__1c = (r+1 <= _numrows) ? _backwardTable[M][r+1][c] : NEGINF;
            double Mrc__1 = (c+1 <= _numcols) ? _backwardTable[M][r][c+1] : NEGINF;
            
            double Ir__1c__1 = (r+1 <= _numrows && c+1 <= _numcols) ? _backwardTable[I][r+1][c+1] : NEGINF;
            double Ir__1c = (r+1 <= _numrows) ? _backwardTable[I][r+1][c] : NEGINF;
            double Irc__1 = (c+1 <= _numcols) ? _backwardTable[I][r][c+1] : NEGINF;
            
            double Dr__1c__1 = (r+1 <= _numrows && c+1 <= _numcols) ? _backwardTable[D][r+1][c+1] : NEGINF;
            double Dr__1c = (r+1 <= _numrows) ? _backwardTable[D][r+1][c] : NEGINF;
            double Drc__1 = (c+1 <= _numcols) ? _backwardTable[D][r][c+1] : NEGINF;

			double logdouble = getLogJointProb(r+1, c+1);
			double logsingle_r = (r+1 <= _numrows) ? _logBkgProbs[_seqTarget->at(r)] : NEGINF;
			double logsingle_c = getLogJointProbSingle(c+1);	

			double logval = logsum(logdouble+_logTransMatrix[M][M]+Mr__1c__1, logsingle_r+_logTransMatrix[M][I]+Ir__1c);
			logval = logsum(logval, logsingle_c+_logTransMatrix[M][D]+Drc__1);
			_backwardTable[M][r][c] = logval;

			logval = logsum(logdouble+_logTransMatrix[I][M]+Mr__1c__1, logsingle_r+_logTransMatrix[I][I]+Ir__1c);
			logval = logsum(logval, logsingle_c+_logTransMatrix[I][D]+Irc__1);
			_backwardTable[I][r][c] = logval;		
			
			logval = logsum(logdouble+_logTransMatrix[D][M]+Mr__1c__1, logsingle_c+_logTransMatrix[D][D]+Drc__1);
			logval = logsum(logval, logsingle_r+_logTransMatrix[D][I]+Ir__1c);
			_backwardTable[D][r][c] = logval;		
		} // end of for c
	} // end of for r 
}

double Model::logsum(double logx, double logy)
{
	if (logx <= NEGINF && logy <= NEGINF) return NEGINF;
	if (logx <= NEGINF) return logy;
	if (logy <= NEGINF) return logx;
	if (logy >= logx) return (logy + log(1.0 + exp(logx - logy)));
	return (logx + log(1.0 + exp(logy - logx))); 
}

double Model::getLogJointProb(int r, int c)
{
	// Joint probability of an alignment column based on a star topology
	r--;
	c--;
	if (r < 0 || c < 0 || r >= _numrows || c >= _numcols) return NEGINF;

	// collect non-gap characters
	vector<int> bases;
	bases.push_back(_seqTarget->at(r));
	for (int i = 0; i < _familysize; i++) {
		int ch = _seqFamily->at(i).at(c);
		if (ch == GAP) continue;
		bases.push_back(ch);
	} // end of for i 

	// search cache
	stringstream ss;
	for (int i = 0; i < bases.size(); i++) ss << bases.at(i);
	string key = ss.str();
	map<string, double>::iterator iter = _emissionProbs.find(key);
	if (iter != _emissionProbs.end()) {
		return iter->second;
	} // end of if

	double logprob = NEGINF;
	for (int ai = 0; ai < 4; ai++) { // for all possible ancestral characters
		double val = _logBkgProbs[ai];
		for (int di = 0; di < bases.size(); di++) {
			int descChar = bases.at(di);
			if (descChar == N) continue;
			double descBkgProb = _logBkgProbs[descChar];
			val += _logSubstMatrix[ai][descChar];
		} // end of for di
		logprob = logsum(logprob, val);
	} // end of for ai

	_emissionProbs[key] = logprob;
	return logprob;
}

double Model::getLogJointProbSingle(int c)
{
	// Joint probability of an alignment column based on a star topology
	c--;
	if (c < 0 || c >= _numcols) return NEGINF;

	// collect non-gap characters
	vector<int> bases;
	for (int i = 0; i < _familysize; i++) {
		int ch = _seqFamily->at(i).at(c);
		if (ch == GAP) continue;
		bases.push_back(ch);
	} // end of for i 

	// search cache
	stringstream ss;
	for (int i = 0; i < bases.size(); i++) ss << bases.at(i);
	string key = ss.str();
	map<string, double>::iterator iter = _emissionProbs.find(key);
	if (iter != _emissionProbs.end()) {
		return iter->second;
	} // end of if

	double logprob = NEGINF;
	for (int ai = 0; ai < 4; ai++) { // for all possible ancestral characters
		double val = _logBkgProbs[ai];
		for (int di = 0; di < bases.size(); di++) {
			int descChar = bases.at(di);
			if (descChar == N) continue;
			double descBkgProb = _logBkgProbs[descChar];
			val += _logSubstMatrix[ai][descChar];
		} // end of for di
		logprob = logsum(logprob, val);
	} // end of for ai

	_emissionProbs[key] = logprob;
	return logprob;
}

double Model::sample_alignments(string& straln, vector<Sequence>& sequences, double threshold)
{
	string aln_target;
	vector<string> aln_family(_familysize);

	double lls = 0;	// log likelihood score

	int i = _numrows;
	int j = _numcols;
	
	double MtoE = _logTransMatrix[M][E] + _forwardTable[M][i][j];	
	double ItoE = _logTransMatrix[M][E] + _forwardTable[I][i][j];	
	double DtoE = _logTransMatrix[M][E] + _forwardTable[D][i][j];	
	int tableindex = getNextTableIndex(MtoE, ItoE, DtoE, threshold);
	lls += _logTransMatrix[tableindex][E];

	while (i > 0 && j > 0) {
		if (tableindex == 0) {
			aln_target = ALPHS[_seqTarget->at(i-1)] + aln_target;
			for (int fi=0; fi < _familysize; fi++) 
				aln_family[fi] = ALPHS[_seqFamily->at(fi).at(j-1)] + aln_family[fi];

			double logdouble = getLogJointProb(i, j);
			double mlogp = logdouble + _logTransMatrix[M][M] + _forwardTable[M][i-1][j-1];
			double ilogp = logdouble + _logTransMatrix[I][M] + _forwardTable[I][i-1][j-1];
			double dlogp = logdouble + _logTransMatrix[D][M] + _forwardTable[D][i-1][j-1];
			tableindex = getNextTableIndex(mlogp, ilogp, dlogp, threshold);

			i--;
			j--;
			
			lls += logdouble;
			if (i > 0 && j > 0) lls += _logTransMatrix[tableindex][M];
			else lls += _logTransMatrix[S][M];
		} else if (tableindex == 1) {
			aln_target = ALPHS[_seqTarget->at(i-1)] + aln_target;
			for (int fi=0; fi < _familysize; fi++) 
				aln_family[fi] = "-" + aln_family[fi];

			double logsingle = _logBkgProbs[_seqTarget->at(i-1)];
			double mlogp = logsingle + _logTransMatrix[M][I] + _forwardTable[M][i-1][j];
			double ilogp = logsingle + _logTransMatrix[I][I] + _forwardTable[I][i-1][j];
			double dlogp = logsingle + _logTransMatrix[D][I] + _forwardTable[D][i-1][j];
			tableindex = getNextTableIndex(mlogp, ilogp, dlogp, threshold);

			i--;
			
			lls += logsingle;
			if (i > 0) lls += _logTransMatrix[tableindex][I];
			else lls += _logTransMatrix[S][I];
		} else {
			aln_target = "-" + aln_target;
			for (int fi=0; fi < _familysize; fi++) 
				aln_family[fi] = ALPHS[_seqFamily->at(fi).at(j-1)] + aln_family[fi];
			
			double logsingle = getLogJointProbSingle(j);
			double mlogp = logsingle + _logTransMatrix[M][D] + _forwardTable[M][i][j-1];
			double ilogp = logsingle + _logTransMatrix[I][D] + _forwardTable[I][i][j-1];
			double dlogp = logsingle + _logTransMatrix[D][D] + _forwardTable[D][i][j-1];
			tableindex = getNextTableIndex(mlogp, ilogp, dlogp, threshold);

			j--;
			
			lls += logsingle;
			if (j > 0) lls += _logTransMatrix[tableindex][D];
			else lls += _logTransMatrix[S][D];
		}
	} // end of while i & j

	while (i > 0) {
		double logsingle = _logBkgProbs[_seqTarget->at(i-1)];
		aln_target = ALPHS[_seqTarget->at(i-1)] + aln_target;
		for (int fi=0; fi < _familysize; fi++) 
			aln_family[fi] = "-" + aln_family[fi];
		i--;
		
			
		lls += logsingle;
		if (i > 0) lls += _logTransMatrix[I][I];
		else lls += _logTransMatrix[S][I];
	} // end of while i

	while (j > 0) {
		double logsingle = getLogJointProbSingle(j);
		aln_target = "-" + aln_target;
		for (int fi=0; fi < _familysize; fi++) 
			aln_family[fi] = ALPHS[_seqFamily->at(fi).at(j-1)] + aln_family[fi];

		j--;
			
		lls += logsingle;
		if (i > 0) lls += _logTransMatrix[D][D];
		else lls += _logTransMatrix[S][D];
	} // end of while j
	
	// make alignment string
	straln = "";
	int fi = 0;
	for (int si = 0; si < _familysize+1; si++) {
		string seqname = sequences.at(si).name();
		//straln += ">" + seqname	+ "\n";	
		ostringstream os;
        os << lls - _logfullprob;
		straln += ">" + seqname	+ "\t" + os.str() + "\n";	
		if (_targetIndex == si) {
			straln += aln_target + "\n";
		} else {
			straln += aln_family.at(fi) + "\n";
			fi++;
		}
	} // end of for si
	
	return (lls - _logfullprob);
}

double Model::getNextTableIndex(double mlogp, double ilogp, double dlogp, double threshold)
{
	int tableindex = 0;
	double sum = logsum(mlogp, ilogp);
	sum = logsum(sum, dlogp);
    double probM = exp(mlogp - sum);
    double probI = exp(ilogp - sum);
    double probX = exp(dlogp - sum);

	// fine max
	int maxindex = 0;
	double probMax = probM;
	if (probI > probM && probI > probX) {
		probMax = probI;
		maxindex = 1;
	} else if (probX > probM && probX > probI) {
		probMax = probX;
		maxindex = 2;
	}
	if (probMax >= threshold) {
		return maxindex; 
	}

    double rndnumber = rand()/(double)RAND_MAX;
    if (rndnumber < probM) tableindex = 0;
    else if (rndnumber >= probM && rndnumber < (probM+probI)) tableindex = 1;
    else tableindex = 2;

	return tableindex;
}



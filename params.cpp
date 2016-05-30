#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include "params.h"

Params::Params(const char* filename)
{
	ifstream infile;
	infile.open(filename);
	if (!infile) {
		cerr << "\n[ERROR]Unable to open file " << filename << endl;
		exit(1);
	} // end of if

	_MtoI = _MtoD = _ItoD = _DtoI = _ItoI = _DtoD = -1;

	string line;
	while (infile) {
		getline(infile, line);
		if (line.length() == 0 || line[0] == '#') { continue; }

		string key;	
		double value;
		int cnt = 0;
		istringstream iss(line);
		string token;
		while (getline(iss, token, '=')) {
			if (cnt == 0) key = token;
			else value = atof(token.c_str());
			cnt++;
		} // end of while 
		if (key.compare("SUBSTPAR") == 0)  _numsubst = value; 
		else if (key.compare("PA") == 0) _pA = value;
		else if (key.compare("PC") == 0) _pC = value;
		else if (key.compare("PG") == 0) _pG = value;
		else if (key.compare("PT") == 0) _pT = value;
		else if (key.compare("MtoIS") == 0) _MtoI = value;
		else if (key.compare("MtoIA") == 0) _MtoD = value;
		else if (key.compare("IStoIS") == 0) _ItoI = value;
		else if (key.compare("IAtoIA") == 0) _DtoD = value;
		else if (key.compare("IAtoIS") == 0) _DtoI = value;
		else if (key.compare("IStoIA") == 0) _ItoD = value;
		else {
			cerr << "\n[ERROR]Unrecognized parameter entry in the file " << filename << endl;
			exit(1);
		}	
	} // end of while
}

double Params::getNumSubst() 
{
	return _numsubst; 
}

double Params::getProbA()
{
	return _pA;
}

double Params::getProbC()
{
	return _pC;
}

double Params::getProbG()
{
	return _pG;
}

double Params::getProbT()
{
	return _pT;
}

double Params::getMtoI()
{
	return _MtoI;
}

double Params::getMtoD()
{
	return _MtoD;
}

double Params::getItoI()
{
	return _ItoI;
}

double Params::getItoD()
{
	return _ItoD;
}

double Params::getDtoD()
{
	return _DtoD;
}

double Params::getDtoI()
{
	return _DtoI;
}

bool Params::hasTransProb()
{
	if (_MtoI > 0.0 && _MtoD > 0.0 && _ItoI > 0.0 && _ItoD > 0.0 && _DtoD > 0.0 && _DtoI > 0.0) return true;
	return false;
}


#ifndef PARAMS_H
#define PARAMS_H

using namespace std;

class Params 
{
private:
	double _numsubst;
	double _pA, _pC, _pG, _pT;
	double _MtoI, _MtoD, _ItoD, _DtoI, _ItoI, _DtoD;

public:
	Params(const char* filename);
	double getNumSubst();
	double getProbA();
	double getProbC();
	double getProbG();
	double getProbT();
	double getMtoI();
	double getMtoD();
	double getItoI();
	double getDtoD();
	double getItoD();
	double getDtoI();
	bool hasTransProb();
};

#endif

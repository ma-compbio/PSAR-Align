#ifndef SCORES_H
#define SCORES_H

#include <vector>
#include <map>
#include "sequence.h"

using namespace std;

class ResiduePair
{
public: 
	int seq1, seq2;
	int pos1, pos2;

	ResiduePair(int _seq1, int _seq2, int _pos1, int _pos2): seq1(_seq1), seq2(_seq2), pos1(_pos1), pos2(_pos2) {}

	bool operator< (const ResiduePair& rp) const {
		if (seq1 == rp.seq1 && seq2 == rp.seq2 && pos1 == rp.pos1) return pos2 < rp.pos2;
		if (seq1 == rp.seq1 && seq2 == rp.seq2) return pos1 < rp.pos1;
		if (seq1 == rp.seq1) return seq2 < rp.seq2;
		return seq1 < rp.seq1;
	}
};

class Scores
{
private:
	map<ResiduePair, double> _pairscores;	
	map<ResiduePair, double>::iterator _pairiter;	
	
	map<int, vector<ResiduePair> > _columns;
	map<int, vector<ResiduePair> >::iterator _coliter;

public:
	Scores(vector<Sequence>& sequences);
	void compute_scores (string dir);
	void writetofile (string fnamepair, string fnamecol);
};

#endif

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <string>

using namespace std;
   
enum {A,C,G,T,N,GAP};
const string ALPHS[] = {"A","C","G","T","N","-"};

class Sequence 
{
private:
	string _name;
	int _length;
	vector<int> _sequence;

public:
	Sequence(string name, string sequence);	
	// copy constructor
	Sequence(const Sequence& seq);

	Sequence& operator=(const Sequence& seq);

	string name();	
	int length();
	int at(int i);
	void print();
	vector<int> removeGaps();
	void removeGapAt(int i);
	int getGapSize();

	static vector<Sequence>& removeNullColumns(vector<Sequence>& seqs);
	static int getNumofGappedColumns(vector<Sequence>& seqs);
	static void read_alignments(const char* filename, vector<Sequence>& sequences);
};

#endif

#include <iostream>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include "sequence.h"

Sequence::Sequence(string name, string sequence)
{
	_name = name;
	_length = sequence.length();
	for (int i = 0; i < _length; i++) {
		char ch = toupper(sequence.at(i));
		int chnum = 0;
		if (ch == 'A') chnum = 0; 
		else if (ch == 'C') chnum = 1; 
		else if (ch == 'G') chnum = 2; 
		else if (ch == 'T') chnum = 3; 
		else if (ch == 'N') chnum = 4; 
		else if (ch == '-') chnum = 5; 
		else {
			cerr << "[ERROR] Unsupported nucleotide symbol " << ch << endl;
			exit(1);
		} // end of else
		_sequence.push_back(chnum);
	}
}

Sequence::Sequence(const Sequence& seq)
{
	_name = seq._name;
	_length = seq._length;
	_sequence = seq._sequence;
}

Sequence& Sequence::operator=(const Sequence& seq)
{
	if (this == &seq) return *this;

	_name = seq._name;
	_length = seq._length;
	_sequence = seq._sequence;

	return *this;
}

string Sequence::name()
{
	return _name;
}

int Sequence::length()
{
	return _length;
}

int Sequence::at(int i)
{
	return _sequence.at(i);
}

void Sequence::print()
{
	cout << ">" << _name << " len=" << _length << endl;
	string sequence = "";
	for (int i = 0; i < _length; i++) sequence += ALPHS[_sequence.at(i)];
	cout << sequence << endl;
}

vector<int> Sequence::removeGaps()
{
	int gapsize = 0;
	int gapcnt = 0;
	int arraysize = _length;
	int i = 0;
	bool ingap = false;
	while (i < arraysize) {
		if (_sequence.at(i) == GAP) {
			gapsize++;
			if (ingap == false) {
				ingap = true;
				gapcnt++;	
			}
			_sequence.erase(_sequence.begin() + i);
			arraysize--;
		} else {
			i++;
			ingap = false;
		}
	} // end of while
	_length = i;

	vector<int> gapstats;
	gapstats.push_back(gapsize);
	gapstats.push_back(gapcnt);
	return gapstats;	
}

void Sequence::removeGapAt(int i)
{
	_sequence.erase(_sequence.begin() + i);
	_length--;
}

vector<Sequence>& Sequence::removeNullColumns(vector<Sequence>& seqs)
{
	int numseq = seqs.size();
	int seqlen = seqs.at(0).length();

	int i = 0;
	while (i < seqlen) {
		// count the number of gaps
		int gapcnt = 0;
		for (int j = 0; j < numseq; j++) {
			if (seqs.at(j).at(i) == GAP) gapcnt++;
			else break;
		} // end of for j

		if (gapcnt == numseq) {
			for (int j = 0; j < numseq; j++) seqs.at(j).removeGapAt(i);
			seqlen--;
		} else i++;
	}

	return seqs;
}

int Sequence::getGapSize() 
{
	int cnt = 0;
	for (int i = 0; i < _length; i++) {
		if (_sequence.at(i) == GAP) cnt++;
	}
	return cnt;
}

int Sequence::getNumofGappedColumns(vector<Sequence>& seqs)
{
	int numseq = seqs.size();
	int seqlen = seqs.at(0).length();

	int cnt = 0;
	for (int i = 0; i < seqlen; i++) {	
		for (int j = 0; j < numseq; j++) {
			if (seqs.at(j).at(i) == GAP) { cnt++; break; }
		} // end of for j
	}

	return cnt;
}

void Sequence::read_alignments(const char* filename, vector<Sequence>& sequences)
{
    ifstream infile;
    infile.open(filename);
    if (!infile) {
        cerr << "\n[ERROR]Unable to open file " << filename << endl;
        exit(1);
    } // end of if

    string line;
    string spc = "";
    string seq = "";
    while (infile) {
        getline(infile, line);
        if (line.length() == 0) { continue; }
        if (line[0] == '>') {
            if (seq.length() > 0) {
                sequences.push_back(Sequence(spc, seq));
            } // end of if  

            spc = line.substr(1);
            seq = "";
        } else {
            seq += line;
        } // end of if else
    } // end of while
    infile.close();
    sequences.push_back(Sequence(spc, seq));
}


/* Code generated by HMMoC version 1.3, Copyright (C) 2006 Gerton Lunter */
/* Generated from file nucleotide.xml (author:  Robert K. Bradley ) on Tue Dec 23 01:04:16 CST 2008 */

/*
This file is a work based on HMMoC 1.3, a hidden Markov model compiler.
Copyright (C) 2006 by Gerton Lunter, Oxford University.

HMMoC and works based on it are free software; you can redistribute 
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

HMMOC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HMMoC; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef _nucleotidedp_h_
#define _nucleotidedp_h_


#include "dptables.h"
#include "algebras.h"
#include <string>

#include <map>

using std::map;


// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesnucBlock2;
typedef States<bfloat,1> StatesnucBlock1;
typedef States<bfloat,1> StatesnucBlock3;

class NucleotideAlignDPTable {
    public:
    // If true, this class' destructor will delete the DP arrays
    bool isInCharge;
    // Pointers to arrays containing ids of states and transitions
    const string* const stateId;
    const string* const emissionId;
    const string* const transitionId;
    const string* const transitionFrom;
    const string* const transitionTo;
    const string* const transitionProb;
    const string* const transitionEmit;
    const string* const outputId;
    // The actual DP tables, and total sequence lengths (which determine size of DP arrays) follow:
    int iLen1;
    int iLen2;
    DPTable<StatesnucBlock2,2> StateMemorynucBlock2;
    DPTable<StatesnucBlock1,0> StateMemorynucBlock1;
    DPTable<StatesnucBlock3,0> StateMemorynucBlock3;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    NucleotideAlignDPTable(int iLen1,int iLen2);
    ~NucleotideAlignDPTable();
    // returns probability from DP table, given position and int or string state identifier
    bfloat getProb(int iState ,int ,int ) const;
    bfloat getProb(const string sState ,int ,int ) const;
    // converts string identifier (for state, transition or emission) into integer id
    static int getId(const string& sState);
    static const string& getTransitionId(int id);
    static const string& getEmissionId(int id);
    static const string& getStateId(int id);
    static const string& getOutputId(int id);
    static void _cleanup() { getId("_cleanup_"); }
};

// give a name to the real type used for this HMM
typedef bfloat NucleotideAlignReal;
// define type for a 'short' real -- usually double, but can be logspace for efficiency
typedef double NucleotideAlignShortReal;



// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesnucBlock2;
typedef States<bfloat,1> StatesnucBlock3;
typedef States<bfloat,1> StatesnucBlock1;

class NucleotideAlignFoldedDPTable {
    public:
    // If true, this class' destructor will delete the DP arrays
    bool isInCharge;
    // Pointers to arrays containing ids of states and transitions
    const string* const stateId;
    const string* const emissionId;
    const string* const transitionId;
    const string* const transitionFrom;
    const string* const transitionTo;
    const string* const transitionProb;
    const string* const transitionEmit;
    const string* const outputId;
    // The actual DP tables, and total sequence lengths (which determine size of DP arrays) follow:
    int iLen1;
    int iLen2;
    FoldedTable<DPTable,StatesnucBlock2,2> StateMemorynucBlock2;
    DPTable<StatesnucBlock3,0> StateMemorynucBlock3;
    DPTable<StatesnucBlock1,0> StateMemorynucBlock1;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    NucleotideAlignFoldedDPTable(int iLen1,int iLen2);
    ~NucleotideAlignFoldedDPTable();
    // returns probability from DP table, given position and int or string state identifier
    bfloat getProb(int iState ,int ,int ) const;
    bfloat getProb(const string sState ,int ,int ) const;
    // converts string identifier (for state, transition or emission) into integer id
    static int getId(const string& sState);
    static const string& getTransitionId(int id);
    static const string& getEmissionId(int id);
    static const string& getStateId(int id);
    static const string& getOutputId(int id);
    static void _cleanup() { getId("_cleanup_"); }
};



class NucleotideAlignBaumWelch {
    public:
    // Default copy constructor is used.
    // Void constructor:
    NucleotideAlignBaumWelch() { resetCounts(); }
    // Not calling resetCounts() across calls allows to aggregate results over multiple datasets
    void resetCounts();
    void scaleCounts(bfloat scale);
    // Translate an identifier (string or integer) to the index into their corresponding Baum-Welch counter array (below)
    // Which array is used for any particular emission/transition depends on its order signature - see documentation for details
    int transitionIndex(int intId) const { return atransitionIdx[intId]; }
    int transitionIndex(string strId) const;
    int emissionIndex(int intId) const { return aemissionIdx[intId]; }
    int emissionIndex(string strId) const;
    // Now follow, in triplets (one for each order signature):
    //  Transition or emission counters;
    //  Array of identifiers; and
    //  Dimension of array (number of counters).
    bfloat transitionBaumWelchCount00[13];
    static int transitionIdentifier00[13];   
    static const int transitionDimension00 = 13;
    bfloat emissionBaumWelchCount00[1];
    static int emissionIdentifier00[1];   
    static const int emissionDimension00 = 1;
    bfloat emissionBaumWelchCount01[4][1];
    static int emissionIdentifier01[1];   
    static const int emissionDimension01 = 1;
    bfloat emissionBaumWelchCount10[4][1];
    static int emissionIdentifier10[1];   
    static const int emissionDimension10 = 1;
    bfloat emissionBaumWelchCount11[4][4][1];
    static int emissionIdentifier11[1];   
    static const int emissionDimension11 = 1;
    private:
    static int atransitionIdx[13];
    static int aemissionIdx[4];
    static map<const string,int> mId;
};




// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesnucBlock2withbanding;
typedef States<bfloat,1> StatesnucBlock1;
typedef States<bfloat,1> StatesnucBlock3;

class NucleotideAlignWithBandingDPTable {
    public:
    // If true, this class' destructor will delete the DP arrays
    bool isInCharge;
    // Pointers to arrays containing ids of states and transitions
    const string* const stateId;
    const string* const emissionId;
    const string* const transitionId;
    const string* const transitionFrom;
    const string* const transitionTo;
    const string* const transitionProb;
    const string* const transitionEmit;
    const string* const outputId;
    // The actual DP tables, and total sequence lengths (which determine size of DP arrays) follow:
    int iLen1;
    int iLen2;
    DPTable<StatesnucBlock2withbanding,2> StateMemorynucBlock2withbanding;
    DPTable<StatesnucBlock1,0> StateMemorynucBlock1;
    DPTable<StatesnucBlock3,0> StateMemorynucBlock3;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    NucleotideAlignWithBandingDPTable(int iLen1,int iLen2);
    ~NucleotideAlignWithBandingDPTable();
    // returns probability from DP table, given position and int or string state identifier
    bfloat getProb(int iState ,int ,int ) const;
    bfloat getProb(const string sState ,int ,int ) const;
    // converts string identifier (for state, transition or emission) into integer id
    static int getId(const string& sState);
    static const string& getTransitionId(int id);
    static const string& getEmissionId(int id);
    static const string& getStateId(int id);
    static const string& getOutputId(int id);
    static void _cleanup() { getId("_cleanup_"); }
};

// give a name to the real type used for this HMM
typedef bfloat NucleotideAlignWithBandingReal;
// define type for a 'short' real -- usually double, but can be logspace for efficiency
typedef double NucleotideAlignWithBandingShortReal;



// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesnucBlock2withbanding;
typedef States<bfloat,1> StatesnucBlock3;
typedef States<bfloat,1> StatesnucBlock1;

class NucleotideAlignWithBandingFoldedDPTable {
    public:
    // If true, this class' destructor will delete the DP arrays
    bool isInCharge;
    // Pointers to arrays containing ids of states and transitions
    const string* const stateId;
    const string* const emissionId;
    const string* const transitionId;
    const string* const transitionFrom;
    const string* const transitionTo;
    const string* const transitionProb;
    const string* const transitionEmit;
    const string* const outputId;
    // The actual DP tables, and total sequence lengths (which determine size of DP arrays) follow:
    int iLen1;
    int iLen2;
    FoldedTable<DPTable,StatesnucBlock2withbanding,2> StateMemorynucBlock2withbanding;
    DPTable<StatesnucBlock3,0> StateMemorynucBlock3;
    DPTable<StatesnucBlock1,0> StateMemorynucBlock1;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    NucleotideAlignWithBandingFoldedDPTable(int iLen1,int iLen2);
    ~NucleotideAlignWithBandingFoldedDPTable();
    // returns probability from DP table, given position and int or string state identifier
    bfloat getProb(int iState ,int ,int ) const;
    bfloat getProb(const string sState ,int ,int ) const;
    // converts string identifier (for state, transition or emission) into integer id
    static int getId(const string& sState);
    static const string& getTransitionId(int id);
    static const string& getEmissionId(int id);
    static const string& getStateId(int id);
    static const string& getOutputId(int id);
    static void _cleanup() { getId("_cleanup_"); }
};



class NucleotideAlignWithBandingBaumWelch {
    public:
    // Default copy constructor is used.
    // Void constructor:
    NucleotideAlignWithBandingBaumWelch() { resetCounts(); }
    // Not calling resetCounts() across calls allows to aggregate results over multiple datasets
    void resetCounts();
    void scaleCounts(bfloat scale);
    // Translate an identifier (string or integer) to the index into their corresponding Baum-Welch counter array (below)
    // Which array is used for any particular emission/transition depends on its order signature - see documentation for details
    int transitionIndex(int intId) const { return atransitionIdx[intId]; }
    int transitionIndex(string strId) const;
    int emissionIndex(int intId) const { return aemissionIdx[intId]; }
    int emissionIndex(string strId) const;
    // Now follow, in triplets (one for each order signature):
    //  Transition or emission counters;
    //  Array of identifiers; and
    //  Dimension of array (number of counters).
    bfloat transitionBaumWelchCount00[13];
    static int transitionIdentifier00[13];   
    static const int transitionDimension00 = 13;
    bfloat emissionBaumWelchCount00[1];
    static int emissionIdentifier00[1];   
    static const int emissionDimension00 = 1;
    bfloat emissionBaumWelchCount01[4][1];
    static int emissionIdentifier01[1];   
    static const int emissionDimension01 = 1;
    bfloat emissionBaumWelchCount10[4][1];
    static int emissionIdentifier10[1];   
    static const int emissionDimension10 = 1;
    bfloat emissionBaumWelchCount11[4][4][1];
    static int emissionIdentifier11[1];   
    static const int emissionDimension11 = 1;
    private:
    static int atransitionIdx[13];
    static int aemissionIdx[4];
    static map<const string,int> mId;
};



bfloat Forward(NucleotideAlignDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT);

bfloat BackwardBaumWelch(NucleotideAlignBaumWelch& bw,NucleotideAlignDPTable* pInTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT);

bfloat Backward(NucleotideAlignDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT);

bfloat ForwardBanding(NucleotideAlignWithBandingDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT,int iWidth);

bfloat BackwardBaumWelchBanding(NucleotideAlignWithBandingBaumWelch& bw,NucleotideAlignWithBandingDPTable* pInTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT,int iWidth);

bfloat BackwardBanding(NucleotideAlignWithBandingDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT,int iWidth);

#endif // _nucleotidedp_h_

/* --- end of HMMoC-generated file --- */
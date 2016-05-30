/* Code generated by HMMoC version 1.3, Copyright (C) 2006 Gerton Lunter */
/* Generated from file aminoacid.xml (author:  Robert K. Bradley ) on Tue Dec 23 01:04:18 CST 2008 */

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

#ifndef _aminoaciddp_h_
#define _aminoaciddp_h_


#include "dptables.h"
#include "algebras.h"
#include <string>

#include <map>

using std::map;


// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesaaBlock2;
typedef States<bfloat,1> StatesaaBlock1;
typedef States<bfloat,1> StatesaaBlock3;

class AminoAcidAlignDPTable {
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
    DPTable<StatesaaBlock2,2> StateMemoryaaBlock2;
    DPTable<StatesaaBlock1,0> StateMemoryaaBlock1;
    DPTable<StatesaaBlock3,0> StateMemoryaaBlock3;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    AminoAcidAlignDPTable(int iLen1,int iLen2);
    ~AminoAcidAlignDPTable();
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
typedef bfloat AminoAcidAlignReal;
// define type for a 'short' real -- usually double, but can be logspace for efficiency
typedef double AminoAcidAlignShortReal;



// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesaaBlock2;
typedef States<bfloat,1> StatesaaBlock3;
typedef States<bfloat,1> StatesaaBlock1;

class AminoAcidAlignFoldedDPTable {
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
    FoldedTable<DPTable,StatesaaBlock2,2> StateMemoryaaBlock2;
    DPTable<StatesaaBlock3,0> StateMemoryaaBlock3;
    DPTable<StatesaaBlock1,0> StateMemoryaaBlock1;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    AminoAcidAlignFoldedDPTable(int iLen1,int iLen2);
    ~AminoAcidAlignFoldedDPTable();
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



class AminoAcidAlignBaumWelch {
    public:
    // Default copy constructor is used.
    // Void constructor:
    AminoAcidAlignBaumWelch() { resetCounts(); }
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
    bfloat emissionBaumWelchCount01[20][1];
    static int emissionIdentifier01[1];   
    static const int emissionDimension01 = 1;
    bfloat emissionBaumWelchCount10[20][1];
    static int emissionIdentifier10[1];   
    static const int emissionDimension10 = 1;
    bfloat emissionBaumWelchCount11[20][20][1];
    static int emissionIdentifier11[1];   
    static const int emissionDimension11 = 1;
    private:
    static int atransitionIdx[13];
    static int aemissionIdx[4];
    static map<const string,int> mId;
};




// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesaaBlock2withbanding;
typedef States<bfloat,1> StatesaaBlock1;
typedef States<bfloat,1> StatesaaBlock3;

class AminoAcidAlignWithBandingDPTable {
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
    DPTable<StatesaaBlock2withbanding,2> StateMemoryaaBlock2withbanding;
    DPTable<StatesaaBlock1,0> StateMemoryaaBlock1;
    DPTable<StatesaaBlock3,0> StateMemoryaaBlock3;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    AminoAcidAlignWithBandingDPTable(int iLen1,int iLen2);
    ~AminoAcidAlignWithBandingDPTable();
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
typedef bfloat AminoAcidAlignWithBandingReal;
// define type for a 'short' real -- usually double, but can be logspace for efficiency
typedef double AminoAcidAlignWithBandingShortReal;



// Here go the state memory clique typedefs:
typedef States<bfloat,3> StatesaaBlock2withbanding;
typedef States<bfloat,1> StatesaaBlock3;
typedef States<bfloat,1> StatesaaBlock1;

class AminoAcidAlignWithBandingFoldedDPTable {
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
    FoldedTable<DPTable,StatesaaBlock2withbanding,2> StateMemoryaaBlock2withbanding;
    DPTable<StatesaaBlock3,0> StateMemoryaaBlock3;
    DPTable<StatesaaBlock1,0> StateMemoryaaBlock1;
    // Member functions:
    public:
    // Default copy constructor is used; user has to set isInCharge appropriately afterwards!
    AminoAcidAlignWithBandingFoldedDPTable(int iLen1,int iLen2);
    ~AminoAcidAlignWithBandingFoldedDPTable();
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



class AminoAcidAlignWithBandingBaumWelch {
    public:
    // Default copy constructor is used.
    // Void constructor:
    AminoAcidAlignWithBandingBaumWelch() { resetCounts(); }
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
    bfloat emissionBaumWelchCount01[20][1];
    static int emissionIdentifier01[1];   
    static const int emissionDimension01 = 1;
    bfloat emissionBaumWelchCount10[20][1];
    static int emissionIdentifier10[1];   
    static const int emissionDimension10 = 1;
    bfloat emissionBaumWelchCount11[20][20][1];
    static int emissionIdentifier11[1];   
    static const int emissionDimension11 = 1;
    private:
    static int atransitionIdx[13];
    static int aemissionIdx[4];
    static map<const string,int> mId;
};



bfloat Forward(AminoAcidAlignDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT);

bfloat BackwardBaumWelch(AminoAcidAlignBaumWelch& bw,AminoAcidAlignDPTable* pInTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT);

bfloat Backward(AminoAcidAlignDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT);

bfloat ForwardBanding(AminoAcidAlignWithBandingDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT,int iWidth);

bfloat BackwardBaumWelchBanding(AminoAcidAlignWithBandingBaumWelch& bw,AminoAcidAlignWithBandingDPTable* pInTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT,int iWidth);

bfloat BackwardBanding(AminoAcidAlignWithBandingDPTable** ppOutTable,const std::string& iSequence1,const std::string& iSequence2,const vector<double>& iSingleDistribution,const vector<vector<double> >& iPairDistribution,const vector<vector<double> >& iT,int iWidth);

#endif // _aminoaciddp_h_

/* --- end of HMMoC-generated file --- */

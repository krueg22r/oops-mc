//The ForceField is in charge of holding potential objects
//and mediates everything that requires force field parameters

#ifndef FORCEFIELD_H
#define FORCEFIELD_H
#include <string>
#include <vector>
#include "Bead.h"
#include "PairPot.h"
#include "BondPot.h" 
#include "ExtPot.h" 
#include "Spring.h" 
#include "CoulLjCut.h"
#include "LjCut.h" 
#include "LjWall2.h" 
#include <fstream>
#include <iomanip>
#include <random> 
using namespace std;

class ForceField{
	//basic run params
	double m_beta; 
	int m_dimension; 
	double m_deltaLen; 
	//data structs for virtual volume changes
        double pPlus[3];
        double pMinus[3];
	bool extDef; //has an external potential been defined? 
	//grand canonical params
	bool m_gc; 
	string m_chainSymbol; 
	int m_chainLen; 
	int m_chainCharge; 
	int m_chargeFrac; 
	double m_qT; 
	double m_chemPot; 
	int m_gcEvery; 
	int m_nVec;
	//data structures for gc objects
	vector < Bead > trialVecs; 
	vector < Bead > trialChain;  
	double* wFactors; //dynamic array will have size m_nVec; 
	//potential objects
	PairPot* pairPot; //first two are non-optional 
	BondPot* bondPot;
	ExtPot* extPot; 
	//variables to track gc move success
	int delAccept; 
	int insAccept; 	
	int delAttempt; 
	int insAttempt; 
	int enCounter; //count how many times en has been calculated
	double totDense; //running total density value
	double totForce; //running total force value
	public: 
	double length[3]; 
	ForceField(double, double, int, int); 
	void readParam(int); 
	void initFF(vector < Molecule >&); 
	double processMove(vector < Molecule >&, int, int); 
	void resetEnergies(vector < Molecule >&, bool, int); 
	void printHeader(ofstream&); 
	void printEnergies(vector < Molecule >&, ofstream&, int, int, int); 
	//gc funcs 
	double beadSysEn(Bead&, vector < Molecule >&, int); 
	double genVecs(Bead&, vector < Molecule >&, double, int, mt19937&); 
	bool trialIns(vector < Molecule >&, double, mt19937&); 
	int trialDel(vector < Molecule >&, int, double, mt19937&);
	void addMol(vector < Molecule >&); 
	//for virtual volume changes
	double getDeltaV(int);
        void calcPress(vector < Molecule >&, ofstream&, int, double);
	void enStats(ofstream&); 
	int gcEvery(); 
	bool gc(); 
}; 


#endif

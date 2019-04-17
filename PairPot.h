//Base object for pair potential. This must be defined for every simulation. 
//To avoid recalculating after each move (very expensive)
//exisiting energies are kept in a map with keys consisting of a pair containing 
//the two bead IDs. 
#ifndef PAIRPOT_H
#define PAIRPOT_H
#include <string> 
#include <vector> 
#include <map> 
#include <iostream> 
#include <numeric> 
#include <iomanip> 
#include "Molecule.h" 
#include "Bead.h"
#include "Small.h"
using namespace std; 

class PairPot{
      	map<pair<int,int>, double> eOldMap; 
	map<pair<int,int>, double> eTrialMap;   
	double m_totOld; 
	double m_deTrial; 
	public: 
	PairPot(); //default constructor 
        PairPot(int);
	void setVal(int, int, int, double);
	void setAllVal(int, int, double); 
	double getVal(int, int, int);
	double deltaEnMove(vector < Molecule >&, int, double[], int);
	void resetMaps(vector < Molecule >&, int, bool); 
	void delMol(vector < Molecule >&, int); 
	virtual void readParam() = 0; 
	//virtual string name() = 0; 
        void initialize(vector < Molecule >&, double[], int); 
	void initMol(vector < Molecule >&, double[], int); 
        double totEn(); 
        virtual double pairEn(Bead&, Bead&, double[], int) = 0;
	double molMolEn(Molecule&, Molecule&, double, int); 
	double intermolEn(Molecule&, double, int); 
	virtual void printHeader(ostream&) = 0; 
	void printEn(ostream&); 
	double allTrialEn(vector < Molecule >&, double[], int);
}; 


#endif

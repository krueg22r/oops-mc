//Base class for external potentials. 
//"external potential" is any pot where the energy for a bead
//depends only on the bead's own coordinates. 
//useful for hard walls, electric fields, etc. 
//the extpot is optional; only one may be defined per simulation
#ifndef EXTPOT_H
#define EXTPOT_H

#include <vector>
#include <map> 
#include <iostream>
#include <iomanip> 
#include <numeric> 
#include "Bead.h"
#include "Molecule.h" 
using namespace std;

class ExtPot{
	map<int, double> eOldMap; 
	map<int, double> eTrialMap; 
	double m_totOld;
        double m_deTrial;
	public: 
	ExtPot(); 
	void setVal(int, int, double); 
	double getVal(int, int); 
	void setAllVal(int, double); 
	void initialize(vector < Molecule >& mols, double[], int); 
	void initMol(vector < Molecule >& mols, double[], int);
	void delMol(vector < Molecule >& mols, int); 
	double deltaEnMove(vector < Molecule >& mols, int, double[], int);	
	void resetMaps(vector < Molecule >& mols, int, bool); 
	double totEn(); 
	double allTrialEn(vector < Molecule >& mols, double[], int);
	void printEn(ostream&); 
	virtual void readParam() = 0; 
	virtual double beadEn(Bead&, double[]) = 0; 
	virtual double calcForce(vector < Molecule >&, double[]) = 0;
	virtual void printHeader(ostream&) = 0;
}; 

#endif

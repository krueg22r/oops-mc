//A base class for bonded potentials
//Derived classes may contain any combination of bond, angle, and dihedral energy funcs
//May be any potential that depends only upon the positions of beads within a single molecule
//Exactly one BondPot must be defined in a given simulation, even if no bonds are listed in the 
//topology file 

#ifndef BONDPOT_H 
#define BONDPOT_H

#include <vector>
#include <iostream>
#include <iomanip> 
#include <numeric> 
#include <random> 
#include "Molecule.h" 
using namespace std; 



class BondPot{
	vector < double > arr_eOld; //stores energies in vectors 
        vector < double > arr_eTrial; //to avoid necessity of recalculating 
	protected: 
	double m_totOld; //total "accepted" energy
        double m_deTrial; //total trial energy 
	public:
        BondPot(int); //arg: initial # molecules 
	void initialize(vector < Molecule >&, double[], int); //args are mol array, length, dimension
	void initMol(vector < Molecule >&, double[], int); //init one new mol being added 
	void delMol(int); //delete stored energies for a molecule being deleted 
	virtual void printHeader(ostream&) = 0; 
	void setVal(int, int, double); //same array, index, value args 
        double getVal(int, int);
	void resetEn(int, bool); 
	double totEn(); 
	virtual double enMol(Molecule&, double[], int) = 0; 
	virtual double deltaEnMove(vector < Molecule >&, double[], int, int, int) = 0; 
	virtual void readParam() = 0;
	virtual double ranBond(double, mt19937&) = 0; //for growing mols for gc insertion
	double allTrialEn(vector < Molecule >&, double[], int); //for use in pressure calculation 
	void printEn(ostream&); 	
}; 

#endif

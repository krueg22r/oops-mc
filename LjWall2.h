//One possible external potential 
//places two walls perpendicular to z-axis
//one at 0, one at zLength. 
//each wall has the same LJ parameters and 
//cutoff. potentials truncated at cutoff distance 
//from the relevant wall. mixing rules same as for 
//LJ pair potential. 
#ifndef LJWALL2_H
#define LJWALL2_H

#include <iostream> 
#include <iomanip> 
#include <numeric> 
#include <string> 
#include <map> 
#include <vector> 
#include "Bead.h"
#include "Molecule.h" 
#include "ExtPot.h" 

using namespace std; 

class LjWall2 : public ExtPot{
	double m_cut; 
	map < string, double> sigmas; 
	map < string, double> epsilons; 
	double m_sigWall; 
	double m_epWall; 
	public: 
	LjWall2(); 
	void readParam(); 
	double beadEn(Bead&, double[]);
	double calcForce(vector < Molecule >&, double[]); 
	void printHeader(ostream&); 
}; 

#endif

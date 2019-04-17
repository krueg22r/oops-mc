//A simple bonded potential with only harmonic bonds between specified atoms
//no angle or dihedral potentials 
//currently the only BondPot supported for gc simulations 
#ifndef SPRING_H
#define SPRING_H
#include <stdlib.h> 
#include "BondPot.h" 
#include "Molecule.h" 
#include "Bead.h" 
#include "Small.h" 
#include <iomanip> 
#include <iostream> 
#include <random> 

class Spring: public BondPot{
	double m_kBond;
        double m_r0;
        public:
	Spring(int); 
	void readParam();
	double enMol(Molecule&, double[], int);
	double ranBond(double, mt19937&); 
	double deltaEnMove(vector < Molecule >&, double[], int, int, int);
	void printHeader(ostream&); 
}; 

#endif


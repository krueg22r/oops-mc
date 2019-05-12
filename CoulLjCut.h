//a pair potential that includes standard LJ/Coulomb interactions
//each one being truncated at the given cutoffs (may be different 
//for LJ and Coul). Each element gets its own LJ params. 
//we use arithmetic mixing for sigmas, geometric mixing for epsilons
#ifndef COULLJCUT_H
#define COULLJCUT_H

#include <string> 
#include <vector> 
#include <iostream> 
#include <iomanip> 
#include <map> 
#include "Bead.h" 
#include "Molecule.h" 
#include "PairPot.h" 
using namespace std; 

class CoulLjCut : public PairPot{
	double m_ljCut; 
	double m_coulCut; 
	double m_dielectric; // This just scales the interaction strength 
	map<string, double> sigmas; //potential keeps track of LJ params
	map<string, double> epsilons; //in constrast to charges, which are bead attributes
	public: 
	CoulLjCut(); 
	CoulLjCut(int); 
	void readParam(); 
	// Calculates energy between two beads 
	double pairEn(Bead&, Bead&, double[], int); 	
	void printHeader(ostream&); 
}; 

#endif

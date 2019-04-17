//Basic LJ interaction pair potential 
//Potential is truncated at given cutoff 
//Uses arithmetic mixing for sigmas, geometric 
//mixing for epsilons 

#ifndef LJCUT_H
#define LJCUT_H

#include <string> 
#include <vector> 
#include <iostream> 
#include <iomanip> 
#include <map> 
#include "Bead.h" 
#include "Molecule.h" 
#include "PairPot.h" 
using namespace std; 

class LjCut : public PairPot{
	double m_ljCut; 
	map<string, double> sigmas; 
	map<string, double> epsilons; 
	public: 
	LjCut(); 
	LjCut(int); 
	void readParam(); 
	double pairEn(Bead&, Bead&, double[], int); 	
	void printHeader(ostream&); 
}; 

#endif

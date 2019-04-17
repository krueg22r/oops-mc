//The SimBox object contains the molecule vector
//as well as info about how many things are in the box 
//and general MC params related to moving mols (beta, delta, etc) 
#ifndef SIMBOX_H
#define SIMBOX_H
#include "ForceField.h"
#include <stdlib.h>
#include "Molecule.h"
#include "Bead.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <set> 
#include <vector> 
#include <string>
#include <random> 
#include <chrono> 

using namespace std; 

class SimBox{
	int m_nTot; //total number of BEADS
	int m_nMol; //total number Molecules
	int m_nChain; //total number CHAIN Molecules 
	int m_rgCounter; 
	double m_rgTot; 
	double m_delta; 
	double m_beta; 
	int m_every; 
	int m_equilLen; 
	int m_idCounter; 
	vector < Molecule > mols;
	set<int> boxIDs;  
	double moveProbs[4]; 
	int accepted[5];
        int attempted[5];
	double length[3]; 
	mt19937 ranGen; 
	public: 
	//constructor with some basic info about species
	SimBox(double, double, int); 
	void readCoords(ifstream&); 
	void printCoords(ofstream&, int); //prints to traj file in xyz format
	void printLastFrame(ofstream&); //prints in my input format
	void printFinalTop(ofstream&); 
	void readTop(ifstream&); 
	void setLength(ForceField&); 
	void initEnergy(ForceField&); 
	double radGy(); 
	void printEnergy(ForceField&, ofstream&, int);
	int newID(); 
	void makeMove(ForceField&); 
	void gcMove(ForceField&);
	void runStats(); 
	int nMol(); 
	int nTot(); 
}; 
#endif

//A molecule contains a vector of monomer beads
//For single particles, this vector just contains one bead. 
//It also contains 2-d arrays giving indices of beads involved in 
//bonds/angles/dihedrals
//The molecule is in charge of bead positions, so it has functions to do 
//translational MC moves, scale coordinates for pressure calculation, and 
//print bead positions 

#ifndef MOLECULE_H
#define MOLECULE_H
#include "Bead.h"
#include <cmath> 
#include <vector> 
#include <iomanip> 
#include <random> 
using namespace std; 

class Molecule{
	//Bead* m_monomers; 
	int m_size; 
	int m_nBond; 
	int m_nAngle; 
	int m_nDihedral; 
	bool m_moved; 
	public: 
	Molecule(); 
	Molecule(const Molecule&); 
	vector < Bead > monomers; 	
	vector < vector < int > > bonds; //bonds are stored in a data structure of size nBonds x 2, where the 
	//two entries are the indices of the two involved beads from the monomer vector. for each mol, these include 0-nBead 
	vector < vector < int > > angles; 
	vector < vector < int > > dihedrals; 
	//add beads/topological features
	void addBead(Bead);
	void addBead(string, int, int, double, double, double); 
	void addBond(int, int); 
	void addAngle(int, int, int); 
	void addDihedral(int, int, int, int); 
	//get info about molecule
	int size(); 
	int nBond(); 
	int nAngle(); 
	int nDihedral(); 
	//do MC translation moves
	void beadTranslate(double, mt19937&); 
	void comTranslate(double delta, mt19937&);
	void pivot(double delta, mt19937&);
	void crankshaft(double delta, mt19937&); 
	void movedTrue(); //turns the m_moved flag to true, if bead moved this time
	bool moved(); //returns the m_moved flag in whatever state
        void resetMoved(); //returns m_moved flag to false;
	void printMol(ostream&); //print xyz coords of molecule's beads
	//funcs for virtual volume moves
	void scaleCoord(double[], double, int);
        void unscale(int);
}; 
#endif

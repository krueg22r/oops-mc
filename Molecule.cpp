#include "Molecule.h"
#include "Small.h"
#include "Bead.h"
#include <iostream> 
#include <iomanip> 
#include <stdlib.h> 
#include <Eigen/Dense>
using namespace std; 
using namespace Eigen; 

//standard constructor; 
Molecule::Molecule(){
	m_size = 0;
	m_nBond = 0; 
	m_nAngle = 0; 
	m_nDihedral = 0; 
	m_moved = false; 
}

//copy constructor
Molecule::Molecule(const Molecule& mol){
	m_size = 0; //these values will be incremented 
	m_nBond = 0; //as the add functions are called. 
	m_nAngle = 0;
        m_nDihedral = 0;
	//cout<<"invoking mol copy constructor!"<<endl; 
	m_moved = mol.m_moved;
	for(int i = 0; i < (int)mol.bonds.size(); i++){
                addBond(mol.bonds[i][0], mol.bonds[i][1]); 
        }
        for(int i = 0; i < (int)mol.angles.size(); i++){
                addAngle(mol.angles[i][0], mol.angles[i][1], mol.angles[i][2]); 
        }
        for(int i = 0; i < (int)mol.dihedrals.size(); i++){
                addDihedral(mol.dihedrals[i][0], mol.dihedrals[i][1], mol.dihedrals[i][2], mol.dihedrals[i][3]); 
        }
	for(int i = 0; i < (int)mol.monomers.size(); i++){
                addBead(mol.monomers[i]); 
        }
}

//adds a bead to the first empty slot in the bead array
void Molecule::addBead(Bead b){
	monomers.push_back(b); 
	m_size++; 
}

void Molecule::addBead(string symbol, int id, int charge, double x, double y, double z){
	monomers.push_back(Bead(symbol, id, charge, x, y, z)); 
	monomers[m_size].setCoord(0,0,x); 
	monomers[m_size].setCoord(0,1,y); 
	monomers[m_size].setCoord(0,2,z); 
	m_size++; 
}

void Molecule::addBond(int ind1, int ind2){
	vector < int > bond = {ind1, ind2}; 
	bonds.push_back(bond); 
	m_nBond++; 
}

int Molecule::nBond(){
	return m_nBond; 
}

void Molecule::addAngle(int ind1, int ind2, int ind3){
	vector < int > angle = {ind1, ind2, ind3}; 
	angles.push_back(angle); 
	m_nAngle++; 
}

void Molecule::addDihedral(int ind1, int ind2, int ind3, int ind4){
	vector < int > di = {ind1, ind2, ind3, ind4};
	dihedrals.push_back(di); 
	m_nAngle++;  
}

//returns number of beads in molecule
int Molecule::size(){
	return m_size; 
}

int Molecule::nAngle(){
	return m_nAngle; 
}

int Molecule::nDihedral(){
	return m_nDihedral; 
}

//move functions! 

//translate a random monomer (update new coords array of one bead)
void Molecule::beadTranslate(double delta, mt19937& ranGen){
	double vec[3]; 
	randSphere(vec, ranGen); 
	//cout<<" randomly chosen displacement is "<<vec[0]<<" "<< vec[1] << " " <<vec[2]<<endl; 
	int toMove = floor(m_size * (double)ranGen()/ranGen.max()); //pick monomer to move
	monomers[toMove].movedTrue(); //toggle moved flag on monomer to true
	for(int i = 0; i < 3; i++){
		double old = monomers[toMove].getCoord(0,i); 
		monomers[toMove].setCoord(1, i, old + delta * vec[i]); 
	}
}

//translate a whole fucking molecule
void Molecule::comTranslate(double delta, mt19937& ranGen){
	double vec[3];
        randSphere(vec, ranGen);
	for(int i = 0; i < m_size; i++){
		monomers[i].movedTrue(); 
		for(int k = 0; k < 3; k++){
			double old = monomers[i].getCoord(0,k);
			monomers[i].setCoord(1, k, old + delta * vec[k]); 
		}
	}
}

//pivot part of a molecule
void Molecule::pivot(double delta, mt19937& ranGen){
	int pivotPoint = floor((m_size - 2) * (double)ranGen()/ranGen.max() + 1); //cannot pivot from first bead (need bond vector) and useless to pivot from last 
	double angle = delta * ((double)ranGen()/ranGen.max() * 2 * M_PI - M_PI); 
	Vector3d bondVec; 
	for(int i = 0; i < 3; i++) bondVec(i) = monomers[pivotPoint].getCoord(0,i) - monomers[pivotPoint - 1].getCoord(0,i); 
	bondVec.normalize(); 
	AngleAxisd a = AngleAxisd(angle, bondVec);
	Matrix3d rot; 
	rot = a.toRotationMatrix(); //made rotation matrix 
	for(int i = pivotPoint + 1; i < m_size; i++){
		monomers[i].movedTrue(); 
		Vector3d v; //put coords into eigen vector object
		//and translate them so that the pivot point is the origin
		for(int k = 0; k < 3; k++){
			 v(k) = monomers[i].getCoord(0,k) - monomers[pivotPoint].getCoord(0,k); 
		}
		v = rot * v; //do rotation & put them back 
		for(int k = 0; k < 3; k++){
			//put in new coords, adding back the position of the pivot point. 
			 monomers[i].setCoord(1,k,v(k) + monomers[pivotPoint].getCoord(0,k)); 
		}
	}
}

//crankshaft. basically pivot for an internal section of the molecule
void Molecule::crankshaft(double delta, mt19937& ranGen){
	int first = floor((m_size-2) * (double)ranGen()/ranGen.max()); 
	int last = floor((m_size - first) * (double)ranGen()/ranGen.max()) + first; 
	Vector3d axis; 
	for(int i = 0; i < 3; i++) axis(i) = monomers[last].getCoord(0,i) - monomers[first].getCoord(0,i); 
	axis.normalize();
	double angle = delta * ((double)ranGen()/ranGen.max() * 2 * M_PI - M_PI);
	AngleAxisd a = AngleAxisd(angle, axis); 
	Matrix3d rot;
	rot = a.toRotationMatrix(); 
	for(int i = first + 1; i < last; i++){
		monomers[i].movedTrue(); //toggle moved flags for affected beads to true
		Vector3d v; //put coords into eigen vector object, make the first bead the origin
                for(int k = 0; k < 3; k++) v(k) = monomers[i].getCoord(0,k) - monomers[first].getCoord(0,k);
                v = rot * v; //do rotation & put them back, restoring normal origin 
                for(int k = 0; k < 3; k++) monomers[i].setCoord(1, k, v(k) + monomers[first].getCoord(0,k));
	}
	
}

void Molecule::movedTrue(){
	m_moved = true; 
}

void Molecule::resetMoved(){
	m_moved = false; 
}

bool Molecule::moved(){
	return m_moved; 
}

//prints accepted coordinates for a mol to either cout or a chosen ofstream
//xyz format 
void Molecule::printMol(ostream& out){
	for(int i = 0; i < m_size; i++){
		monomers[i].printXYZ(out); 
	}
}

//for virtual volume changes
//scales selected coord (0=x, 1=y, 2=z) of molecule's trial position
//trial coords are put back in original box before this happens
void Molecule::scaleCoord(double length[], double deltaL, int coord){
        for(int i = 0; i < m_size; i++){
		double oldCoord = monomers[i].getCoord(0,coord); 
		oldCoord -= length[coord]*floor(oldCoord/length[coord]); //put it back in the box
		oldCoord += oldCoord * deltaL; 
		monomers[i].setCoord(1, coord, oldCoord); 
        }
}

//reset appropriate coord for all molecules after scaling. 
void Molecule::unscale(int coord){
        for(int i = 0; i < m_size; i++){
                monomers[i].newToOld();
        }
}

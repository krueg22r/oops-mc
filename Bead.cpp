#include "Bead.h"
#include <cmath> 
#include <iomanip> 
#include <string> 
#include <iostream> 
#include <fstream> 
using namespace std; 

//default constructor exists because it's necessary, not because 
//you should use it
Bead::Bead(){
	for(int i = 0; i < 3; i++){
		posTrial[i] = 0; 
		posOld[i] = 0; 
	}
	m_charge  = 0; 
	m_symbol = "!!!"; 
	m_moved = false; 
}

//constructor meant for use during initial coordinate read
Bead::Bead(string symbol, int id, int charge, double x, double y, double z){
	posTrial[0] = x; 
	posTrial[1] = y; 
	posTrial[2] = z; 
	for(int i = 0; i < 3; i++) posOld[i] = posTrial[i]; 
	m_symbol = symbol; 
	m_charge = charge;
	m_moved = false;  
	m_id = id; 
}

//copy constructor
Bead::Bead(const Bead& b){
	for(int i = 0; i < 3; i++){
		posTrial[i] = b.posTrial[i]; 
	}
	for(int i = 0; i < 3; i++){
		posOld[i] = b.posOld[i]; 
	}	
	m_charge = b.m_charge; 
	m_type = b.m_type; 
	m_symbol = b.m_symbol;
	m_moved = b.m_moved;  
	m_id = b.m_id; 
}

//take the coordinates of a bead off another bead
//these coords are the only fields changed! (unlike copy constructor) 
void Bead::setAllCoord(const Bead& b){
	for(int i = 0; i < 3; i++){
                posTrial[i] = b.posTrial[i];
        }
        for(int i = 0; i < 3; i++){
                posOld[i] = b.posOld[i];
        }
}

//get bead attributes  
string Bead::symbol(){
	return m_symbol; 
}

int Bead::type(){
	return m_type; 
}

int Bead::charge(){
	return m_charge; 
}

int Bead::id(){
	return m_id; 
}

//get a coordinate from the old or new array
double Bead::getCoord(int arr, int index){
	if(arr == 0){
		 return posOld[index]; 
	}else if(arr == 1){
		 return posTrial[index]; 
	}else cout << "invalid array!"<<endl; 
	return 1; 
}

//set a coordinate in the old or new array
void Bead::setCoord(int arr, int index, double value){
	if(arr == 0) posOld[index] = value; 
	else if(arr == 1) posTrial[index] = value; 
	else cout << "invalid array!"<<endl; 
}

//set all coords to the coords given in array 
void Bead::setAllCoord(double pos[]){
	for(int i = 0; i < 3; i++){
		posOld[i] = pos[i]; 
		posTrial[i] = pos[i]; 
	}
}

//functions related to whether the bead was moved or not 

void Bead::movedTrue(){
	m_moved = true; 
}

void Bead::resetMoved(){
	m_moved = false; 
}

bool Bead::moved(){
	return m_moved; 
}

// particle information 

void Bead::setType(int type){
	m_type = type; 
}

void Bead::setID(int id){
	m_id = id; 
}

void Bead::setSymbol(string symbol){
	m_symbol = symbol; 
}

void Bead::setCharge(int charge){
	m_charge = charge; 
}

// adjust coordinates after move is accepted / rejected 

void Bead::newToOld(){
	for(int i = 0; i < 3; i++) posTrial[i] = posOld[i]; 
}

void Bead::oldToNew(){
	for(int i = 0; i < 3; i++) posOld[i] = posTrial[i]; 
}

void Bead::printXYZ(ostream& out){
	out << m_symbol; 
	out << setw(15)<<setprecision(7)<< posOld[0] << setw(15)<<setprecision(7)<< posOld[1] << setw(15)<<setprecision(7)<<posOld[2]<<endl; 
}

void Bead::printFull(ostream& out){
	out <<setw(6)<< m_symbol;
        out << setw(15)<<setprecision(7)<< posOld[0] << setw(15)<<setprecision(7)<< posOld[1] << setw(15)<<setprecision(7)<<posOld[2]<<setw(10)<<setprecision(3)<<m_charge<<endl; 
}

//calculates distance to the arg bead using TRIAL coords
double Bead::beadDist(Bead& b, double length[], int dimension){
	double dist = 0;
        for(int i = 0; i < 3; i++){
		double di = posTrial[i] - b.posTrial[i]; 
		if(i < dimension) di -= length[i] * round(di / length[i]); //periodic boundary conditions 
		dist += (di * di); 
	}
	return sqrt(dist); 
}

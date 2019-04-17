//Basic object describing a single particle
//contains position along with other attributes
#ifndef BEAD_H
#define BEAD_H
#include <iomanip> 
#include <string>
#include <vector> 
#include <fstream> 
#include <iostream> 
using namespace std; 

class Bead{
	string m_symbol; 
	int m_charge; 
	double posTrial[3]; //trial coordinates
	double posOld[3]; //accepted coordinates 
	int m_type; 
	int m_id; //unique identifier (not reused) 
	bool m_moved; //did the bead move this time? 
	public: 
	Bead(); 
	//constructor args are symbol, id, charge, x, y, z
	Bead(string, int, int, double, double, double); 
	Bead(const Bead&); //copy constructor 
	bool moved(); 
	void movedTrue(); //turns the m_moved flag to true, if bead moved this time
	void resetMoved(); //returns m_moved flag to false; 
	string symbol(); 
	int charge();
	int id(); 
	void setID(int); 
	void setSymbol(string); 
	void setCharge(int); 
	void setType(int); 
	//args are 0 or 1 (old or trial position) and 0, 1, or 2 (x, y, or z) 
	double getCoord(int, int); //args are 0 or 1 (old or trial position)
	void setCoord(int, int, double); //first 2 args same as getCoord, last one is desired new value
	void setAllCoord(double[3]);  //set xyz using a coord array
	void setAllCoord(const Bead&); //or another bead
	void oldToNew(); //set old coords to new ones or visa versa
	void newToOld(); 
	int type(); 
	void printXYZ(ostream&); //print just this bead's coords
	void printFull(ostream&); //print coords and charge 
	double beadDist(Bead&, double[], int); //bead calculates its distance to another bead
}; 
#endif

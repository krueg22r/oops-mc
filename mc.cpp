//main program for OOPS 
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <fstream>
#include <map>
#include <chrono> 
#include <random>  
#include "SimBox.h"

using namespace std; 



int main(){
	cout<<"starting program!"<<endl; 
	int every, equilLen, steps, dimension; 
	double delta, deltaLen, beta; 
	string label, runName, inFile, topFileName;
	//read in simple run parameters 
	cin >> label >> runName >> label >> inFile >> label >> topFileName; 
	cin >> label >> steps >> label >> every >> label >> equilLen; 
	cin >> label >> deltaLen >> label >> delta >> label >> beta >> label >> dimension; 
	if(beta <= 0){
		cout<<beta<<" is not an acceptable beta value. ending program!"<<endl; 
		exit(1); 
	}
	if(dimension > 3 || dimension < 2){
		cout<<dimension<<" is an invalid number of periodic dimensions! ending program!"<<endl; 
		exit(1); 
	}
	//initialize main prog's rand generator using time
	mt19937 ranGen; 
	unsigned s = chrono::system_clock::now().time_since_epoch().count();
        ranGen.seed(s);
	//open necessary input & output files
	ofstream energies, traj, frame, top; 
	ifstream coordFile, topFile; 
	string enFile = runName + "Energies.dat"; 
	string trajFile = runName + "Traj.xyz"; 
	string frameFile = runName + "LastFrame.dat"; 
	string topName = runName + "FinalTop.dat"; 
	energies.open(enFile.c_str()); 
	traj.open(trajFile.c_str());
	coordFile.open(inFile.c_str()); 
	topFile.open(topFileName.c_str());
	frame.open(frameFile.c_str());  
	top.open(topName.c_str()); 
	cout<<"all files opened!"<<endl; 
	//if they can't be opened, frown & quit program. 
	if(!coordFile){
		cout << "Initial coordinate file could not be opened! Ending program. :-("<<endl;
                exit(1);
        }
	if(!topFile){
		cout<< "Topology file could not be opened! Ending program. :-("<<endl;
                exit(1);
	}
	//declare SimBox object. 
	SimBox box(delta, beta, equilLen); 
	box.readTop(topFile); //reads bead connectivity info & initial bead/molecule numbers
	cout<<"read topology!"<<endl; 
	box.readCoords(coordFile); //reads coords off separate file
	cout<<"read coordinates!"<<endl; 
	box.printCoords(traj, 0); 
	//declare ForceField object
	ForceField field(beta, deltaLen, dimension, box.nMol()); 
	box.setLength(field); 
	//initialize all energies
	box.initEnergy(field);
	cout<<"field has initialized energies!"<<endl; 
	field.printHeader(energies); 
	box.printEnergy(field, energies, 0); //print initial "step 0" energies 
	//exit(1); 
	//start main MC loop
	for(int s = 1; s <= steps; s++){
		int r = ranGen(); 
		if(field.gc() && (r % field.gcEvery() == 0)){ //randomly choose whether to do a gc move, if doing that
//			cout<<"about to make gc move!"<<endl; 
			box.gcMove(field); 
		}else{
			box.makeMove(field); //try a translational move
		}
		if(s % every == 0){
			box.printCoords(traj, s);
			box.printEnergy(field, energies, s); 
		}
	}
	box.printLastFrame(frame); 
	box.printFinalTop(top); 
	box.runStats(); //prints translational move acceptance stats
	field.enStats(energies); //prints gc move acceptance stats (if applicable) 
	//cout<<"exiting program!"<<endl; 
	energies.close(); 
	return 0; 
}

#include "PairPot.h"
#include <cmath>
#include <vector> 

using namespace std; 

//default constructor
PairPot::PairPot(){
        m_deTrial = 0;
        m_totOld = 0;
}

//constructor w / max num of beads 
PairPot::PairPot(int nMol){
        m_deTrial = 0;
        m_totOld = 0;
}

//sets a pairwise energy value in either array
void PairPot::setVal(int arr, int key1, int key2, double val){
	if(arr == 0){
		eOldMap[make_pair(min(key1, key2), max(key1, key2))] = val; 
	}else if(arr == 1){
		eTrialMap[make_pair(min(key1, key2), max(key1, key2))] = val;
	}else{
		cout<<"invalid map requested!"<<endl; 
	}
}

//puts the value into all pair slots in both arrays
void PairPot::setAllVal(int key1, int key2, double val){
	eOldMap[make_pair(min(key1, key2), max(key1, key2))] = val;
	eTrialMap[make_pair(min(key1, key2), max(key1, key2))] = val;
}

//gets a pairwise energy value from either array
double PairPot::getVal(int arr, int key1, int key2){
	if(arr == 0){
		return  eOldMap[make_pair(min(key1, key2), max(key1, key2))]; 
	}else if(arr == 1){
		return  eTrialMap[make_pair(min(key1, key2), max(key1, key2))];
	}else{
		cout<<"invalid array!"<<endl;
		return 1; 
	}
	return 0; 
}

//fill in all energy maps with values for initial configuration
void PairPot::initialize(vector < Molecule >& mols, double length[], int dimension){
	cout<<"in pairPot initializer! "<<mols.size()<<" molecules. lengths "<<length[0]<<" "<<length[1]<<" "<<length[2]<<endl; 
	//first put in intramolecular pair energies
	for(int i = 0; i < (int)mols.size(); i++){
		for(int j = 0; j < mols[i].size()-1; j++){
			for(int k = j+1; k < mols[i].size(); k++){
				double en = pairEn(mols[i].monomers[j], mols[i].monomers[k], length, dimension); 
				m_totOld += en; 
				setAllVal(mols[i].monomers[j].id(), mols[i].monomers[k].id(), en); 
			}
		}
	}
	//now intermolecular pair energies
	for(int i = 0; i < (int)mols.size()-1; i++){
		for(int j = i+1; j < (int)mols.size(); j++){
			for(int k = 0; k < mols[i].size(); k++){
				for(int l = 0; l < mols[j].size(); l++){
					double en = pairEn(mols[i].monomers[k], mols[j].monomers[l], length, dimension); 
					m_totOld += en; 
					setAllVal(mols[i].monomers[k].id(), mols[j].monomers[l].id(), en); 
				}
			}
		}
	}
}

//initialize a new molecule (requires that the beads have proper IDs) Assumes mol is inserted at end of mol array! 
void PairPot::initMol(vector < Molecule >& mols, double length[], int dimension){
	for(int i = 0; i < mols[mols.size() - 1].size()-1; i++){
		for(int j = i+1; j < mols[mols.size() - 1].size(); j++){
			double en = pairEn(mols[mols.size() - 1].monomers[i], mols[mols.size() - 1].monomers[j], length, dimension); 
			m_totOld += en; 
			setAllVal(mols[mols.size() - 1].monomers[i].id(), mols[mols.size() - 1].monomers[j].id(), en); 
		}
	}
	//cout<<"now intermolecular energies."<<endl; 
	for(int i = 0; i < mols[mols.size() - 1].size()-1; i++){
		for(int j = 0; j < (int)mols.size()-1; j++){
			for(int k = 0; k < mols[j].size(); k++){	
				double en = pairEn(mols[mols.size()-1].monomers[i], mols[j].monomers[k], length, dimension); 
				m_totOld += en;
				setAllVal(mols[mols.size()-1].monomers[i].id(), mols[j].monomers[k].id(), en); 
			}
		}
	}
}

double PairPot::deltaEnMove(vector < Molecule >& mols, int activeMol, double length[], int dimension){
	m_deTrial = 0;
	//all pairs in mol that moved 
	for(int i = 0; i < (int) mols[activeMol].size()-1; i++){
		for(int j = i+1; j < mols[activeMol].size(); j++){
			if(mols[activeMol].monomers[j].moved() || mols[activeMol].monomers[i].moved()){
				double eNew = pairEn(mols[activeMol].monomers[j], mols[activeMol].monomers[i], length, dimension); 
				setVal(1, mols[activeMol].monomers[j].id(), mols[activeMol].monomers[i].id(), eNew); 
				m_deTrial += (eNew - getVal(0, mols[activeMol].monomers[j].id(), mols[activeMol].monomers[i].id())); 
			}
		}
	}
	//all pairs with other mols
	for(int i = 0; i < (int)mols.size(); i++){
		if(i != activeMol){
			for(int j = 0; j < mols[i].size(); j++){
				for(int k = 0; k < mols[activeMol].size(); k++){
					if(mols[activeMol].monomers[k].moved()){
						double eNew = pairEn(mols[activeMol].monomers[k], mols[i].monomers[j], length, dimension); 
						setVal(1, mols[activeMol].monomers[k].id(), mols[i].monomers[j].id(), eNew); 
						m_deTrial += (eNew - getVal(0, mols[activeMol].monomers[k].id(), mols[i].monomers[j].id())); 
					}
				}
			}
		}
	}
	return m_deTrial; 
}

void PairPot::resetMaps(vector < Molecule > & mols, int activeMol, bool accepted){
	//adjust totEn variable
	if(accepted){
		 m_totOld += m_deTrial;
	}
	//intramolecular
	for(int i = 0; i < (int) mols[activeMol].size()-1; i++){
		for(int j = i+1; j < mols[activeMol].size(); j++){
			if(accepted){
				setVal(0, mols[activeMol].monomers[i].id(), mols[activeMol].monomers[j].id(), getVal(1, mols[activeMol].monomers[i].id(), mols[activeMol].monomers[j].id())); 
			}else{
				setVal(1, mols[activeMol].monomers[i].id(), mols[activeMol].monomers[j].id(), getVal(0, mols[activeMol].monomers[i].id(), mols[activeMol].monomers[j].id())); 
			}
		}
	}
	for(int i = 0; i < (int)mols.size(); i++){
                if(i != activeMol){
                        for(int j = 0; j < mols[i].size(); j++){
                                for(int k = 0; k < mols[activeMol].size(); k++){
                                        if(mols[activeMol].monomers[k].moved()){
						if(accepted){
							setVal(0, mols[activeMol].monomers[k].id(), mols[i].monomers[j].id(), getVal(1, mols[activeMol].monomers[k].id(), mols[i].monomers[j].id())); 
						}else{
							setVal(1, mols[activeMol].monomers[k].id(), mols[i].monomers[j].id(), getVal(0, mols[activeMol].monomers[k].id(), mols[i].monomers[j].id()));
						}
					}
				}
			}
		}
	}
}

void PairPot::delMol(vector < Molecule > & mols, int toDel){
	for(int i = 0; i < mols[toDel].size(); i++){ //we're not worried about double counting, so do intra/inter mol pairs at same time
		for(int j = 0; j < (int)mols.size(); j++){
			for(int k = 0; k < mols[j].size(); k++){
				//adjust total energy 
				m_totOld -= eOldMap[make_pair(min(mols[toDel].monomers[i].id(), mols[j].monomers[k].id()), max(mols[toDel].monomers[i].id(), mols[j].monomers[k].id()))]; 
				//taken en value out of array
				eOldMap.erase(make_pair(min(mols[toDel].monomers[i].id(), mols[j].monomers[k].id()), max(mols[toDel].monomers[i].id(), mols[j].monomers[k].id()))); 	
				eTrialMap.erase(make_pair(min(mols[toDel].monomers[i].id(), mols[j].monomers[k].id()), max(mols[toDel].monomers[i].id(), mols[j].monomers[k].id())));
			}
		}
	}
}

double PairPot::totEn(){
	//cout<<"in PairPot totEn, En = "<<m_totOld<<endl; 
	return m_totOld; 
}


//calculate en of all pair interactions in box, using the trial coordinates. 
//do not update any arrays
double PairPot::allTrialEn(vector < Molecule >& mols, double length[], int dimension){
        double totEn = 0;
        //first intramolecular energies
        for(int i = 0; i < (int)mols.size(); i++){
                for(int j = 0; j < mols[i].size()-1; j++){
                        for(int k = j+1; k < mols[i].size(); k++){
                                totEn += pairEn(mols[i].monomers[j], mols[i].monomers[k], length, dimension);
                        }
                }
        }
        //now intermolecular pair energies
        for(int i = 0; i < (int)mols.size()-1; i++){
                for(int j = i+1; j < (int)mols.size(); j++){
                        for(int k = 0; k < mols[i].size(); k++){
                                for(int l = 0; l < mols[j].size(); l++){
                                        totEn += pairEn(mols[i].monomers[k], mols[j].monomers[l], length, dimension);
                                }
                        }
                }
        }
	//cout<<"totEn = "<<totEn<<endl; 
        return totEn;
}

void PairPot::printEn(ostream& out){
	out << setw(17) << setprecision(6) <<m_totOld; 
}

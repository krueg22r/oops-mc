#include "ExtPot.h"

using namespace std;

//default constructor
ExtPot::ExtPot(){
	m_deTrial = 0;
        m_totOld = 0;
}

void ExtPot::setVal(int arr, int key, double val){
	if(arr == 0){
		eOldMap[key] = val; 
	}else if(arr == 1){
		eTrialMap[key] = val; 
	}else{
		cout<<"invalid map requested!"<<endl;
        }
}

void ExtPot::setAllVal(int key, double val){
	eOldMap[key] = val; 
	eTrialMap[key] = val; 
}

double ExtPot::getVal(int arr, int key){
	if(arr == 0){
		return eOldMap[key]; 
	}else if(arr == 1){
		return eTrialMap[key]; 
	}else{
		cout<<"invalid map requested!"<<endl; 
	}
	return 1; 
}

void ExtPot::initialize(vector < Molecule >& mols, double length[], int dimension){
	for(int i = 0; i < (int)mols.size(); i++){
		for(int j = 0; j < mols[i].size(); j++){
			double en = beadEn(mols[i].monomers[j], length); 
			m_totOld += en; 
			setAllVal(mols[i].monomers[j].id(), en); 
		}
	}
}

//assume it's last mol being initialized
void ExtPot::initMol(vector < Molecule >& mols, double length[], int dimension){
	for(int i = 0; i < mols[mols.size()-1].size(); i++){
		double en = beadEn(mols[mols.size()-1].monomers[i], length); 
		m_totOld += en; 
		setAllVal(mols[mols.size()-1].monomers[i].id(), en); 
	}
}

//remove a molecule 
void ExtPot::delMol(vector < Molecule >& mols, int toDel){
	for(int i = 0; i < mols[toDel].size(); i++){
		m_totOld -= eOldMap[mols[toDel].monomers[i].id()]; 
		eOldMap.erase(mols[toDel].monomers[i].id()); 
		eTrialMap.erase(mols[toDel].monomers[i].id());
	}
}

double ExtPot::deltaEnMove(vector < Molecule >& mols, int activeMol, double length[], int dimension){
	m_deTrial = 0; 
	for(int i = 0; i < (int) mols[activeMol].size(); i++){
		if(mols[activeMol].monomers[i].moved()){
			double eNew = beadEn(mols[activeMol].monomers[i], length); 
			setVal(1, mols[activeMol].monomers[i].id(), eNew); 
			m_deTrial += eNew - getVal(0, mols[activeMol].monomers[i].id()); 
		}
	}
	return m_deTrial; 
}

void ExtPot::resetMaps(vector < Molecule >& mols, int activeMol, bool accepted){
	if(accepted) m_totOld += m_deTrial; 
	for(int i = 0; i < (int) mols[activeMol].size(); i++){
		if(mols[activeMol].monomers[i].moved()){
			if(accepted){
				setVal(0, mols[activeMol].monomers[i].id(), getVal(1, mols[activeMol].monomers[i].id())); 
			}else{
				setVal(1, mols[activeMol].monomers[i].id(), getVal(0, mols[activeMol].monomers[i].id()));
			}
		}
	}
}

double ExtPot::totEn(){
	return m_totOld; 
}

//total box energy due to external potential for trial positions
double ExtPot::allTrialEn(vector < Molecule >& mols, double length[], int dimension){
        double totEn = 0;
        for(int i = 0; i < (int)mols.size(); i++){
                for(int j = 0; j < mols[i].size(); j++){
                        totEn += beadEn(mols[i].monomers[j], length);
                }
        }
        return totEn;
}

void ExtPot::printEn(ostream& out){
        out << setw(17) << setprecision(8) <<m_totOld;
}

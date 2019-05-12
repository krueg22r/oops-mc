#include "BondPot.h" 

using namespace std; 

//construct with initial # of mols
BondPot::BondPot(int nMol){
	arr_eOld.resize(nMol); 
	arr_eTrial.resize(nMol); 
	m_deTrial = 0; 
	m_totOld = 0; 
}

//calc all initial molecular energies, fill in both vectors. 
void BondPot::initialize(vector <Molecule>& mols, double length[], int dimension){
	cout<<"initializing bondPot!"<<endl; 
	for(int i = 0; i < (int) mols.size(); i++){
		double en = enMol(mols[i], length, dimension); 
		arr_eOld[i] = en; 
		arr_eTrial[i] = en; 
		m_totOld += en; 
	}
	cout<<"done initializing bondPot!"<<endl; 
}

//initializes the last mol in the array (presumed to be new one) 
void BondPot::initMol(vector <Molecule>& mols, double length[], int dimension){
	double en = enMol(mols[mols.size()-1], length, dimension);
	arr_eOld.push_back(en); 
	arr_eTrial.push_back(en); 
	m_totOld += en;  
}

// if we removed a molecule, get rid of its energies 
void BondPot::delMol(int toDel){
	m_totOld -= arr_eOld[toDel]; 
	arr_eOld.erase(arr_eOld.begin() + toDel); 
	arr_eTrial.erase(arr_eTrial.begin() + toDel); 
}

void BondPot::setVal(int arr, int ind, double val){
        if(arr == 0){
                arr_eOld[ind] = val;
        }else if(arr == 1){
                arr_eTrial[ind] = val;
        }else{
                cout<<"invalid array!"<<endl;
        }
}

double BondPot::getVal(int arr, int ind){
        if(arr == 0){
              return  arr_eOld[ind];
        }else if(arr == 1){
              return arr_eTrial[ind];
        }else{
                cout<<"invalid array!"<<endl;
                return 1;
        }
}

// after a move has occured, adjust energy according to whether 
// accepted or not  
void BondPot::resetEn(int activeMol, bool accepted){
	if(accepted){
		arr_eOld[activeMol] = arr_eTrial[activeMol]; 
		m_totOld += m_deTrial; 
	}else{
		arr_eTrial[activeMol] = arr_eOld[activeMol]; 
	}
	m_deTrial = 0; //energy difference has now been processed. 
}

//returns energy for all bonds in box, based on trial position
double BondPot::allTrialEn(vector <Molecule>& mols, double length[], int dimension){
        double totEn = 0;
        for(int i = 0; i < (int) mols.size(); i++){
                totEn += enMol(mols[i], length, dimension);
        }
        return totEn;
}

double BondPot::totEn(){
	return m_totOld; 
}

// print this energy contribution 
void BondPot::printEn(ostream& out){
        out << setw(17) << setprecision(6) <<m_totOld;
}

#include "Spring.h" 

using namespace std; 

Spring::Spring(int nMol) : BondPot(nMol){
	cout<<"in spring constructor!"<<endl; 
	readParam(); 
}

void Spring::readParam(){
	double kBond, r0; 
	string label; 
	cin >> label >> kBond >> label >> r0; 
	m_kBond = kBond; 
	m_r0 = r0; 
}

double Spring::enMol(Molecule& mol, double length[], int dimension){
	double en = 0; 
	for(int i = 0; i < mol.nBond(); i++){
		int ind1 = mol.bonds[i][0]; //indices of 1st and 2nd bond beads 
		int ind2 = mol.bonds[i][1]; //in monomer vec
		double r = mol.monomers[ind1].beadDist(mol.monomers[ind2], length, dimension);//current bond dist
		en += 0.5 * m_kBond * (r - m_r0) * (r - m_r0); 
	}
	return en; 
}

//returns a bond length from appropriate distribution for harmonic bond
double Spring::ranBond(double beta, mt19937& ranGen){
	double len = 0; 
	double sigma = sqrt(1/(beta * m_kBond));
        double a = (m_r0 + 3*sigma) * (m_r0 + 3*sigma);
	bool ready = false;
        while(!ready){
		len = gasdev(m_r0, sigma, ranGen); 
		ready = ((double)ranGen()/ranGen.max() < (len * len/a)); 
	}
	return len; 
}

double Spring::deltaEnMove(vector < Molecule >& mols, double length[], int dimension, int molInd, int moveType){
	if(moveType==1 && mols[molInd].size() > 1){
		double newEn = enMol(mols[molInd], length, dimension); 
		setVal(1, molInd, newEn);
                m_deTrial = (newEn - getVal(0, molInd)); 
		return m_deTrial; 
	}else{
		return 0; 
	}
	return 0; 
}

void Spring::printHeader(ostream& out){
	out << setw(17) << "SpringBonds"; 
}

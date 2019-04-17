#include "CoulLjCut.h"

CoulLjCut::CoulLjCut(int nMol) : PairPot(nMol){
	readParam(); 
}

CoulLjCut::CoulLjCut() : PairPot(){
	readParam(); 
}

void CoulLjCut::readParam(){
	string label; 
	double ljCut, coulCut, dielectric; 
	//read simple pot parameters
	cin >> label >> ljCut >> label >> coulCut >> label >> dielectric; 
	m_ljCut = ljCut; 	
	m_dielectric = dielectric; 
	m_coulCut = coulCut; 
	//read sigma + epsilon associated with each symbol 
	string symbol; 
	double sigma, epsilon; 
	cin >> symbol; 
	while(symbol != "end"){
		cin >> label >> sigma >> label >> epsilon; 
		sigmas[symbol] = sigma; 
		epsilons[symbol] = epsilon; 
		cin >> symbol; 
	}
}

//energy between two beads
double CoulLjCut::pairEn(Bead& b1, Bead& b2, double length[], int dimension){
	double en = 0; 
	double r = b1.beadDist(b2, length, dimension); 
	if(r < m_ljCut){
		if(r < 0.00000001){
			cout<<"dist between bead "<<b1.id()<<" and bead "<<b2.id()<<" is "<<r<<endl; 
			cout<<"bead coords: "<<endl; 
			b1.printXYZ(cout); 
			b2.printXYZ(cout); 
			cout<<"pair energy requested between identical beads! ending program."<<endl; 
			exit(1); 
		}
		double sig = (sigmas[b1.symbol()] + sigmas[b2.symbol()]) / 2; //arithmetic mixing
		double ep = sqrt((epsilons[b1.symbol()] * epsilons[b2.symbol()])); //geometric mixing 
		double r6 = pow((sig/r),6); 
		en += 4 * ep * (r6*r6 - r6); 
	}
	if(r < m_coulCut){
		en += b1.charge() * b2.charge() / (r * m_dielectric) * 14.36; //const converts value to eV 
	}
	return en; 
}

void CoulLjCut::printHeader(ostream& out){
	out<<setw(19)<<"CoulLjCut";  
} 

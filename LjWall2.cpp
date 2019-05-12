#include "LjWall2.h" 
// This uses the potential obtained by integrating over 
// a 2-d surface made up of LJ particles. 
// Rather than 6-12, it ends up being 3-9. 

using namespace std; 

LjWall2::LjWall2(){
	readParam(); 
}

void LjWall2::readParam(){
	string label; 
	double sigWall, epWall; 
	double cut; 
	string symbol; 
	double sigma, epsilon; 
	cin >> label >> cut; 
	cin >> label >> sigWall >> label >> epWall; 
	cin >> symbol;
        while(symbol != "end"){
                cin >> label >> sigma >> label >> epsilon;
                //cout<<sigma<<" "<<epsilon<<endl; 
                sigmas[symbol] = sigma;
                epsilons[symbol] = epsilon;
                cin >> symbol;
        }
	m_cut = cut; 
	m_sigWall = sigWall; 
	m_epWall = epWall; 
	
}

//lj walls are perpendicular to the z-direction! 

double LjWall2::beadEn(Bead& b, double length[]){
	double zCoord = b.getCoord(1, 2); //z coord from trial coord array
	double en = 0;
	if(zCoord < 0 || zCoord > length[2]){ //if this z coord falls out of the box
		return 900000; //return a really big energy
	} 
	if(zCoord < m_cut){
		double sig = sigmas[b.symbol()]; 
		double ep = epsilons[b.symbol()]; 
		double r3 = pow((sig/zCoord),3);
                en +=  4 * ep * (r3*r3 - r3);
	}
	double rRight = length[2] - zCoord; //dist from right-hand wall
        if(rRight < m_cut){
		double sig = sigmas[b.symbol()]; 
                double ep = epsilons[b.symbol()];
		double r3 = pow((sig/rRight),3);
		en +=  4 * ep * (r3*r3 - r3);		
	}
	return en; 
}

//calculate force PER UNIT AREA due to molecule-wall interaction 
double LjWall2::calcForce(vector < Molecule >& mols, double length[]){
	double force = 0; 
	for(int i = 0; i < (int)mols.size(); i++){
		for(int j = 0; j < mols[i].size(); j++){
			double zCoord = mols[i].monomers[j].getCoord(1,2); 
			if(zCoord < 0 || zCoord > length[2]){
				cout<<"bead outside the box for force calculation! endng program!"<<endl; 
				exit(1); 
			}
			if(zCoord < m_cut){
				double sig = sigmas[mols[i].monomers[j].symbol()];
				double ep = epsilons[mols[i].monomers[j].symbol()];
				force +=  4 * ep * (3*pow(sig,3)*(1/pow(zCoord,4)) - 9 * pow(sig,9)*(1/pow(zCoord,10)));
			}
			double rRight = length[2] - zCoord; //dist from right-hand wall
                	if(rRight < m_cut){
				double sig = sigmas[mols[i].monomers[j].symbol()];
                                double ep = epsilons[mols[i].monomers[j].symbol()];
                                force += 4 * ep * (3*pow(sig,3)*(1/pow(rRight,4)) - 9 * pow(sig,9)*(1/pow(rRight,10)));
			}
		}
	}
	force /= (length[0]*length[1]*2); //force per unit area 
	force *= (16.02); //convert from messy units to units of 10^10 Pa
	return (-1*force); 
}

void LjWall2::printHeader(ostream& out){
        out<<setw(17)<<"LjWall2"<<setw(17)<<"LiquidForce"<<setw(17)<<"AvgForce"; 
}



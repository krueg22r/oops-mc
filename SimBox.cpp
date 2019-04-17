#include "SimBox.h" 
using namespace std; 

//construct SimBox with some initial parameters
SimBox::SimBox(double delta, double beta, int equilLen){
	m_delta = delta; 
	m_beta = beta; 
	m_equilLen = equilLen; 
	m_idCounter = 0; 
	m_nChain = 0; 
	m_rgTot = 0; 
	m_rgCounter = 0; 
	for(int i = 0; i < 5; i++){
		accepted[i] = 0; 
		attempted[i] = 0; 
	}
	for(int i = 0; i < 3; i++){
		length[i] = 0; 
	}
	//initialize random number generator using time-based seed 
	unsigned s = chrono::system_clock::now().time_since_epoch().count(); 
	//unsigned s = 1987449896; 
	ranGen.seed(s); 
	//cout<<"initialized ranGen! seed = "<<s<<endl; 
	//if there's a bug that only happens sometimes, it's a good idea 
	//to use the same seed each time or print the seed to use it again.  
}

//to be called from main to initialize FF, because that requires mol array access
void SimBox::initEnergy(ForceField& field){
	field.initFF(mols); 
}

//calculates average value of sqrt(Rg^2) for ALL CHAINS (mols with nBead > 1)
//based on positions just accepted in last step
double SimBox::radGy(){
	double rgAvg = 0; 
	for(int i = 0; i < (int)mols.size(); i++){
		if(mols[i].size() > 1){ //if this is a chain molecule
			double rgMol = 0; 
			for(int j = 0; j < mols[i].size()-1; j++){
				for(int k = j+1; k < mols[i].size(); k++){
					double r = mols[i].monomers[j].beadDist(mols[i].monomers[k], length, 0); 
					//effect of setting dimension to 0: this doesn't apply
					//so we'll always calc distances in same molecule image
					rgMol += (r*r); 
				}
			}
			rgMol /= (mols[i].size() * mols[i].size()); 
			rgAvg += rgMol; 
		}
	}
	rgAvg = sqrt(rgAvg / m_nChain); 
	return rgAvg; 
}


void SimBox::printEnergy(ForceField& field, ofstream& energies, int step){
	field.printEnergies(mols, energies, m_equilLen, step, m_nTot); 
	if(step > m_equilLen){
		m_rgCounter++; 
		m_rgTot += radGy(); 	
	}
	energies<<setw(14)<<setprecision(6)<< m_rgTot / m_rgCounter<<endl; 
}

//returns an ID suitable for use with a new bead
int SimBox::newID(){
	pair<set<int>::iterator, bool> ret; 
	bool inserted = false; 
	while(!inserted){ //until a new ID has been successfully inserted
		ret = boxIDs.insert(m_idCounter); //try to insert idCounter
		inserted = ret.second; //check if it was inserted
		m_idCounter++; //increment idCounter to give the next id to try
	}
	//at end of loop, idCounter should be greater than the greatest number in set
	return (m_idCounter-1); 
	//return value is the number just inserted into set, is now ready to give to a bead
}

//takes an ifstream hooked up to a coordinate file (slight variation on xyz format) 
void SimBox::readCoords(ifstream& coordFile){
	int charge, molNum; 
	string dummy, symbol; 
	getline(coordFile, dummy); 
	getline(coordFile, dummy); 
	double x, y, z; 
	while(coordFile >> molNum >> symbol >> x >> y >> z >> charge){
		cout<<molNum<<" "<<symbol<<" "<<x<<" "<<y<<" "<<z<<" "<<charge<<endl; 
		if(molNum < 0 || molNum > m_nMol){ //if given molID is greater than nMol or negative, quit
			cout<<"Molecule ID "<<molNum<<" is out of bounds!"<<endl; 
			cout<<"Ending program. :-("<<endl; 
			exit(1); 
		}else{ //given mol ID is valid, add bead to that mol
			int id = newID(); 
			mols[molNum].addBead(Bead(symbol, id, charge, x, y, z)); 
		}
	}
	//count number of chain molecules
	for(int i = 0; i < (int)mols.size(); i++){
		if(mols[i].size() > 1) m_nChain++; 
	}
}

//reads in initial bond, angle, dihedral information, preceded by some species info and move probabilities
void SimBox::readTop(ifstream& topFile){
	string label; 
	int nMol, nTot; 
	double lx, ly, lz; 
	topFile >> label >> nTot >> label >> nMol; 	
	cout<<"nMol = "<<nMol<<endl; 
	//read in box dimensions
	topFile >> label >> lx >> label >> ly >> label >> lz; 
	length[0] = lx; 
	length[1] = ly; 
	length[2] = lz; 
	double pBeadTrans, pComTrans, pPivot, pCrankshaft; 
	topFile >> label >> pBeadTrans >> label >> pComTrans >> label >> pPivot >> label >> pCrankshaft; 
	moveProbs[0] = pBeadTrans; 
	moveProbs[1] = pComTrans; 
	moveProbs[2] = pPivot; 
	moveProbs[3] = pCrankshaft; 
	m_nTot = nTot; 
	m_nMol = nMol; 
	int ind1, ind2, ind3, ind4, molID; 
	topFile >> label >> molID; 
	//populate mol array with initial number of empty mols
	for(int i = 0; i < nMol; i++){
		mols.push_back(Molecule()); 
	}
	while(molID != -1){
		topFile >> ind1 >> ind2; 
		//cout<<molID<<" "<<ind1<<" "<<ind2<<endl; 
		mols[molID].addBond(ind1, ind2); 
		topFile >> molID; 
	}
	topFile >> label >> molID; 
	while(molID != -1){
                topFile >> ind1 >> ind2 >> ind3; 
		mols[molID].addAngle(ind1, ind2, ind3); 
		topFile >> molID; 
	}
	topFile >> label >> molID; 
	while(molID != -1){
                topFile >> ind1 >> ind2 >> ind3 >> ind4; 
		mols[molID].addDihedral(ind1, ind2, ind3, ind4); 
		topFile >> molID; 
	}
}

void SimBox::setLength(ForceField& field){
	for(int i = 0; i < 3; i++){
		field.length[i] = length[i]; 
	}
}

//Attempts to perform one of the available translational moves on a random molecule 
//Accepts/rejects, adjusts all coords & energy arrays appropriately, updates attempted/accepted counters
void SimBox::makeMove(ForceField& field){
	int species = (double)ranGen() / ranGen.max() * m_nMol; //randomly pick mol to move
	mols[species].movedTrue(); //toggle mol's moved flag to true
	int moveType = 0; 
	if(mols[species].size() == 1){ //single-particle molecule
		mols[species].beadTranslate(m_delta, ranGen); 
	}else{ //chain molecule
		double r = (double)ranGen() / ranGen.max(); 
		double cumProb = moveProbs[0];
                while(cumProb < r){
                        moveType++;
                        cumProb += moveProbs[moveType];
                }
                moveType += 1; //movetypes for polymers run from 1 (monomer translate) to 4 (crankshaft)
                switch (moveType){
                        case 1:
                                mols[species].beadTranslate(m_delta, ranGen);
                                break;
                        case 2:
                                mols[species].comTranslate(m_delta, ranGen);
                                break;
                        case 3:
                                mols[species].pivot(m_delta, ranGen);
                                break;
                        case 4:
                                mols[species].crankshaft(m_delta, ranGen);
                                break;
                }
        }      
        attempted[moveType]++;
	double dE = field.processMove(mols, species, moveType); 
	bool accept = ((double)ranGen() / ranGen.max() < (exp(-1*m_beta*dE))); //accept or reject based on metropolis criterion 
	field.resetEnergies(mols, accept, species); //make trial energies old energies
	if(accept){
		for(int i = 0; i < mols[species].size(); i++){ //make trial positions old positions
			mols[species].monomers[i].oldToNew(); 	
		}
		accepted[moveType]++; 
	}else{
		for(int i = 0; i < mols[species].size(); i++){ //reset trial positions to old positions
			mols[species].monomers[i].newToOld(); 
		}	
	}
	for(int i = 0; i < mols[species].size(); i++){
		mols[species].monomers[i].resetMoved(); 
	}
}

void SimBox::gcMove(ForceField& field){
	bool ins = ranGen() % 2; //choose whether to attempt insertion or deletion
	if(ins){ //attempt insertion
		bool accept = field.trialIns(mols, m_beta, ranGen); //makes mol w/right coords. adds to end of mol vec
		if(accept){
			m_nTot += mols[mols.size()-1].size(); 
			//cout<<"in SimBox, insertion accepted! num mols now = "<<mols.size()<<endl;  
			for(int i = 0; i < mols[mols.size()-1].size(); i++){ //give new beads proper IDs 
				mols[mols.size()-1].monomers[i].setID(newID()); 
			}
			m_nMol++; //adjust internal variable
			m_nChain++; 
			field.addMol(mols); //do all energy initialization (also adjusts energy totals).  
		}
	}else{ //attempt deletion
		int delMol = field.trialDel(mols, m_nChain, m_beta, ranGen); //returns -1 if deletion unsuccessful, index of mol deleted if successful 
		//also does all adjustment to energy arrays + total energy trackers
		if(delMol != -1){
	//		}
			m_nTot -= mols[delMol].size(); 
			m_nChain--; 
			m_nMol--; //adjust internal variable
			mols.erase(mols.begin()+delMol); //delete chosen mol from mol vector  
		}
	}
}

//takes ofstream hooked up to bead trajectory file, print a frame of coords in exact xyz format
void SimBox::printCoords(ofstream& traj, int step){
	if(step > m_equilLen){ //we don't print coords during equilibration! 
		traj << m_nTot<<endl;
		traj <<"STEP: "<<step<<endl;
		for(int i = 0; i < m_nMol; i++){
			mols[i].printMol(traj); 
		}
	}
}

void SimBox::printLastFrame(ofstream& frame){
	frame << m_nTot<<endl;
	frame<<" "<<endl; 
	for(int i = 0; i < m_nMol; i++){
		for(int j = 0; j < mols[i].size(); j++){
			frame<<i; 
			mols[i].monomers[j].printFull(frame); 
		}
	}
}

//prints a topology file based on final config suitable for use in another run 
void SimBox::printFinalTop(ofstream& top){
	top <<"nTot: "<<m_nTot<<endl; 
	top <<"nMol: "<<m_nMol<<endl; 
	top << "xLen: "<<length[0]<<endl; 
	top << "yLen: "<<length[1]<<endl; 
	top << "zLen: "<<length[2]<<endl; 
	top <<"pBeadTranslate: "<<moveProbs[0] <<endl; 
	top <<"pComTranslate: "<<moveProbs[1] << endl; 
	top << "pPivot: "<<moveProbs[2]<<endl; 
	top << "pCrankshaft: "<<moveProbs[3]<<endl; 
	top<<"bonds:"<<endl; 
	for(int i = 0; i < m_nMol; i++){
		for(int j = 0; j < (int)mols[i].bonds.size(); j++){
			top<<i<<" "<<mols[i].bonds[j][0]<<" "<<mols[i].bonds[j][1]<<endl; 
		}
	}
	top<<"-1"<<endl; 
	top<<"angles:"<<endl; 
	for(int i = 0; i < m_nMol; i++){
                for(int j = 0; j < (int)mols[i].angles.size(); j++){
			top<<i<<" "<<mols[i].angles[j][0]<<" "<<mols[i].angles[j][1]<<" "<<mols[i].angles[j][2]<<endl; 
		}
	}	
	top<<"-1"<<endl; 
	top<<"dihedrals"<<endl; 
	for(int i = 0; i < m_nMol; i++){
                for(int j = 0; j < (int)mols[i].dihedrals.size(); j++){
			top<<i<<" "<<mols[i].dihedrals[j][0]<<" "<<mols[i].dihedrals[j][1]<<" "<<mols[i].dihedrals[j][2]<<" "<<mols[i].dihedrals[j][3]<<endl;
		}
	}
	top<<"-1"<<endl; 
}

//print stats about acceptance of translational moves to stdout. 
void SimBox::runStats(){
        cout<<"run complete! "<<accepted[0] / (double) attempted[0] * 100<<" \% of ion bead translate moves accepted,"<<endl;
        cout<<accepted[1] / (double) attempted[1] * 100<<" \% of polymer bead translate moves accepted"<<endl;
        cout<<accepted[2] / (double) attempted[2] * 100<<" \% of com translate moves accepted, ";
        cout<<accepted[3] / (double) attempted[3] * 100<<" \% of pivot moves accepted,"<<endl;
        cout<<"and "<<accepted[4] / (double) attempted[4] * 100<<" \% of crankshaft moves accepted"<<endl;
}
 
//self-explanatory get functions
int SimBox::nMol(){
	return m_nMol; 
}

int SimBox::nTot(){
	return m_nTot; 
}

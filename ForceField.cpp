#include "ForceField.h" 
using namespace std; 

ForceField::ForceField(double beta, double deltaLen, int dimension, int nMol){
        m_dimension = dimension;
	m_beta = beta; 
	m_deltaLen = deltaLen; 
	enCounter = 0; 
	totDense = 0; 
	totForce = 0; 
	extDef = false; 
	for(int i = 0; i < 3; i++){
		pPlus[i] = 0; 
		pMinus[i] = 0; 
	}
	readParam(nMol); 	
}

//reads in gc params if needed for run + potentials, adds potentials
void ForceField::readParam(int nMol){
	string label, chainSymbol; 
	bool gc; 
	int chainLen, chainCharge, chargeFrac, gcEvery, nVec; 
	double qT, chemPot; 
	cin >> label >> gc; 
	m_gc = gc; 
	if(m_gc){ //read in gc stuff if needed  
		cin >> label >> chainLen >> label >> chainCharge >> label >> chargeFrac >> label >> chainSymbol; 
		cin >> label >> gcEvery >> label >> nVec; 
		cin >> label >> qT >> label >> chemPot; 
		m_chainLen = chainLen; 
		m_chainCharge = chainCharge; 
		m_chargeFrac = chargeFrac; 
		m_chainSymbol = chainSymbol; 
		m_gcEvery = gcEvery; 
		m_nVec  = nVec; 
		m_qT = qT; 
		m_chemPot = chemPot; 
		delAccept = 0; 
		insAccept = 0; 
		delAttempt = 0; 
		insAttempt = 0; 
	}
	//set potentials
	string potName; 
	//first pairwise one
	cin >> label >> potName; 
	cout<<potName<<endl; 
	if(potName == "CoulLjCut"){
		pairPot = new CoulLjCut(); 
	}else if(potName == "LjCut"){
		pairPot = new LjCut(); 
	}else{
		cout<<"Undefined pair potential! Exiting."<<endl; 
		exit(1); 
	}
	//now the bonded one
	cin >> label >> potName;
	cout<<potName<<endl;  
	if(potName == "Spring"){
		cout<<"attempting to add spring pot!"<<endl; 
		bondPot = new Spring(nMol); 
	}else{
		cout<<"Undefined bonded potential! Exiting."<<endl; 
                exit(1); 
	}
	cin >> label >> potName;
	cout<<"label and potName: "<<endl; 
	cout<<label<<" "<<potName<<endl; 
	//now external potential
	if(potName == "LjWall2"){
		cout<<"attempting to add LjWall2 extPot!"<<endl; 
		extDef = true; 
		extPot = new LjWall2(); 
	}
}

//initialize all energy maps / vectors in potentials, set up gc if necessary
void ForceField::initFF(vector < Molecule >& mols){
	pairPot->initialize(mols, length, m_dimension); 
	cout<<"initialized pairPot!"<<endl; 
	bondPot->initialize(mols, length, m_dimension);
	cout<<"initialized bondPot!!"<<endl; 
	if(extDef){
		 extPot->initialize(mols, length, m_dimension);
		cout<<"initialized extPot"<<endl; 
	}
	//now add beads to gc trial structures. no proper IDs for beads not in SimBox! 
	if(m_gc){
		for(int i = 0; i < m_chainLen; i++){
			Bead b(m_chainSymbol, -1, m_chainCharge, 0, 0, 0); 
			trialChain.push_back(b); 
		}
		for(int i = 0; i < m_nVec; i++){
			Bead b(m_chainSymbol, -1, m_chainCharge, 0, 0, 0);
			trialVecs.push_back(b); 
		}
		wFactors = new double[m_nVec]; //will be initialized in func before use
	}
	cout<<"FF initialized! size of trialChain is "<<trialChain.size()<<" size of trialVecs is "<<trialVecs.size()<<endl; 
}

double ForceField::processMove(vector < Molecule >& mols, int activeMol, int moveType){
	double deTot = 0; 
	deTot += pairPot->deltaEnMove(mols, activeMol, length, m_dimension);
	deTot += bondPot->deltaEnMove(mols, length, m_dimension, activeMol, moveType); 
	if(extDef) deTot += extPot->deltaEnMove(mols, activeMol, length, m_dimension);
	return deTot; 
}

void ForceField::resetEnergies(vector < Molecule >& mols, bool accepted, int activeMol){
	bondPot->resetEn(activeMol, accepted); 
	pairPot->resetMaps(mols, activeMol, accepted); 
	if(extDef) extPot->resetMaps(mols, activeMol, accepted);
}

void ForceField::printHeader(ofstream& energies){
	energies <<setw(9)<< "#step";
	pairPot->printHeader(energies); 
	bondPot->printHeader(energies); 
	if(extDef){
		 extPot->printHeader(energies); 
	}
	energies<<setw(17)<<"TotalEnergy"; 	
	if(m_gc){
		energies<<setw(10)<<"nMol"; 	
		 energies<<setw(14)<<"MolDensity"; 
		energies<<setw(14)<<"AvgDensity"; 
	}
	energies<<setw(14)<<"Pxx"<<setw(14)<<"Pyy"<<setw(14)<<"Pzz"; 
	energies<<setw(14)<<"|Rg|"; 
	energies<<endl; 
}

double ForceField::getDeltaV(int coord){
	double area = 1;
	for(int q = 0; q < 3; q++){
		if(coord != q) area *= length[q]; //calculate area of perpendicular face
	}
	double deltaV = area * length[coord]*m_deltaLen;  
	return deltaV;
}

//calculate all diagonal pressure compontents using method from de Miguel 2006
void ForceField::calcPress(vector < Molecule >& mols, ofstream& energies, int nTot, double currEn){
        double origVol = length[0]*length[1]*length[2];
        for(int i = 0; i < 3; i++){ //loop over coords 
                for(int j = 0; j < (int)mols.size(); j++){
                        mols[j].scaleCoord(length, m_deltaLen, i);
                }
		double lengthScaled[3]; 
		for(int j = 0; j < 3; j++){ //initialize array w/scaled length values
			lengthScaled[j] = length[j];  
			if(j == i){ //if this is the coord we're dealing with this loop, scale it. 
				lengthScaled[j] += length[j]*m_deltaLen;
			}
		}
                double deltaEnPlus = 0;
                double deltaPair = pairPot->allTrialEn(mols, lengthScaled, m_dimension);
                double deltaBond = bondPot->allTrialEn(mols, lengthScaled, m_dimension);
                if(extDef){
			 double deltaExt = extPot->allTrialEn(mols, lengthScaled, m_dimension);
			deltaEnPlus += deltaExt; 
		}
		deltaEnPlus += (deltaPair + deltaBond); 
		//subtract off current energy
		deltaEnPlus -= currEn; 
                double deltaV = getDeltaV(i);
		double innerTermPlus = (pow((1 + deltaV/origVol), nTot)*exp(-1*m_beta*deltaEnPlus)); 
                pPlus[i] += innerTermPlus;
                for(int j = 0; j < (int)mols.size(); j++){
                        mols[j].scaleCoord(length, -1*m_deltaLen, i);
                }
		for(int j = 0; j < 3; j++){ 
                        if(j == i){ //if this is the coord we're dealing with this loop, scale it. 
                                lengthScaled[j] = length[j] - length[j]*m_deltaLen;
                        }
                }
                double deltaEnMinus = 0;
		deltaPair = pairPot->allTrialEn(mols, lengthScaled, m_dimension);
		deltaBond = bondPot->allTrialEn(mols, lengthScaled, m_dimension);
		if(extDef){
                         double deltaExt = extPot->allTrialEn(mols, lengthScaled, m_dimension);
                        deltaEnMinus += deltaExt;
                }
		deltaEnMinus += (deltaPair + deltaBond); 
		//subtract off current energy
		deltaEnMinus -= currEn; 
		double innerTermMinus = (pow((1 - deltaV/origVol), nTot)*exp(-1*m_beta*deltaEnMinus));
                pMinus[i] += innerTermMinus;
		//now undo scaling
		for(int j = 0; j < (int)mols.size(); j++){
			mols[j].unscale(i); 
		}
		//instantaneous pressure tensor value for this coord in units of 10^10 Pa
		double pII = (1/m_beta)*16.02*0.5*(1/deltaV * log(innerTermPlus) + -1/deltaV * log(innerTermMinus)); 
		energies<<setw(14)<<setprecision(6)<<pII; 
        }
}

//print energy contributions from elements of FF, calculate pressure tensor components
// if GC run print stats about density
void ForceField::printEnergies(vector < Molecule >& mols, ofstream& energies, int equilLen, int step, int nTot){
	if(step > equilLen){
        	enCounter++;
	}
        energies << setw(9) << step;
        double tot = 0;
        pairPot->printEn(energies); //print pairwise energies 
        tot += pairPot->totEn();
        bondPot->printEn(energies); //print bonded energies
        tot += bondPot->totEn();
        if(extDef){ //if we've defined an external potential that has associated force
                tot += extPot->totEn();
                extPot->printEn(energies);
                //calculate force due to this potential
                double force = extPot->calcForce(mols, length);
                totForce += force;
                double forceToNow = totForce / enCounter;
                energies << setw(17) << setprecision(7)<< force << setw(17) << setprecision(7)<< forceToNow;
                //calculate perpendicular component of pressure using virtual volume scaling. 
        }
        energies << setw(17) << tot;
        if(m_gc){
                energies << setw(10) << mols.size();
                double density = mols.size() / (length[0]*length[1]*length[2]);
		if(step > equilLen){
	                totDense += density;
		}
                double densToNow = totDense / enCounter;
                energies << setw(14) <<setprecision(7) << density << setw(14) << setprecision(7)<< densToNow;
        }
        //calculate pressure contribs and add into arrays but don't print pressure tensor components
	if(step > equilLen){ //only does this after equil period. numbers are liable to be a little crazy  
		//if configuration is not decently equilibrated. 
	        calcPress(mols, energies, nTot, tot);
	}
}

//calculate energy of all of a bead's "external" interactions w/other beads (pairwise + external) 
double ForceField::beadSysEn(Bead& b, vector < Molecule >& mols, int lenNow){
	double en = 0; 
	//first deal with pairwise potential 
	for(int i = 0; i < (int)mols.size(); i++){ //loop over all beads in box
		for(int j = 0; j < mols[i].size(); j++){
			if(mols[i].monomers[j].id() != b.id()){ //exclude bead's interaction w/itself
				en += pairPot->pairEn(b, mols[i].monomers[j], length, m_dimension); 
			}
		}
	}
	//cout<<" in beadSysEn, calculating pot w/beads selected for trial chain"<<endl; 
	for(int i = 0; i < lenNow; i++){ //loop over beads in trial chain so far
		en += pairPot->pairEn(b, trialChain[i], length, m_dimension);
	}
	if(extDef){ //now interaction w/external potential
		en += extPot->beadEn(b, length); 
	}
	return en; 
}

//grand canonical fun!!!
//lenNow is trial chian length BEFORE this segment is generated 
double ForceField::genVecs(Bead& end, vector < Molecule >& mols, double beta, int lenNow, mt19937& ranGen){
	for(int i = 0; i < m_nVec; i++){ //loop to fill in trial vecs being used to pick this segment
		double bLen = bondPot->ranBond(beta, ranGen); 
		double pos[3]; //array for coords of trial bead
		randSphere(pos, ranGen); //fill it in w/random vec on a unit sphere
		for(int q = 0; q < 3; q++){
                        pos[q] *= bLen; //make it the proper length
                        pos[q] += end.getCoord(0,q); //stick it onto existing chain
                }
		trialVecs[i].setAllCoord(pos); //give the appropriate trial bead these coordinates
		double sysEn = beadSysEn(trialVecs[i], mols, lenNow); 
		wFactors[i] = exp(-1 * beta * sysEn);		
		//trialVecs[i].printXYZ(cout); //print trial vec 
	}
	double wTot = 0; 
	for(int i = 0; i < m_nVec; i++){
		wTot += wFactors[i]; 
	}
	return wTot; //return total w factor for this set of trial vecs
}

bool ForceField::trialIns(vector < Molecule >& mols, double beta, mt19937& ranGen){
	//cout<<" in FF insertion func"<<endl; 
	bool accepted = false;
	insAttempt++;  
	double wAll = 1.0; //wFactor for whole chain will equal prod. of segment wTot values
	//set coords of 1st bead to random point
	double initCoords[3]; 
	for(int i = 0; i < 3; i++) initCoords[i] = (double)ranGen() / ranGen.max() * length[i]; 
	trialChain[0].setAllCoord(initCoords); 
	wAll *= exp(beadSysEn(trialChain[0], mols, 0) * beta * -1); //calculate init bead's w factor
	//cout<<"w for first bead is "<<wAll<<endl; 
	//now loop to make the rest of the chain
	//cout<<" printing trial chain: "<<endl; 
	for(int i = 1; i < m_chainLen; i++){
	//	cout<<"generating trial segment "<<i<<endl; 
		double wVecs = genVecs(trialChain[i-1], mols, beta, i, ranGen); 
		wAll *= wVecs; 
		//now pick a segment in proportion with its w factor
		double r = (double)ranGen() / ranGen.max() * wVecs; 
		int nextBead = 0; //will index in trialVecs of next bead to be added
		double cumW = wFactors[0];
                while(cumW < r){
                        nextBead++;
                        cumW += wFactors[nextBead];
                }
		trialChain[i].setAllCoord(trialVecs[nextBead]); 
		//bead in trial chain now has everything it needs except a valid ID number 	
	}
	//chain has been made; now accept/reject
	double r = (double)ranGen() / ranGen.max(); 
	//normalize wAll 
	wAll /= m_nVec; 
	if(r < m_qT * exp(beta * m_chemPot) * wAll / (mols.size() + 1)){
		//cout<<"  chain accepted! r = "<<r<<endl; 
		accepted = true; 
		insAccept++; 
		//chain accepted! now add an empty mol to  mol vector
		mols.push_back(Molecule()); //gets added to the latter (chain mol) part of the vector
		//now copy data in trialChain beads into new beads in empty mol
		//using copy constructor
		for(int i = 0; i < m_chainLen; i++){
			mols[mols.size() - 1].addBead(trialChain[i]); 
		}
		//add bonds to new mol, ASSUMING LINEAR MOLECULE! this thing cannot build non-linear mols! 
		for(int i = 0; i < mols[mols.size() - 1].size()-1; i++){
			mols[mols.size() - 1].addBond(i, i+1); 
		}
	}
	//at this point, new beads are in w/correct coords. still need proper IDs
	//still need to initialize its energies, update SimBox nParticle fields
	return accepted; 
}

int ForceField::trialDel(vector < Molecule >& mols, int nChain, double beta, mt19937& ranGen){
	delAttempt++; 
	//randomly choose a polymer to try deleting. right now kind of naff b/c we don't have a way of 
        //selecting just from the pool of polymers
	int toDel = floor((double)ranGen() / ranGen.max() * mols.size());
	while(mols[toDel].size() == 1){ //while the one you've selected is not polymer, select another
		toDel = floor((double)rand() / RAND_MAX * mols.size());
	}
	double wAll = 1; 
	wAll *= exp(beadSysEn(mols[toDel].monomers[0], mols, 0) * beta * -1); //w for 1st bead in chain toDel
	//check that trialChain is long enough to retrace the molecule
	if((int)trialChain.size() < mols[toDel].size()){
		cout<<trialChain.size()<<" beads insufficient to retrace chain of length "<<mols[toDel].size()<<"!"<<endl; 
		cout<<"ending program!"<<endl; 
		exit(1); 
	}
	for(int i = 1; i < mols[toDel].size(); i++){ //retrace beads in toDel
		genVecs(mols[toDel].monomers[i-1], mols, beta, i, ranGen);//doesn't include the actual bead, so don't want wtot 	
		wFactors[0] = exp(-1 * beta * beadSysEn(mols[toDel].monomers[i], mols, i)); //now it includes real bead wFactor
		double wSeg = 0;
                for(int j = 0; j < m_nVec; j++){ //now we total wFactors 
                        wSeg += wFactors[j];
                }
                wAll *= wSeg;
	}
	double r = (double)ranGen() / ranGen.max(); 
	//normalize wAll
	wAll /= m_nVec; 
	if(r < mols.size() / (m_qT * exp(beta * m_chemPot) * wAll)){ //deletion successful
		//cout<<" deletion accepted! r = "<<r<<endl; 
		if((int)mols.size() == 1){
			cout<<"About to delete last molecule! Ending program. :-( "<<endl; 
			exit(1); 
		}
		delAccept++; 
		//delete molecule in all energy arrays
		pairPot->delMol(mols, toDel); 
		bondPot->delMol(toDel); 
		if(extDef) extPot->delMol(mols, toDel); 
		//delete molecule from mol vector. maybe have SimBox do this?  
		return toDel; 
	}
	return -1; //signals unsuccessful deletion

}

//does all energy initializing for a new mol. requires that all IDs are properly assigned. 
void ForceField::addMol(vector < Molecule >& mols){
//	cout<<"in FF, initializing new mol!"<<endl; 
	pairPot->initMol(mols, length, m_dimension); 
//	cout<<"pairPot initialized"<<endl; 
	bondPot->initMol(mols, length, m_dimension);
//	cout<<"bondPot initialized"<<endl; 
	if(extDef){
		 extPot->initMol(mols, length, m_dimension);
//		cout<<"extPot initialized"<<endl; 
	}
}

void ForceField::enStats(ofstream& energies){
	if(m_gc){
		cout<<setprecision(4)<<(double)insAccept/insAttempt*100<<" \% of trial insertions accepted, and "; 
		cout<<setprecision(4)<<(double)delAccept/delAttempt*100<<" \% of trial deletions accepted"; 
	}
	cout<<" "<<endl; 
	energies<<"#Pressure tensor components: "<<endl;
	//cout<<"enCoutner = "<<enCounter<<endl; 
        for(int i = 0; i < 3; i++){
                double deltaV = getDeltaV(i);
                double expPlus = pPlus[i] / enCounter;
                double expMinus = pMinus[i] / enCounter;
                double p = (1/m_beta) * (1/deltaV * log(expPlus) + -1/deltaV * log(expMinus)) / 2;
		p *= 16.02; //convert from weird units to units of 10^10 Pa 
                energies<<"#  "<<p<<endl;
        }
}

int ForceField::gcEvery(){
	return m_gcEvery; 
}

bool ForceField::gc(){
	return m_gc; 
}

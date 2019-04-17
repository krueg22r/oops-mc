#include <stdlib.h>
#include <cmath> 
#include <random> 
#include "Bead.h"
using namespace std; 

//bead-bead dist func that takes bead ptrs. from MC12 & earlier
double getDist(Bead* b1, Bead* b2, double length[], int dimension){
        double dist = 0;
        for(int i = 0; i < 3; i++){
                double di = b1->getCoord(1, i) - b2->getCoord(1,i);
                if(i < dimension) di -= length[i] * round(di / length[i]); //periodic boundary conditions 
                dist += (di * di);
        }
        return sqrt(dist);
}

//bead-bead dist func that takes refs to beads. this is the MC13 way. 
double getDist(Bead& b1, Bead& b2, double length[], int dimension){
	double dist = 0;
        for(int i = 0; i < 3; i++){
		double di = b1.getCoord(1, i) - b2.getCoord(1,i);
                if(i < dimension) di -= length[i] * round(di / length[i]); //periodic boundary conditions 
                dist += (di * di);
	}
	return sqrt(dist);
}


//fills in a vector with random point on unit sphere
void randSphere(double vec[], mt19937& ranGen){
	        double ransq = 2;
        double r1, r2;
        while(ransq >= 1){
		r1= 1 - 2 * ((double)ranGen() / ranGen.max()); 
		r2 = 1 - 2 * ((double)ranGen() / ranGen.max());
                ransq = r1*r1 + r2 * r2;
        }
        double ranh = 2 * sqrt(1 - ransq);
	vec[0] = r1 * ranh; 
	vec[1] = r2 * ranh; 
	vec[2] = (1 - 2*ransq); 
}

//returns rand deviate from gaussian distribution with mean 0 and stdev 1. 
//multiplied /shifted to give the requested distribution
double gasdev(double mean, double stdev, mt19937& ranGen){
        static int iset = 0;
        static double gset;
        double fac, rsq, v1, v2;
        if (iset == 0){
                do{
                        v1=2.0*((double)ranGen() / ranGen.max()) - 1.0;
                        v2=2.0*((double)ranGen() / ranGen.max()) - 1.0;
                        rsq = v1*v1 + v2*v2;
                }while (rsq >= 1.0 || rsq == 0.0 );
                fac = sqrt(-2.0 * log(rsq)/rsq);

                gset = v1*fac;
                iset = 1;
                return v2*fac * stdev + mean;
        }else{
                iset = 0;
                return gset * stdev + mean;
        }
	return 0; 
}

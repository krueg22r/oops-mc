//Useful small functions used multiple places in program. 
#ifndef SMALL_H
#define SMALL_H
#include <cmath>
#include <stdlib.h> 
#include <random> 
#include "ran1.h"
#include "Bead.h"
#include<Eigen/Dense>
using namespace Eigen; 
using namespace std; 

double getDist(Bead&, Bead&, double[], int); 
double getDist(Bead*, Bead*, double[], int);
void randSphere(double[], mt19937&); 
double gasdev(double, double, mt19937&);
double pGauss(double, double, double);
#endif 

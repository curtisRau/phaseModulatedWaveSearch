//
//  besselFunction.cpp
//  phaseModulatedGravitationalWaveFinder
//
//  Created by Curtis Rau on 3/22/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include <iostream>
#include "besselFunction.hpp"
#include <math.h>

using namespace std;


// Calculate bessel function of first kind using integral
// (1\pi) int_{0}^{pi} Cos[ z Sin (theta) - n theta ] d theta
double BesselJ (unsigned int n, double z) {
    unsigned int Npts = 100 * (n + 1) * floor((2 * z / M_PI) + 1);
    double dx = M_PI / Npts;
    double sum = 0.0;
    
    for (double x = 0.0; x < M_PI; x += dx) {
        sum += cos(z * sin(x) - n * x);
    }
    
    return sum / Npts;
}


//    numOfPts is the number of points returned, which includes endpoints.
double* besselEqationSolve (unsigned int n, double xMax, unsigned int numOfPts) {
    double       dx          = xMax / static_cast<double>(numOfPts - 1.0);
    double       dx2         = dx * dx;
    unsigned int n2          = n * n;
    
    
    // Generate the "a" array:
    double* a = new double[numOfPts - 2];
    for (unsigned int i = 0; i < (numOfPts - 2); i++) {
        a[i] = (i + 1.0) * (i + 0.5);
    }
    
    // Generate the "b" array:
    double* b = new double[numOfPts - 2];
    for (unsigned int i = 0; i < (numOfPts - 2); i++) {
        b[i] = (i + 1.0) * (i + 1.0) * (dx2 - 2.0) - n2;
    }
    
    // Generate the "c" array:
    double* c = new double[numOfPts - 2];
    for (unsigned int i = 0; i < (numOfPts - 2); i++) {
        c[i] = (i + 1.0) * (i + 1.5);
    }
    
    // Generate the solution array:
    double* J = new double[numOfPts];
    
    // Boundary condition at x = 0:
    if (n == 0) {
        J[0] = 1.0;
    } else {
        J[0] = 0.0;
    }
    
    // Boundary condition at x = xMax:
    // This section should make a decision weather to use asymptoticBesselJ or BesselJ
    J[numOfPts - 1] = BesselJ(n, xMax);
    
    // Generate the source array:
    double* S = new double[numOfPts - 2];
    S[0] = - a[0] * J[0];
    for (unsigned int i = 1; i < (numOfPts - 3); i++) {
        S[i] = 0.0;
    }
    S[numOfPts - 3] = - c[numOfPts - 3] * J[numOfPts - 1];
    
    
    // Forward substitution
    for (unsigned int i = 1; i <= (numOfPts - 3); i++) {
        b[i] = b[i] - (a[i]/b[i-1])*c[i-1];
        S[i] = S[i] - (a[i]/b[i-1])*S[i-1];
    }
    
    J[numOfPts - 2] = S[numOfPts - 3] / b[numOfPts - 3];
    
    // Backwards substitution
    for (int i = (numOfPts - 4); i >= 0; i--) {
        S[i]   = S[i] - (c[i]/b[i+1])*S[i+1];
        J[i+1] = S[i]/b[i];
    }
    
    
    
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] S;
    
    return J;
}

// This function returns a matrix where the rows contain "J_n(x)" for differing values of n,
// and the columns contain differing values of x.  "GammaMax" is the maximum value of Gamma,
// that is returned, and 0 is the minimum value of Gamma.  "nGamma" is the number of Gamma
// returned including Gamma = 0 and Gamma = GammaMax.  "nMax" is the maximum order Bessel
// function that is returned.
double** generateBesselJMatrix (double GammaMax, unsigned int nGamma, unsigned int nMax) {
    nMax += 1;      // We want to include 0 ... nMax, so there will be nMax+1 pts.
    
    cout << "Gennerating Bessel Array.  This will take " << sizeof(double) * nGamma * nMax * pow(10.0, -6.0) << "MB of disk space.  (Plus an extra 25% for delimeters)" << endl;
    
    double** J = new double* [nGamma];
    for (unsigned int i = 0; i < nGamma; i++) {
        J[i] = new double [nMax];
    }
    
    double* column = new double [nGamma];
    for (unsigned int i = 0; i < nMax; i++) {
        column = besselEqationSolve(i, GammaMax, nGamma);
        for (unsigned int j = 0; j < nGamma; j++) {
            J[j][i] = column[j];
        }
    }
    
    return J;
}
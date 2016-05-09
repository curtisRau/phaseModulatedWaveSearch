//
//  besselFunction.hpp
//  phaseModulatedGravitationalWaveFinder
//
//  Created by Curtis Rau on 3/22/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#ifndef besselFunction_hpp
#define besselFunction_hpp

#include <stdio.h>

#endif /* besselFunction_hpp */

//double asymptoticBesselJ (unsigned int n, double x);
double BesselJ (unsigned int n, double z);
double* besselEqationSolve (unsigned int n, double xMax, unsigned int numOfPts);
double** generateBesselJMatrix (double GammaMax, unsigned int nGamma, unsigned int nMax);
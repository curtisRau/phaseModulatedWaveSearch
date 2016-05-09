//  main.cpp
//  phaseModulatedGravitationalWaveFinder
//
//  Created by Curtis Rau on 1/31/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include <iostream>
#include <complex>              // For working with complex numbers.
#include <fstream>              // for working with files.
#include <algorithm>            // std::min
#include "stdlib.h"             // Necessary for rand?
//#include "omp.h"                // For Open MP
#include "besselFunction.hpp"
#include "functions.hpp"

using namespace std;
typedef numeric_limits<double> dbl;


//// F is an array of complex points which represents a discrete fourier transform.
//// fMax is the maximum frequency in F.
//// df is the step size in frequency of F.
//// besselArray is an array of bessel functions "J_n(x)" all evaluated at the same
////      value of "x", but different "n".
//// omega is the value of "omega0" we are trying.
//// f1 CANNOT BE 0 BECAUSE DIVISION BY 0!
//// fMax CANNOT BE LESS THAN f!
//double gustafsonAlgorithm (complex<double>* F, double nyquistF, double df, double* besselArray, int besselArrayLength, double f1, double f, double phi1) {
//    
//    unsigned int NFS = floor(2.0 * nyquistF / df) + 1;      // Number of points in fourier series.
//    unsigned int NF0 = 0;                       // F[NF0] = F[f = 0.0 Hz]
////        unsigned int NF0 = (NFS - 1) / 2;                       // F[NF0] = F[f = 0.0 Hz]
//    
//    if (f > nyquistF) {
//        cout << "Fatal Error in calling gustafsonAlgorithm: f > fMax" << endl;
//        return 0.0;
//    }
//    
//    if (f1 == 0.0) {
//        cout << "Fatal Error in calling gustafsonAlgorithm: f1 = 0 => division by zero" << endl;
//    }
//    
//    const complex<double> I   (0.0, -(phi1 + M_PI_2));
//    const complex<double> J   (0.0,  (phi1 + M_PI_2));
//          complex<double> sum (0.0, 0.0);
//    
//    unsigned int nMax = floor((nyquistF - f) / f1);
//    if (nMax > besselArrayLength) {
//        cout << "Non-Fatal Error in calling gustafsonAlgorithm: nMax > besselArrayLength" << endl;
//        nMax = besselArrayLength;
//    }
//    for (int n = 0; n < nMax; n++) {
//        sum += exp(I * static_cast<double>(n)) * besselArray[n] * F[NF0 + static_cast<int>((f + n * f1) / df)];
//    }
//    
//    
////    nMax = floor((nyquistF + f) / f1);
//    nMax = floor(f / f1);
//    if (nMax > besselArrayLength) {
//        cout << "Non-Fatal Error in calling gustafsonAlgorithm: nMax > besselArrayLength" << endl;
//        nMax = besselArrayLength;
//    }
//    for (int n = 0; n < nMax; n++) {
//        sum += exp(J * static_cast<double>(n)) * pow(-1.0,n) * besselArray[n] * F[NF0 + static_cast<int>((f - n * f1) / df)];
//    }
//    
//    return 2.0 * abs(sum);
//}


unsigned int calcNumberSidebands (double Gamma) {
    return max(static_cast<unsigned int>(1),
               static_cast<unsigned int>(
                                     2.0 * (Gamma - log(Gamma)/log(M_PI))
                                     ));
}

// Fold a frequency.
// "f"   is the frequency to be folded
// "fnq" is the nyquist or folding frequency
// "f"   is folded between 0 and fnq
unsigned int foldIntFrequency (unsigned int f, unsigned int fnq) {
    unsigned int ffolded = f;
    if (f > fnq) {
        unsigned int n = static_cast<unsigned int>(ceil(static_cast<double>(f)/static_cast<double>(fnq))) % 2;
        switch (n) {
            case 0:                         // Even folding zone
                ffolded = fnq - f % fnq;    // Mirror Frequency
                break;
            case 1:                         // Odd folding zone
                ffolded = f % fnq;          // Translate Frequency
                break;
        }
    }
    return ffolded;
}

//// F is an array of complex points which represents a discrete fourier transform.
//// fMax is the maximum frequency in F.
//// df is the step size in frequency of F.
//// besselArray is an array of bessel functions "J_n(x)" all evaluated at the same
////      value of "x", but different "n".
//// omega is the value of "omega0" we are trying.
//// f1 CANNOT BE 0 BECAUSE DIVISION BY 0!
//// fMax CANNOT BE LESS THAN f!
//double gustafsonAlgorithm (complex<double>* F, double nyquistF, double df, double* besselArray, int besselArrayLength, unsigned int nSidebandsSum,double f1, double f, double phi1) {
//    
//    const complex<double> I   (0.0, -(phi1 + M_PI_2));
//    const complex<double> J   (0.0,  (phi1 + M_PI_2));
//    complex<double> sum (0.0, 0.0);
//    
//    unsigned int nMax = min(nSidebandsSum, static_cast<unsigned int>(floor((nyquistF - f) / f1)));
//    if (nMax > besselArrayLength) {
//        cout << "Non-Fatal Error in calling gustafsonAlgorithm: nMax > besselArrayLength" << endl;
//        nMax = besselArrayLength;
//    }
//    
//    
//    for (int n = 0; n < nMax; n++) {
//        sum += exp(I * static_cast<double>(n)) * besselArray[n] * F[static_cast<int>((f + n * f1) / df)];
//    }
//    
//    
//    nMax = floor(f / f1);
//    if (nMax > besselArrayLength) {
//        cout << "Non-Fatal Error in calling gustafsonAlgorithm: nMax > besselArrayLength" << endl;
//        nMax = besselArrayLength;
//    }
//    for (int n = 0; n < nMax; n++) {
//        sum += exp(J * static_cast<double>(n)) * pow(-1.0,n) * besselArray[n] * F[static_cast<int>((f - n * f1) / df)];
//    }
//    
//    return 2.0 * abs(sum);
//}



// F is an array of complex points which represents a discrete fourier transform.
// fMax is the maximum frequency in F.
// df is the step size in frequency of F.
// besselArray is an array of bessel functions "J_n(x)" all evaluated at the same
//      value of "x", but different "n".
// omega is the value of "omega0" we are trying.
// f1 CANNOT BE 0 BECAUSE DIVISION BY 0!
// fMax CANNOT BE LESS THAN f!
double gustafsonAlgorithmNEW (complex<double>* F, unsigned int nyquistF, double* besselArray, int besselArrayLength, unsigned int nSidebandsSum, unsigned int f0, unsigned int f1, double phi1) {
    
    const complex<double> I   (0.0, -(phi1 + M_PI_2));
    const complex<double> J   (0.0,  (phi1 + M_PI_2));
    complex<double> sum (0.0, 0.0);

    
    for (int n = -nSidebandsSum; n <= nSidebandsSum; n++) {
        sum += exp(I * static_cast<double>(n)) * besselArray[n] * F[foldIntFrequency(f0 + n*f1, nyquistF)];
    }
    
//    for (int n = 0; n < nMax; n++) {
//        sum += exp(J * static_cast<double>(n)) * pow(-1.0,n) * besselArray[n] * F[static_cast<int>((f - n * f1) / df)];
//    }
    
    return 2.0 * abs(sum);
}


double* generateData (double t0, double t1, unsigned int N, double Gamma, double f0, double f1, double phi1) {
    double dt = (t1-t0) / static_cast<double> (N);
    double* data = new double[N];
    
    f0 *= 2.0 * M_PI;   // Convert to angular frequency.
    f1 *= 2.0 * M_PI;   // Convert to angular frequency.
    
    double t;
    for (unsigned int i = 0; i < N; i++) {
        t = t0 + i * dt;
        data[i] = cos(f0 * t + Gamma * cos(f1 * t) + phi1);
    }
    
    return data;
}



void addNoise (double* data, unsigned int lengthData, double range) {
    for (unsigned int i = 0; i < lengthData; i++) {
        data[i] += range * ((rand() % 1000) - 500) / 1000.0;
    }
}



//// NTS         = The number of points in the time series.
//// nyquistFreq = The Nyquist Frequency of the time series.
//// df          = step size in frequency.
//complex<double>* fourierTransform (double* timeSeries, double startTime, unsigned int NTS, double nyquistFreq, double df) {
//    unsigned int NFS = floor(2.0 * nyquistFreq / df) + 1;                               // Number of points to return in the Fourier Series
//    const complex<double> I (0.0, -2.0 * M_PI);
//    complex<double>* data = new complex<double> [NFS];
//    
//    double t;
//    double f;
//    for (unsigned int i = 0; i < NFS; i++) {
//        f = - nyquistFreq + i * df;
//        for (unsigned int j = 0; j < NTS; j++) {
//            t = startTime + j * 0.5 / nyquistFreq;
//            data[i] += timeSeries[j] * exp(I * f * t);
//        }
//        data[i] /= sqrt(NTS);
//    }
//    
//    return data;
//}



// NTS         = The number of points in the time series.
// nyquistFreq = The Nyquist Frequency of the time series.
// df          = step size in frequency.
complex<double>* fourierTransform (double* timeSeries, double startTime, unsigned int NTS, double nyquistFreq, double df) {
    unsigned int NFS = floor(2.0 * nyquistFreq / df) + 1;                               // Number of points to return in the Fourier Series
    const complex<double> I (0.0, -2.0 * M_PI);
    complex<double>* data = new complex<double> [NFS];
    
    double t;
    double f;
    for (unsigned int i = 0; i < NFS; i++) {
        f = i * df;
        for (unsigned int j = 0; j < NTS; j++) {
            t = startTime + j * 0.5 / nyquistFreq;
            data[i] += timeSeries[j] * exp(I * f * t);
        }
        data[i] /= sqrt(NTS);
    }
    
    return data;
}




int main(int argc, const char * argv[]) {
    cout.precision(dbl::max_digits10);
    
    // Program Parameters:
    const char* filename    = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/out.csv";       // File for storing output data.
    double sampleRate       = 300.0;        // [Hz] -- The sampleing rate of the time series.
    double initialTime      = 0.0;          // [s]  -- The start time of the trial.  First pt. in time series will be here.
    double finalTime        = 10.0;         // [s]  -- The final time of the time series.
    double df               = 0.01;         // [Hz] -- The step size in frequency of the fourier transform.
    
    // Parameters the User will input when the fourier transform is already calculated:
    // File containing the fourier transform data
    
    
    // Search Parameters:
//    unsigned int nSidebands = 500;          // [ ]   -- The maximum number of sidebands we want to sum over on each side of the carrier.
    double gammaMax         = 20.0;         // [?]   -- The maximum value for Gamma over which we want to perform our search.
    unsigned int nSidebands = calcNumberSidebands(gammaMax);  // [ ]   -- The maximum number of sidebands we want to sum over on each side
    double dGamma           = 0.1;         // [?]   -- The step size in Gamma which we take in performing our search.
    double dPhi1            = 0.01;         // [rad] -- The step size in Phi1 which we take in performing our search.
    
    // Data Parameters:
    double noiseLevel = 10.0;
    double Gamma = 10.0;
    double f0    = 20.0;
    double f1    = 6.0;
    double phi1  = 0.0;
    
    // Derived Parameters:
    double duration     = finalTime - initialTime;          // [s]  -- The duration of the time series.
    unsigned int NTS    = floor(duration * sampleRate);     // [ ]  -- The number of sample points in the time series.
    double nyquistF     = 0.5 * sampleRate;                 // [Hz] -- The Nyquist Frequency.
//    double df = sampleRate / static_cast<double>(NTS-1);
//    unsigned int NFS = NTS;
//    unsigned int NFS    = floor(2.0 * nyquistF / df) + 1;   // [ ]  -- The number of points in the Fourier Series.
    unsigned int NFS    = floor(nyquistF / df) + 1;         // [ ]  -- The number of points in the Fourier Series.
    unsigned int nGamma = floor(gammaMax / dGamma);         // [ ]  -- The number of steps in Gamma we will take.
    unsigned int nPhi1  = floor(2.0 * M_PI / dPhi1);        // [ ]  -- The number of steps in Phi1 we will take.
    
    
    // Print Useful Parameters:
    cout << "Sample Rate       = " << sampleRate << " Hz"  << endl;
    cout << "Nyquist Frequency = " << nyquistF   << " Hz"  << endl;
    cout << "df                = " << df         << " Hz"  << endl;
    
    cout << "Gamma Max         = " << gammaMax   << " rad" << endl;
    cout << "# of Sidebands    = " << nSidebands << " "    << endl;
    
    // TEST THE FOURIER TRANSFORM FUNCTION:
    if (false) {
        clock_t begin_time = clock();
        double* sol = generateData(initialTime, finalTime, NTS, Gamma, f0, f1, phi1);
        
        addNoise(sol, NTS, noiseLevel);
        complex<double>* trans = fourierTransform(sol, initialTime, NTS, nyquistF, df);
        cout << "Generating Data and Calculating Fourier Transform took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        begin_time = clock();

        double* out = new double [NFS];
        for (unsigned int i = 0; i < NFS; i++) {
            out[i] = abs(trans[i]);
        }
        
        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/out2.csv";
        saveArray4Mathematica(filename, out, NFS);
        
        cout << "Writing the solution to disk took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
    }
    
    

// CALCULATING BESSEL FUNCTIONS FOR DIFFERENT GAMMA BUT SAME N:
    if (false) {
        double* sol = besselEqationSolve(0, gammaMax, nGamma);
        
        saveArray4Mathematica(filename, sol, nGamma);

    }
    
    
    
// CALCULATING HOW MUCH FASTER THE THOMAS ALGORITHM IS THAN POINT BY POINT NUMERIC INTEGRATION
    if (false) {
        clock_t begin_time = clock();
        besselEqationSolve(3, 20.0, 500);
        cout << "Diffeq method time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
        
        begin_time = clock();
        double dG = 20.0 / 500.0;
        for (unsigned int i = 0; i <= 500; i++) {
            BesselJ(3, dG * i);
        }
        cout << "Standard method time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
    }
    
    
    
// IMPLEMENTING THE GUSTAFSON ALGORITHM
    if (true) {
        
        clock_t begin_time = clock();
        double* sol = generateData(initialTime, finalTime, NTS, Gamma, f0, f1, phi1);
        cout << "Generating Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/timeSeries.csv";
        saveArray4Mathematica(filename, sol, NTS);
        cout << "Saving Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        complex<double>* trans = fourierTransform(sol, initialTime, NTS, nyquistF, df);
        cout << "Calculating Fourier Transform of Noisy Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        double* absTrans = new double [NFS];
        for (unsigned int i = 0; i < NFS; i++) {
            absTrans[i] = abs(trans[i]);
        }
        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/frequencySeries.csv";
        saveArray4Mathematica(filename, absTrans, NFS);
        cout << "Saving Fourier Transform of Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        addNoise(sol, NTS, noiseLevel);
        cout << "Adding Noise took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/timeSeriesNoise.csv";
        saveArray4Mathematica(filename, sol, NTS);
        cout << "Saving Noisy Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        trans = fourierTransform(sol, initialTime, NTS, nyquistF, df);
        cout << "Calculating Fourier Transform of Noisy Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
        begin_time = clock();
        for (unsigned int i = 0; i < NFS; i++) {
            absTrans[i] = abs(trans[i]);
        }
        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/frequencySeriesNoise.csv";
        saveArray4Mathematica(filename, absTrans, NFS);
        cout << "Saving Fourier Transform of Noisy Data took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";

        
        begin_time = clock();
        double** J = generateBesselJMatrix(gammaMax, nGamma, nSidebands);
        cout << "Generating Bessel Matrix took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";

        
        begin_time = clock();
        double* out = new double [NFS];
        cout << "Allocating Memory for Solution Matrix took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";

        
// Gustafson Algorithm Steps in omega0
//        begin_time = clock();
//        for (unsigned int i = 0; i < NFS; i++) {
//            out[i] = gustafsonAlgorithm(trans,
//                                        nyquistF,
//                                        df,
//                                        J[static_cast<int>(Gamma / dGamma)],
//                                        nSidebands,
//                                        calcNumberSidebands(Gamma),
//                                        f1,
//                                        df * i,
////                                        df * i - nyquistF,
//                                        phi1);
//        }
        begin_time = clock();
        for (unsigned int i = 0; i < NFS; i++) {
            out[i] = gustafsonAlgorithmNEW(trans,
                                           NFS-1,
                                           J[static_cast<int>(Gamma / dGamma)],
                                           nSidebands,
                                           nSidebands, //calcNumberSidebands(Gamma),
                                           i,
                                           600,
                                           phi1);
        }
        cout << "Performing the Gustafson Algorithm over omega0 took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";

        
        begin_time = clock();
        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/outOmega0.csv";
        saveArray4Mathematica(filename, out, NFS);
        cout << "Writing the solution to disk took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
// Gustafson Algorithm Steps in omega1
//        begin_time = clock();
////        for (unsigned int i = 0; i < (NFS-1)/2; i++) {
//        for (unsigned int i = 0; i < NFS; i++) {
//            out[i] = gustafsonAlgorithm(trans,
//                                        nyquistF,
//                                        df,
//                                        J[static_cast<int>(Gamma / dGamma)],
//                                        nSidebands,
//                                        calcNumberSidebands(Gamma),
//                                        df * i,
//                                        f0,
//                                        phi1);
//        }
//        cout << "Performing the Gustafson Algorithm over omega1 took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
//        
//        
//        begin_time = clock();
//        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/outOmega1.csv";
////        saveArray4Mathematica(filename, out, (NFS-1)/2);
//        saveArray4Mathematica(filename, out, NFS);
//        cout << "Writing the solution to disk took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
        
        
// Steps in Gamma:
//        begin_time = clock();
//        double* out2 = new double [nGamma];
//        for (unsigned int i = 0; i < nGamma; i++) {
//            out2[i] = gustafsonAlgorithm(trans,
//                                        nyquistF,
//                                        df,
//                                        J[i],
//                                        nSidebands,
//                                        calcNumberSidebands(dGamma * i),
//                                        f1,
//                                        f0,
//                                        phi1);
//        }
//        cout << "Performing the Gustafson Algorithm over Gamma took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
//        
//        
//        begin_time = clock();
//        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/outGamma.csv";
//        saveArray4Mathematica(filename, out2, nGamma);
//        cout << "Writing the solution to disk took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
//        
//        
//        begin_time = clock();
//        double* out3 = new double [nPhi1];
//        for (unsigned int i = 0; i < nPhi1; i++) {
//            out3[i] = gustafsonAlgorithm(trans,
//                                         nyquistF,
//                                         df,
//                                         J[static_cast<int>(Gamma / dGamma)],
//                                         nSidebands,
//                                         calcNumberSidebands(Gamma),
//                                         f1,
//                                         f0,
//                                         dPhi1 * i + M_PI);
//        }
//        cout << "Performing the Gustafson Algorithm over Phi1 took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
//        
//        
//        begin_time = clock();
//        filename = "/Volumes/userFilesPartition/Users/curtisrau/Desktop/output/outPhi1.csv";
//        saveArray4Mathematica(filename, out3, nPhi1);
//        cout << "Writing the solution to disk took = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s \r";
//    }
    
    }
    return 0;
}

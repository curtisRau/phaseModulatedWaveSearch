//
//  functions.cpp
//  phaseModulatedGravitationalWaveFinder
//
//  Created by Curtis Rau on 4/15/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include "functions.hpp"
#include <iostream>
#include <fstream>              // for working with files.


// Can be vectorized.
void saveArray4Mathematica (const char* filename, double* array, unsigned int arraySize) {
    std::ofstream outputFile;
    outputFile.open(filename, std::ios::out | std::ios::trunc);         // Open a file for output and overwrite current content if it exists.
    
    if (outputFile.is_open()) {                                         // If the file is open...
        outputFile << array[0];
        for (unsigned int i = 1; i < arraySize; i++) {
            outputFile << "\t" << array[i];
        }
    } else {
        std::cout << "File '" << filename << "' did not open /r";
    }
    outputFile.close();
}


void saveMatrix4Mathematica (const char* filename, double** matrix, unsigned int matrixSizeM, unsigned int matrixSizeN) {
    std::ofstream outputFile;
    outputFile.open(filename, std::ios::out | std::ios::trunc);         // Open a file for output and overwrite current content if it exists.
    
    if (outputFile.is_open()) {                                         // If the file is open...
        for (unsigned int i = 0; i < (matrixSizeM - 1); i++) {
            outputFile << matrix[i][0];
            for (unsigned int j = 1; j < matrixSizeN; j++) {
                outputFile << "\t" << matrix[i][j];
            }
            outputFile << "\r";
        }
        outputFile << matrix[matrixSizeM - 1][0];
        for (unsigned int j = 1; j < matrixSizeN; j++) {
            outputFile << "\t" << matrix[matrixSizeM - 1][j];
        }
    } else {
        std::cout << "File '" << filename << "' did not open /r";
    }
    outputFile.close();
}
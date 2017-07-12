/**
 *   \file BlasHeader.cpp
 *   FiniteDiffGSL
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

//#include <queue>
#include <iostream>
#include <sstream>
#include "ModelRunning.h"
#include "BlasHeader.h"
#include "CVodeHelpers.h"

using namespace std;

extern "C" double pyEntry(double *pIn) {
    double nulll[3];
    
    return calcError(Param(pIn), nulll);
}

extern "C" void pyEntryVec(double *pIn, double *pOut, int n) {
    double nulll[3];
    
    for (int ii = 0; ii < n; ii++) {
        pOut[ii] = calcError(Param(&pIn[ii*Nparams]), nulll);
    }
}


extern "C" void U87pyEntryVec(double *pIn, double *pOut, int n) {
    double nulll[3];
    
    for (int ii = 0; ii < n; ii++) {
        pOut[ii] = U87calcError(Param(&pIn[ii*Nparams]), nulll);
    }
}


extern "C" int calcProfileMatlab(double *pYData, double *totData, double *surfData, double *speciesData, double *params, double *tps, unsigned int nTps, double GasStim, double *convFac) {
    struct rates pInS = Param(params);
    
    try {
        calcProfileSet (pYData, totData, surfData, speciesData, tps, &pInS, nTps, GasStim, convFac);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int diffCalc(double *dataPtrTemp, double *dataPtrpY, double *dataPtrTot, double *dataPtrSurf, double *GasIn, unsigned int gridIn, double *params, double *tps, unsigned int nTps, double *dIn) {
    
    try {
        diffusionSolution(dataPtrTemp, GasIn, gridIn, params, tps, nTps, dIn);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    // Create the parameter structure
    struct rates pInS;
    pInS = Param(params);
    
    N_Vector state = N_VNew_Serial(Nspecies);
    
    for (size_t time = 0; time < nTps; time++) {
        for (size_t gridP = 0; gridP < gridIn; gridP++) {
            for (size_t spec = 0; spec < Nspecies; spec++) {
                Ith(state,spec) = dataPtrTemp[time*gridIn*Nspecies + spec*gridIn + gridP];
            }
            
            dataPtrpY[time*gridIn + gridP] = pYcalc(state, &pInS);
            dataPtrTot[time*gridIn + gridP] = totCalc(state, &pInS);
            dataPtrSurf[time*gridIn + gridP] = surfCalc(state);
        }
    }
    
    
    N_VDestroy_Serial(state);
    return 0;
}

extern "C" int getNspecies () {
    return Nspecies;
}

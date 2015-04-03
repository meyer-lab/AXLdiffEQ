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

extern "C" int calcProfileMatlab(double *pYData, double *totData, double *surfData, double *params, double *tps, unsigned int nTps, double GasStim) {
    struct rates pInS = Param(params);
    
    try {
        calcProfileSet (pYData, totData, surfData, tps, &pInS, nTps, GasStim);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS(double *dataPtr, double *GasIn, unsigned int gridIn, double *params, double *tps, unsigned int nTps, double *dIn) {
    
    try {
        diffusionSolution(dataPtr, GasIn, gridIn, params, tps, nTps, dIn);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS_pY(double *dataPtr, double *GasIn, unsigned int gridIn, double *params, double *tps, unsigned int nTps, double *dIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps*Nspecies];
    
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
            
            if (frac == 0) {
                dataPtr[time*gridIn + gridP] = pYcalc(state, &pInS);
            } else if (frac == 1) {
                dataPtr[time*gridIn + gridP] = pYcalc(state, &pInS) / totCalc(state, &pInS);
            } else {
                dataPtr[time*gridIn + gridP] = totCalc(state, &pInS);
            }
        }
    }
    
    N_VDestroy_Serial(state);
    return 0;
}

extern "C" int matlabDiffTPS_pYavg(double *dataPtr, double *GasIn, unsigned int gridIn, double *params, double *tps, unsigned int nTps, double *dIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps];
    
    int flag = matlabDiffTPS_pY(dataPtrTemp, GasIn, gridIn, params, tps, nTps, dIn, frac);
    if (flag == 1) return(1);
    
    double summ;
    
    for (size_t time = 0; time < nTps; time++) {
        summ = 0;
        
        for (size_t gridP = 0; gridP < gridIn; gridP++) {
            summ += dataPtrTemp[time*gridIn + gridP]*((double) gridP);
        }
        
        
        dataPtr[time] = 2*summ/gridIn/gridIn;
    }
    
    return 0;
}

extern "C" int getNspecies () {
    return Nspecies;
}

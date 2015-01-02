/**
 *   \file BlasHeader.cpp
 *   FiniteDiffGSL
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

//#include <queue>
#include <iostream>
#include <sstream>
#include "CVodeHelpers.h"
#include "ReactionCode.h"
#include "ModelRunning.h"
#include "BlasHeader.h"

using namespace std;

extern "C" double pyEntry(double *pIn) {
    return calcError(Param(pIn));
}

extern "C" double pyEntryFull(double *pIn) {
    return calcErrorFull(Param(pIn));
}


extern "C" int calcProfileMatlab(double *dataPtr, double *params, double *tps, int nTps, double GasStim, int frac) {
    struct rates pInS = Param(params);
    
    try {
        calcProfileSet (dataPtr, tps, &pInS, nTps, GasStim, frac);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS(double *dataPtr, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn) {
    
    try {
        diffusionSolution(dataPtr, GasIn, gridIn, params, tps, nTps, dIn, endoImpairIn, degImpairIn);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS_pY(double *dataPtr, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps*Nspecies];
    
    int flag = matlabDiffTPS(dataPtrTemp, GasIn, gridIn, params, tps, nTps, dIn, endoImpairIn, degImpairIn);
    if (flag == 1) return(1);
    
    // Create the parameter structure
    struct rates pInS;
    pInS = Param(params);
    

    N_Vector state = N_VNew_Serial(Nspecies);
    
    for (size_t time = 0; time < (size_t) abs(nTps); time++) {
        for (size_t gridP = 0; gridP < (size_t) abs(gridIn); gridP++) {
            for (size_t spec = 0; spec < Nspecies; spec++) {
                Ith(state,spec) = dataPtrTemp[time*(((size_t) abs(gridIn))*Nspecies) + spec*((size_t) abs(gridIn)) + gridP];
            }
            
            if (frac == 0) {
                dataPtr[time*((size_t) abs(gridIn)) + gridP] = pYcalc(state, &pInS);
            } else if (frac == 1) {
                dataPtr[time*((size_t) abs(gridIn)) + gridP] = pYcalc(state, &pInS) / totCalc(state, &pInS);
            } else {
                dataPtr[time*((size_t) abs(gridIn)) + gridP] = totCalc(state, &pInS);
            }
        }
    }
    
    N_VDestroy_Serial(state);
    return 0;
}

extern "C" double pyDiffTPS_Activation(double GasC, double *params, double tp, double *dIn) {
    double tps[2] = {0, tp};
    const int gridIn = 30;
    double GasIn[gridIn];
    
    for (size_t ii = 0; ii < gridIn; ii++) {
        if (ii < ((double) gridIn / 15)) {
            GasIn[ii] = GasC;
        } else {
            GasIn[ii] = 0;
        }
    }
    
    double dataPtr[2];
    
    
    int flag = matlabDiffTPS_pYavg(dataPtr, GasIn, gridIn, params, tps, 2, dIn, 1, 1, 0);
    
    if (flag == -1) {
        return -100;
    } else {
        return dataPtr[1] / dataPtr[0];
    }
}




extern "C" int matlabDiffTPS_pYavg(double *dataPtr, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps];
    
    int flag = matlabDiffTPS_pY(dataPtrTemp, GasIn, gridIn, params, tps, nTps, dIn, endoImpairIn, degImpairIn, frac);
    if (flag == 1) return(1);
    
    double summ;
    
    for (size_t time = 0; time < (size_t) abs(nTps); time++) {
        summ = 0;
        
        for (size_t gridP = 0; gridP < (size_t) abs(gridIn); gridP++) {
            summ += dataPtrTemp[time*((size_t) abs(gridIn)) + gridP]*((double) gridP);
        }
        
        
        dataPtr[time] = 2*summ/gridIn/gridIn;
    }
    
    return 0;
}


extern "C" double pyDiff_pYavg(double *params, double *dIn) {
    double tps[2] = {0, 10};
    double GasIn[80];
    
    for (size_t ii = 0; ii < NELEMS(GasIn); ii++) {
        GasIn[ii] = 0;
    }
    
    GasIn[0] = 16;
    GasIn[1] = 16;
    
    double dataPtr[2];
    
    int retVal = matlabDiffTPS_pYavg(dataPtr, GasIn, 80, params, tps, 2, dIn, 1, 1, 0);
    
    if (retVal == 0) {
        return dataPtr[1] / dataPtr[0];
    } else {
        return -1;
    }
}




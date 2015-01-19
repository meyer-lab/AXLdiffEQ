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
    return calcError(Param(pIn));
}

extern "C" void pyEntryVec(double *pIn, double *pOut, int n) {
    for (int ii = 0; ii < n; ii++) {
        if (pIn[ii*11 + 0] < pIn[ii*11 + 2]) {
            pOut[ii] = 1E6;
            continue;
        }
        
        pOut[ii] = calcError(Param(&pIn[ii*11]));
    }
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



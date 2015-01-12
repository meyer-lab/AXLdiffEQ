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

extern "C" double pyEntryNew(double *pIn) {
    return calcError(ParamNew(pIn));
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



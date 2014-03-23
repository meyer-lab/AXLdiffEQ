/**
 *   \file BlasHeader.cpp
 *   FiniteDiffGSL
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

#include <queue>
#include <iostream>
#include <thread>
#include "CVodeHelpers.h"
#include "ReactionCode.h"
#include "ModelRunning.h"
#include "BlasHeader.h"

using namespace std;

/// Queue for multithreaded evaluation of parameter sets
struct queT {
    param_type In; ///< Parameter set to be evaluated
    double *out;   ///< Evaluated parameter set performance
};


/// This is the entry for MatLab calling to get the fitness of a particular parameter set 
/// by running the simulations and matching it to the data. This function runs in parallel
/// with queuing because the optimization process is so computationally intensive. This
/// requires vectorized calling on the MatLab side.
///
/// Input:
/// 	dataPtr: pointer to double array of length nIn for output fitness
///		pIn: pointer to sets of parameters
/// 	nIn: number of parameter sets inputted
/// Output: 0 regardless. Fitness will be 1e10 if simulation fails.
extern "C" int matlabEntry(double *dataPtr, double *pIn, int nIn) {
    const int nThreads = 8;
    queue<queT> runThese;
    thread t[nThreads];
    atomic<bool> done[nThreads];
    queT in;
    
    for (size_t ii = 0; ii < (size_t) abs(nIn); ii++) {
        param_type pInSlice;
        
        for (size_t jj = 0; jj < pInSlice.size(); jj++) {
            pInSlice[jj] = pIn[(size_t) ii*pInSlice.size() + jj];
        }
        
        in.In = pInSlice;
        in.out = &dataPtr[ii];
        runThese.push(in);
     }

    
    for (int ii = 0; ii < nThreads; ii++) {
        
        // End if the queue is empty
        if (runThese.size() == 0) break;
        
        t[ii] = std::thread(calcErrorRef,runThese.front().In,runThese.front().out, &done[ii]);
        runThese.pop();
    }
    
    while (runThese.size() > 0) {
        for (int ii = 0; ii < nThreads; ii++) {
            
            // End if the queue is empty
            if (runThese.size() == 0) break;
            
            
            if (done[ii] == true) {
                t[ii].join();
                done[ii] = false;
                t[ii] = std::thread(calcErrorRef,runThese.front().In,runThese.front().out, &done[ii]);
                runThese.pop();
            }
        }
    }
    
    // Clear out the last running threads
    for (int ii = 0; ii < std::min(nThreads,nIn); ii++) t[ii].join();
    
    return 0;
}

extern "C" double calcErrorOneCellLine (int cellLine, const double *pIn) {
    const size_t numPs = 17;
    
    param_type pInTemp;
    
    for (size_t ii = 0; ii < (numPs - 2); ii++) {
        pInTemp[ii] = pIn[ii];
        //cout << pIn[ii] << endl;
    }
    
    struct rates pInC = Param(pInTemp);
    pInC.expression = pIn[16];
    
    return calcErrorOneLine (pInC, (size_t) cellLine, pIn[15]);
}

/// Function for fitting a parameter set from R
extern "C" void rEntry(double *dataPtr, const double *pIn) {
    param_type pInSlice;
    for (size_t jj = 0; jj < pInSlice.size(); jj++) {
        pInSlice[jj] = pIn[jj];
        cout << pIn[jj] << endl;
    }
    
    *dataPtr = calcError(pInSlice);
}

extern "C" int calcProfileMatlab(double *dataPtr, double *params, double *tps, int nTps, double autocrine, double AXL, double GasStim, int frac) {
    param_type pIn;
    
    for (size_t ii = 0; ii < pIn.size(); ii++) pIn[ii] = params[ii];
    
    try {
        calcProfileSet (dataPtr, tps, Param(pIn), nTps, autocrine, AXL, GasStim, frac);
    } catch (exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn) {
    
    try {
        diffusionSolution(dataPtr, AXLin, GasIn, gridIn, autocrine, params, tps, nTps, dIn, endoImpairIn, degImpairIn);
    } catch (exception &e) {
        cout << "Failed twice." << endl;
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS_pY(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps*Nspecies];
    
    int flag = matlabDiffTPS(dataPtrTemp, AXLin, GasIn, gridIn, autocrine, params, tps, nTps, dIn, endoImpairIn, degImpairIn);
    if (flag == 1) return(1);
    
    // Create the parameter structure
    param_type pIn;
    struct rates pInS;
    for (size_t ii = 0; ii < pIn.size(); ii++) {
        pIn[ii] = params[ii];
    }
    pInS = Param(pIn);
    

    N_Vector state = N_VNew_Serial(Nspecies);
    
    for (size_t time = 0; time < (size_t) abs(nTps); time++) {
        for (size_t gridP = 0; gridP < (size_t) abs(gridIn); gridP++) {
            for (size_t spec = 0; spec < Nspecies; spec++) {
                Ith(state,spec) = dataPtrTemp[time*(((size_t) abs(gridIn))*Nspecies) + spec*((size_t) abs(gridIn)) + gridP];
            }
            
            if (frac == 0) {
                dataPtr[time*((size_t) abs(gridIn)) + gridP] = pYcalc(state, pInS);
            } else if (frac == 1) {
                dataPtr[time*((size_t) abs(gridIn)) + gridP] = pYcalc(state, pInS) / totCalc(state);
            } else {
                dataPtr[time*((size_t) abs(gridIn)) + gridP] = totCalc(state);
            }
        }
    }
    
    N_VDestroy_Serial(state);
    return 0;
}

extern "C" int matlabDiffTPS_pYavg(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps];
    
    int flag = matlabDiffTPS_pY(dataPtrTemp, AXLin, GasIn, gridIn, autocrine, params, tps, nTps, dIn, endoImpairIn, degImpairIn, frac);
    if (flag == 1) return(1);
    
    double summ;
    
    for (size_t time = 0; time < (size_t) abs(nTps); time++) {
        summ = 0;
        
        for (size_t gridP = 0; gridP < (size_t) abs(gridIn); gridP++) {
            summ += dataPtrTemp[time*((size_t) abs(gridIn)) + gridP]/(gridIn)*(gridP);
        }
        
        dataPtr[time] = 2*summ/gridIn;
    }
    
    return 0;
}


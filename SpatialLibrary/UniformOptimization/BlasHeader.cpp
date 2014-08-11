/**
 *   \file BlasHeader.cpp
 *   FiniteDiffGSL
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

#include <queue>
#include <iostream>
#include <thread>

#ifdef __clang__
#include <atomic>
#else
#include <stdatomic.h>
#endif

#include <sstream>
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
    atomic<bool> *done = new atomic<bool>(nThreads);
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
        t[ii] = thread(calcErrorRef,runThese.front().In,runThese.front().out, &done[ii]);
        runThese.pop();
    }
    
    while (runThese.size() > 0) {
        for (int ii = 0; ii < nThreads; ii++) {
            
            // End if the queue is empty
            if (runThese.size() == 0) break;
            
            
            if (done[ii] == true) {
                t[ii].join();
                done[ii] = false;
                t[ii] = thread(calcErrorRef,runThese.front().In,runThese.front().out, &done[ii]);
                runThese.pop();
            }
        }
    }
    
    // Clear out the last running threads
    for (int ii = 0; ii < std::min(nThreads,nIn); ii++) t[ii].join();
    
    delete done;
    return 0;
}

extern "C" int matlabEntryA549(double *dataPtr, double *pIn, int nIn) {
    const int nThreads = 8;
    queue<queT> runThese;
    thread t[nThreads];
    atomic<bool> *done = new atomic<bool>(nThreads);
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
        t[ii] = thread(calcErrorRefA549,runThese.front().In,runThese.front().out, &done[ii]);
        runThese.pop();
    }
    
    while (runThese.size() > 0) {
        for (int ii = 0; ii < nThreads; ii++) {
            
            // End if the queue is empty
            if (runThese.size() == 0) break;
            
            
            if (done[ii] == true) {
                t[ii].join();
                done[ii] = false;
                t[ii] = thread(calcErrorRefA549,runThese.front().In,runThese.front().out, &done[ii]);
                runThese.pop();
            }
        }
    }
    
    // Clear out the last running threads
    for (int ii = 0; ii < std::min(nThreads,nIn); ii++) t[ii].join();
    
    delete done;
    return 0;
}


extern "C" double multiPyEntry(double *pIn) {
    queue<queT> runThese;
    thread t[4];
    
    struct rates params[4];
    double data[4];
    
    for (size_t ii = 0; ii < 4; ii++) params[ii] = Param_multi(pIn);
    params[0].expression = pIn[16];
    params[1].expression = pIn[17];
    params[2].expression = pIn[18];
    params[3].expression = pIn[19];
    params[0].scaleA = pIn[20];
    params[1].scaleA = pIn[21];
    params[2].scaleA = pIn[22];
    params[3].scaleA = pIn[23];

    t[0] = thread(A549Multi,        params[0],    pIn[24], &data[0]); // inP, cellLine, autocrine, data
    t[1] = thread(oneCellLineMulti, params[1], 0, pIn[25], &data[1]); // inP, cellLine, autocrine, data
    t[2] = thread(oneCellLineMulti, params[2], 2, pIn[26], &data[2]); // inP, cellLine, autocrine, data
    t[3] = thread(oneCellLineMulti, params[3], 3, pIn[27], &data[3]); // inP, cellLine, autocrine, data
    
    
    // Clear out the last running threads
    for (int ii = 0; ii < 4; ii++) t[ii].join();

    return data[0] + data[1] + data[2] + data[3];
}


extern "C" int matlabEntryA549VaryEndo(double *dataPtr, double *pIn, int nIn) {
    const int nThreads = 8;
    queue<queT> runThese;
    thread t[nThreads];
    atomic<bool> *done = new atomic<bool>(nThreads);
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
        t[ii] = thread(calcErrorRefA549,runThese.front().In,runThese.front().out, &done[ii]);
        runThese.pop();
    }
    
    while (runThese.size() > 0) {
        for (int ii = 0; ii < nThreads; ii++) {
            
            // End if the queue is empty
            if (runThese.size() == 0) break;
            
            
            if (done[ii] == true) {
                t[ii].join();
                done[ii] = false;
                t[ii] = thread(calcErrorRefA549VaryEndo,runThese.front().In,runThese.front().out, &done[ii]);
                runThese.pop();
            }
        }
    }
    
    // Clear out the last running threads
    for (int ii = 0; ii < std::min(nThreads,nIn); ii++) t[ii].join();
    
    delete done;
    return 0;
}

extern "C" double matlabEntryA549VaryEndoPy(double *pIn) {
    atomic<bool> done;
    param_type pInSlice;
    
    for (size_t jj = 0; jj < pInSlice.size(); jj++) {
        pInSlice[jj] = pIn[jj];
    }
    double error = -1;
    
    calcErrorRefA549VaryEndo(pInSlice, &error, &done);

    return error;
}

extern "C" double matlabEntryA549VaryEndoPyRed(double *pIn) {
    atomic<bool> done;
    param_type pInSlice;
    double error = -1;
    
    pInSlice[0] = 1.2;
    pInSlice[1] = pIn[0];
    pInSlice[2] = 0.042;
    
    for (size_t i = 3; i < 14; i++) pInSlice[i] = pIn[i - 2];
    
    pInSlice[14] = 1; // fD2
    pInSlice[15] = pIn[12];
    pInSlice[16] = pIn[13];
    pInSlice[17] = pIn[14];
    pInSlice[18] = pIn[15];
    
    calcErrorRefA549VaryEndo(pInSlice, &error, &done);
    
    return error;
}

extern "C" double matlabEntryBT549VaryEndoPy(double *pIn) {
    param_type pInSlice;
    
    for (size_t j = 0; j < pInSlice.size(); j++) pInSlice[j] = pIn[j];
    
    struct rates params = Param(pInSlice);
    
    params.expression = pInSlice[16];
    params.internalFrac = pInSlice[17];
    params.internalV = pInSlice[18];
    
    return calcErrorOneLine (params, 3, pInSlice[15]);
}

extern "C" double matlabEntryBT549VaryEndoPyRed(double *pIn) {
    param_type pInS;
    
    pInS[0] = 1.2;
    pInS[1] = pIn[0];
    pInS[2] = 0.042;
    pInS[3] = pIn[1];
    pInS[4] = pIn[2];
    pInS[5] = pIn[3];
    pInS[6] = pIn[4];
    pInS[7] = pIn[5];
    pInS[8] = pIn[6];
    pInS[9] = pIn[7];
    pInS[10] = pIn[8];
    pInS[11] = pIn[9];
    pInS[12] = pIn[10];
    pInS[13] = pIn[11];
    pInS[14] = 1; // fD2
    
    struct rates params = Param(pInS);
    
    params.expression = pIn[13];
    params.internalFrac = pIn[14];
    params.internalV = pIn[15];
    
    return calcErrorOneLine (params, 3, pIn[12]);
}


extern "C" int calcProfileMatlab(double *dataPtr, double *params, double *tps, int nTps, double autocrine, double AXL, double GasStim, int frac) {
    param_type pIn;
    struct rates pInS = Param(pIn);
    
    for (size_t ii = 0; ii < pIn.size(); ii++) pIn[ii] = params[ii];
    
    try {
        calcProfileSet (dataPtr, tps, &pInS, nTps, autocrine, AXL, GasStim, frac);
    } catch (std::exception &e) {
        errorLogger(&e);
        return 1;
    }
    
    return 0;
}

extern "C" int matlabDiffTPS(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn) {
    
    try {
        diffusionSolution(dataPtr, AXLin, GasIn, gridIn, autocrine, params, tps, nTps, dIn, endoImpairIn, degImpairIn);
    } catch (std::exception &e) {
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


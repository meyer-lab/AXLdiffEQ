//
//  ModelRunning.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include "sundials/sundials_nvector.h"
#include "cvode/cvode.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <exception>
#include <sstream>
#include <fstream>
#include <chrono>
#include <cobyla.h>
#include <cmath>
#include "ModelRunning.h"
#include "CVodeHelpers.h"

using namespace std;

#define NVp N_VGetArrayPointer_Serial

struct inData {
    N_Vector fitt;           ///< Calculated model values.
    const double *pYmeas;    ///< pY measurement.
    const double *errorMeas; ///< Error for pY measurement.
};

void calcParallel (double *outData, double *totData, double *surfData, size_t stimuli, N_Vector init_state, int *retVal, struct rates *params) {
    *retVal = 0;
    N_Vector state = N_VClone(init_state);
    void *cvode_mem = solver_setup(state, params, AXL_react);
    // Initialize state based on autocrine ligand
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(state);
        CVodeFree(&cvode_mem);
        *retVal = 1;
        return;
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    for (int xx = 0; xx < Nspecies; xx++) {
        Ith(state,xx) = Ith(init_state,xx);
    }
    
    Ith(state,0) += Gass[stimuli];
    double t = 0;
    
    solverReset(cvode_mem, state);
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    
    for (unsigned int ii = 0; ii < NELEMS(times); ii++) {
        int flag = CVode(cvode_mem, times[ii], state, &t, CV_NORMAL);
        
        if (flag < 0) {
            N_VDestroy_Serial(state);
            CVodeFree(&cvode_mem);
            *retVal = 1;
            return;
        }
        
        outData[stimuli*NELEMS(times) + ii] = pYcalc(state,params);
        totData[stimuli*NELEMS(times) + ii] = totCalc(state,params);
        surfData[stimuli*NELEMS(times) + ii] = surfAXL(state);
    }
    
    N_VDestroy_Serial(state);
    CVodeFree(&cvode_mem);
}


void calcParallelFull (double *outData, double *totData, double *surfData, size_t stimuli, N_Vector init_state, int *retVal, struct rates *params) {
    *retVal = 0;
    N_Vector state = N_VClone(init_state);
    void *cvode_mem = solver_setup(state, params, AXL_react);
    // Initialize state based on autocrine ligand
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(state);
        CVodeFree(&cvode_mem);
        *retVal = 1;
        return;
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    for (int xx = 0; xx < Nspecies; xx++) {
        Ith(state,xx) = Ith(init_state,xx);
    }
    
    Ith(state,0) += GassDoseFull[stimuli];
    double t = 0;
    
    solverReset(cvode_mem, state);
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    
    for (unsigned int ii = 0; ii < NELEMS(timesFull); ii++) {
        int flag = CVode(cvode_mem, timesFull[ii], state, &t, CV_NORMAL);
        
        if (flag < 0) {
            N_VDestroy_Serial(state);
            CVodeFree(&cvode_mem);
            *retVal = 1;
            return;
        }
        
        outData[stimuli*NELEMS(timesFull) + ii] = pYcalc(state,params);
        totData[stimuli*NELEMS(timesFull) + ii] = totCalc(state,params);
        surfData[stimuli*NELEMS(timesFull) + ii] = surfAXL(state);
    }
    
    N_VDestroy_Serial(state);
    CVodeFree(&cvode_mem);
}




double pyTotOnSi (struct rates *pIn, N_Vector init_state, int *retVal) {
    *retVal = 0;
    struct rates params = *pIn;
    N_Vector state = N_VNew_Serial(Nspecies);
    params.autocrine = 0;
    
    void *cvode_mem = initState(state, &params);
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(state);
        *retVal = 1;
        return 0;
    }
    
    double siTot = totCalc(state, &params);
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(state);
    
    return siTot / totCalc(init_state, &params);
}


void calcKinetic (double *outData, double *totData, double *surfData, double *siRatio, struct rates *params) {
    N_Vector init_state = N_VNew_Serial(Nspecies);
    int retVal[7] = {0, 0, 0, 0, 0, 0, 0};
    thread t[6];
    
    void *cvode_mem = initState(init_state, params);
    // Initialize state based on autocrine ligand
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(init_state);
        throw runtime_error(string("Error during solver threads."));
        return;
    }
    
    for (size_t ii = 0; ii < 6; ii++) t[ii] = std::thread(calcParallel, outData, totData, surfData, ii, init_state, &retVal[ii], params);
    
    *siRatio = pyTotOnSi(params, init_state, &retVal[6]);
    
    CVodeFree(&cvode_mem);
    
    for (size_t ii = 0; ii < NELEMS(t); ii++) t[ii].join();
    
    for (size_t ii = 0; ii < NELEMS(retVal); ii++) {
        if (retVal[ii] == 1) {
            N_VDestroy_Serial(init_state);
            throw runtime_error(string("Error during solver threads."));
            return;
        }
    }
    
    N_VDestroy_Serial(init_state);
}



void calcKineticFull (double *outData, double *totData, double *surfData, double *, struct rates *params) {
    N_Vector init_state = N_VNew_Serial(Nspecies);
    int retVal[7] = {0, 0, 0, 0, 0, 0, 0};
    thread t[7];
    
    void *cvode_mem = initState(init_state, params);
    // Initialize state based on autocrine ligand
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(init_state);
        throw runtime_error(string("Error during solver threads."));
        return;
    }
    
    for (size_t ii = 0; ii < 7; ii++) t[ii] = std::thread(calcParallelFull, outData, totData, surfData, ii, init_state, &retVal[ii], params);
    
    //*siRatio = pyTotOnSi(params, init_state, &retVal[6]);
    
    CVodeFree(&cvode_mem);
    
    for (size_t ii = 0; ii < NELEMS(t); ii++) t[ii].join();
    
    for (size_t ii = 0; ii < NELEMS(retVal); ii++) {
        if (retVal[ii] == 1) {
            N_VDestroy_Serial(init_state);
            throw runtime_error(string("Error during solver threads."));
            return;
        }
    }
    
    N_VDestroy_Serial(init_state);
}







double errorFunc (double fitt, double pYmeas, double errorMeas) {
    return pow((((double)fitt) - pYmeas) / errorMeas, 2) / 2;
}

double errorOpt(unsigned, const double *x, double *, void *data) {
    struct inData *dataS = (struct inData *) data;
    double xx = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(dataS->fitt); ii++)
        errorFunc((double) Ith(dataS->fitt,ii) * x[0], dataS->pYmeas[ii], dataS->errorMeas[ii]);
    
    return xx;
}

double initialCondition (struct inData *dataS) {
    double meas = 0;
    double fit = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(dataS->fitt); ii++) {
        meas += Ith(dataS->fitt,ii);
        fit += dataS->pYmeas[ii];
    }
    
    return fit / meas;
}

double errorFuncOpt (N_Vector fitt, const double *pYmeas, const double *errorMeas) {
    struct inData dataS;
    dataS.fitt = fitt;
    dataS.pYmeas = pYmeas;
    dataS.errorMeas = errorMeas;
    
    double ff = 0;
    double xx = initialCondition(&dataS);
    const double lower = xx/2;
    const double upper = xx*2;
    double dx = 3*xx/8;
    double del = 1E-8;
    
    nlopt_stopping stop;
    stop.n = 0;
    stop.minf_max = 0.0;
    stop.ftol_rel = 0;
    stop.ftol_abs = 0;
    stop.xtol_rel = 1E-8;
    stop.xtol_abs = &del;
    stop.nevals = 0;
    stop.maxeval = 1E9;
    stop.force_stop = 0;
    
    int flag = cobyla_minimize(1, errorOpt, &dataS, 0, NULL, 0, NULL, &lower, &upper, &xx, &ff, &stop, &dx);
    
    if (flag < 0) throw runtime_error(string("Error during error optimization step."));
    
    return ff;
}


double errorFuncFix (N_Vector fitt, const double *pYmeas, const double *errorMeas) {
    double xx = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(fitt); ii++)
        xx += errorFunc((double) Ith(fitt,ii), pYmeas[ii], errorMeas[ii]);
    
    return xx;
}

double calcError (struct rates inP) {
    N_Vector outData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector totData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector surfData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    double siRatio = 0;
    
    double error = 0;
    
    try {
        calcKinetic(NVp(outData), NVp(totData), NVp(surfData), &siRatio, &inP);
        
        error += errorFuncOpt (outData, pY[0], pYerror[0]);
        error += errorFuncFix (totData, tot[0], totError[0]);
        error += errorFuncOpt (surfData, surf[0], surfError[0]);
        
        error += errorFunc(siRatio, 1.487, 0.0571);
    } catch (runtime_error &e) {
        errorLogger(&e);
        error = 1E8;
    }
    
    N_VDestroy_Serial(outData);
    N_VDestroy_Serial(totData);
    N_VDestroy_Serial(surfData);
    
    return error;
}


double calcErrorFull (struct rates inP) {
    N_Vector outData = N_VNew_Serial(NELEMS(GassDoseFull)*NELEMS(timesFull));
    N_Vector totData = N_VNew_Serial(NELEMS(GassDoseFull)*NELEMS(timesFull));
    N_Vector surfData = N_VNew_Serial(NELEMS(GassDoseFull)*NELEMS(timesFull));
    double siRatio = 0;
    
    double error = 0;
    
    try {
        calcKineticFull(NVp(outData), NVp(totData), NVp(surfData), &siRatio, &inP);
        
        error += errorFuncOpt (outData, pYdoseFull, pYdoseFullError);
        error += errorFuncFix (totData, DoseTotFull, DoseTotFullErr);
        error += errorFuncOpt (surfData, surfDoseFull, surfDoseFullError);
        
        //error += errorFunc(siRatio, 1.487, 0.0571);
    } catch (runtime_error &e) {
        errorLogger(&e);
        error = 1E8;
    }
    
    N_VDestroy_Serial(outData);
    N_VDestroy_Serial(totData);
    N_VDestroy_Serial(surfData);
    
    return error;
}




// Calculate the initial state by waiting a long time with autocrine Gas
void *initState( N_Vector init, struct rates *params) {
    endoImpair = 1.0;
    degImpair = 1.0;
    double t;
    
    for (int ii = 0; ii < Nspecies ; ii++) Ith(init,ii) = 0;
    
    Ith(init,0) = params->autocrine;
    
    void *cvode_mem = solver_setup (init, params, AXL_react);
    if (cvode_mem == NULL) return NULL;
    
    int flag = CVode(cvode_mem, autocrineT, init, &t, CV_NORMAL);
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        return NULL;
    }
    
    /* Free integrator memory */
    return cvode_mem;
}



/// Calculate phosphorylation at time points measured
void calcProfileSet (double *outData, double *tps, struct rates *params, int nTps, double GasStim, int frac) {
    N_Vector state = N_VNew_Serial(Nspecies);
    
    double t; ///< Time position of the solver.
    
    void *cvode_mem = NULL;
    int flag;
    
    // Initialize state based on autocrine ligand
    try {
        cvode_mem = initState(state, params);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    Ith(state,0) += GasStim;
    t = 0;
    
    try {
        solverReset(cvode_mem, state);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        CVodeFree(&cvode_mem);
        throw;
    }
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    size_t ii = 0;
    
    if (tps[0] == 0) {
        if (frac == 0) {
            outData[ii] = pYcalc(state,params);
        } else if (frac == 1) {
            outData[ii] = pYcalc(state,params) / totCalc(state,params);
        } else if (frac == 2) {
            outData[ii] = totCalc(state,params);
        } else if (frac == 3) {
            outData[ii] = surfpY(state,params);
        } else if (frac == 4) {
            outData[ii] = surfAXL(state);
        } else {
            throw runtime_error(string("Bad option."));
        }
        
        ii = 1;
    }
    
    for (; ii < (size_t) abs(nTps); ii++) {
        flag = CVode(cvode_mem, tps[ii], state, &t, CV_NORMAL);
        if (flag < 0) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(state);
            throw runtime_error(string("Error at CVode Time Course."));
        }
        
        if (frac == 0) {
            outData[ii] = pYcalc(state,params);
        } else if (frac == 1) {
            outData[ii] = pYcalc(state,params) / totCalc(state,params);
        } else if (frac == 2) {
            outData[ii] = totCalc(state,params);
        } else if (frac == 3) {
            outData[ii] = surfpY(state,params);
        } else if (frac == 4) {
            outData[ii] = surfAXL(state);
        } else {
            throw runtime_error(string("Bad option."));
        }
    }
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(state);
}

void errorLogger (exception *e) {
    ofstream errOut;
    
    if (print_CV_err == 0) return;
    else if (print_CV_err == 1) {
        cout << e->what() << endl;
    } else if (print_CV_err == 2) {
        errOut.open ("error.txt", ios::app);
        errOut << e->what() << endl;
        errOut.close();
    }
}

void errorLogger (stringstream &e) {
    ofstream errOut;
    
    if (print_CV_err == 0) return;
    else if (print_CV_err == 1) {
        cout << e.str() << endl;
    } else if (print_CV_err == 2) {
        errOut.open ("error.txt", ios::app);
        errOut << e.str() << endl;
        errOut.close();
    }
}



void diffusionSolution(double *dataPtr, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn) {
    
    // Common
    double t = 0;
    void *cvode_mem = NULL;
    int flag;
    
    for (int ii = 0; ii < Nspecies; ii++) diffD[ii] = dIn[ii];
    
    // Create the parameter structure
    struct rates pInS = Param(params);
    
    // Get initial state
    N_Vector init_state = N_VNew_Serial(Nspecies);
    
    // Initialize state based on autocrine ligand
    try {
        cvode_mem = initState(init_state, &pInS);
    } catch (exception &e) {
        N_VDestroy_Serial(init_state);
        throw;
    }
    
    CVodeFree(&cvode_mem);
    
    endoImpair = endoImpairIn;
    degImpair = degImpairIn;
    
    
    struct diffRates pInD;
    pInD.params = pInS;
    pInD.reactIn = N_VNew_Serial(Nspecies);
    pInD.reactOut = N_VNew_Serial(Nspecies);
    
    
    // Initialize full diffusion model
    N_Vector state = N_VNew_Serial(Nspecies * gridIn);
    
    for (size_t ii = 0; ii < (size_t) abs(gridIn); ii++) {
        for (size_t spec = 0; spec < Nspecies; spec++) Ith(state,spec*((size_t) gridIn) + ii) = Ith(init_state,spec);
    }
    
    N_VDestroy_Serial(init_state);
    
    for (size_t ii = 0; ii < (size_t) abs(gridIn); ii++) Ith(state,ii) = GasIn[ii];
    // Done initializing diffusion model
    
    try {
        cvode_mem = solver_setup (state, &pInD, AXL_react_diff);
        
        if (cvode_mem == NULL) throw runtime_error(string("Uncaught NULL return from solver_setup in diffusion setup."));
    } catch (exception &e) {
        N_VDestroy_Serial(pInD.reactIn);
        N_VDestroy_Serial(pInD.reactOut);
        N_VDestroy_Serial(state);
        throw;
    }
    
    size_t tIDX = 0;
    
    if (tps[0] == 0) {
        for (size_t jj = 0; jj < (size_t) NV_LENGTH_S(state); jj++) dataPtr[jj] = Ith(state,jj);
        tIDX = 1;
    }
    
    for (; tIDX < (size_t) abs(nTps); tIDX++) {
        flag = CVode(cvode_mem, tps[tIDX], state, &t, CV_NORMAL);
        
        if (flag < 0) {
            N_VDestroy_Serial(pInD.reactIn);
            N_VDestroy_Serial(pInD.reactOut);
            N_VDestroy_Serial(state);
            CVodeFree(&cvode_mem);
            throw runtime_error(string("CVode error on diffusion solution."));
        }
        
        for (size_t jj = 0; jj < (size_t) NV_LENGTH_S(state); jj++) {
            dataPtr[tIDX*((size_t) NV_LENGTH_S(state)) + jj] = Ith(state,jj);
        }
    }
    
    /* Free vectors */
    N_VDestroy_Serial(pInD.reactIn);
    N_VDestroy_Serial(pInD.reactOut);
    N_VDestroy_Serial(state);
    CVodeFree(&cvode_mem);
}

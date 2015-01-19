//
//  ModelRunning.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include "sundials_nvector.h"
#include "cvode.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <exception>
#include "cobyla.h"
#include <cmath>
#include "cvode_impl.h"
#include "ModelRunning.h"
#include "CVodeHelpers.h"

using namespace std;

const double fgMgConv = 135.2;

static int AXL_react(double, double *x_d, double *dxdt_d, void *user_data) {
    struct rates *r = (struct rates *) user_data;
    
    // 0 AXL   // 1 A1    // 2 A2
    // 3 A12    // 4 D1    // 5 D2    // 6 AXLi
    // 7 A1i    // 8 A2i   // 9 A12i // 10 D1i    // 11 D2i   // 12 Gasi
    
    const double dR1 = r->Binding1 * x_d[0] * r->gasCur - r->Unbinding1 * x_d[1];
    const double dR2 = r->Binding2 * x_d[0] * r->gasCur - r->Unbinding2 * x_d[2];
    const double dR3 = r->Binding2 * x_d[1] * r->gasCur - r->Unbinding2 * x_d[3];
    const double dR4 = r->Binding1 * x_d[2] * r->gasCur - r->Unbinding1 * x_d[3];
    const double dR5 = r->xFwd1 * x_d[0] * x_d[1] - r->xRev1 * x_d[4];
    const double dR6 = r->xFwd2 * x_d[0] * x_d[2] - r->xRev2 * x_d[4];
    const double dR7 = r->xFwd3 * x_d[0] * x_d[3] - r->xRev3 * x_d[5];
    const double dR8 = r->xFwd4 * x_d[1] * x_d[1] - r->xRev4 * x_d[5];
    const double dR9 = r->xFwd5 * x_d[2] * x_d[2] - r->xRev5 * x_d[5];
    const double dR11 = r->xFwd6 * r->gasCur * x_d[4] - r->xRev6 * x_d[5];
    
    const double dR32 = r->Binding1 * x_d[6] * x_d[12] / r->internalV - r->Unbinding1 * x_d[7];
    const double dR33 = r->Binding2 * x_d[6] * x_d[12] / r->internalV - r->Unbinding2 * x_d[8];
    const double dR34 = r->Binding2 * x_d[7] * x_d[12] / r->internalV - r->Unbinding2 * x_d[9];
    const double dR35 = r->Binding1 * x_d[8] * x_d[12] / r->internalV - r->Unbinding1 * x_d[9];
    const double dR36 = r->xFwd1 * x_d[6] * x_d[7] - r->xRev1 * x_d[10];
    const double dR37 = r->xFwd2 * x_d[6] * x_d[8] - r->xRev2 * x_d[10];
    const double dR38 = r->xFwd3 * x_d[6] * x_d[9] - r->xRev3 * x_d[11];
    const double dR39 = r->xFwd4 * x_d[7] * x_d[7] - r->xRev4 * x_d[11]; // Checked
    const double dR40 = r->xFwd5 * x_d[8] * x_d[8] - r->xRev5 * x_d[11]; // Checked
    const double dR41 = r->xFwd6 * x_d[12] * x_d[10] / r->internalV - r->xRev6 * x_d[11]; // Checked
    
    dxdt_d[0] = - dR7 - dR6 - dR5 - dR1 - dR2 + r->expression; // AXL
    dxdt_d[1] = -2*(dR8) - dR5 + dR1 - dR3                   ; // AXLgas1
    dxdt_d[2] = -2*(dR9) - dR6 + dR2 - dR4                   ; // AXLgas2
    dxdt_d[3] = -dR7 + dR3 + dR4                             ; // AXLgas12
    dxdt_d[4] = -dR11 + dR6 + dR5                            ; // AXLdimer1
    dxdt_d[5] = dR11 + dR9 + dR8 + dR7                       ; // AXLdimer2
    
    dxdt_d[6]  = - dR38 - dR37 - dR36 - dR32 - dR33          ; // AXLi
    dxdt_d[7]  = -2*(dR39) - dR36 + dR32 - dR34              ; // AXLgas1i
    dxdt_d[8]  = -2*(dR40) - dR37 + dR33 - dR35              ; // AXLgas2i
    dxdt_d[9] = -dR38 + dR34 + dR35                         ; // AXLgas12i
    dxdt_d[10] = -dR41 + dR37 + dR36                          ; // AXLdimer1i
    dxdt_d[11] = dR41 + dR40 + dR39 + dR38                   ; // AXLdimer2i
    
    dxdt_d[12] = -dR41 - dR32 - dR33 - dR34 - dR35 - r->kDeg*x_d[12];
    
    
    dxdt_d[0] += -x_d[0]*(r->internalize) + r->kRec*(1-r->fElse)*x_d[6]*r->internalFrac; // Endocytosis, recycling
    dxdt_d[6] += x_d[0]*(r->internalize)/r->internalFrac - r->kRec*(1-r->fElse)*x_d[6] - r->kDeg*r->fElse*x_d[6]; // Endocytosis, recycling, degradation
    
    for (int ii = 1; ii < 4; ii++) {
        dxdt_d[ii]  += -x_d[ii]*r->internalize + r->kRec*(1-r->fElse)*x_d[ii+6]*r->internalFrac; // Endocytosis, recycling
        dxdt_d[ii+6] += x_d[ii]*r->internalize/r->internalFrac - r->kRec*(1-r->fElse)*x_d[ii+6] // Endocytosis, recycling
        - r->kDeg*r->fElse*x_d[ii+6]; // Degradation
    }
    
    // D1 trafficking
    if (r->pD1 == 1) {
        dxdt_d[4]  += -x_d[4]*(r->internalize + r->pYinternalize) + r->kRec*(1-r->fD2)*x_d[10]*r->internalFrac; // Endocytosis, recycling
        dxdt_d[10] += x_d[4]*(r->internalize + r->pYinternalize)/r->internalFrac - r->kRec*(1-r->fD2)*x_d[10] - r->kDeg*r->fD2*x_d[10]; // Endocytosis, recycling, degradation
    } else {
        dxdt_d[4]  += -x_d[4]*r->internalize + r->kRec*(1-r->fElse)*x_d[10]*r->internalFrac; // Endocytosis, recycling
        dxdt_d[10] += x_d[4]*r->internalize/r->internalFrac - r->kRec*(1-r->fElse)*x_d[10] // Endocytosis, recycling
        - r->kDeg*r->fElse*x_d[10]; // Degradation
    }
    
    dxdt_d[5]  += -x_d[5]*(r->internalize + r->pYinternalize) + r->kRec*(1-r->fD2)*x_d[11]*r->internalFrac; // Endocytosis, recycling
    dxdt_d[11] += x_d[5]*(r->internalize + r->pYinternalize)/r->internalFrac - r->kRec*(1-r->fD2)*x_d[11] - r->kDeg*r->fD2*x_d[11]; // Endocytosis, recycling, degradation
    
    
    return 0;
}

int AXL_react(double t, N_Vector xIn, N_Vector dxdtIn, void *user_data) {
    double* x_d = NV_DATA_S(xIn);
    double* dxdt_d = NV_DATA_S(dxdtIn);
    
    return AXL_react(t, x_d, dxdt_d, user_data);
}

static double surfAXL (N_Vector state) {
    return Ith(state,0) + Ith(state,1) + Ith(state,2) + Ith(state,3) + 2*Ith(state,4) + 2*Ith(state,5);
}


// This takes the model state and calculates the amount of phosphorylated species
static double pYcalc (N_Vector state, struct rates *p) {
    if (p->pD1 == 1) {
        return 2*Ith(state,5) + 2*Ith(state,4) + p->internalFrac*(2*Ith(state,11) + 2*Ith(state,10));
    } else {
        return 2*Ith(state,5) + p->internalFrac*(2*Ith(state,11));
    }
}

// This takes the model state and calculates the amount of phosphorylated species
static double surfpY (N_Vector state, struct rates *p) {
    if (p->pD1 == 1) {
        return 2*Ith(state,5) + 2*Ith(state,4);
    } else {
        return 2*Ith(state,5);
    }
    
}

// This takes the model state and calculates the total amount of receptor in a cell
static double totCalc (N_Vector state, struct rates *p) {
    double total = 0;
    
    for (int ii = 0; ii < 6; ii++) total += Ith(state,ii);
    for (int ii = 6; ii < 12; ii++) total += Ith(state,ii)*p->internalFrac;
    
    total += Ith(state,4);
    total += Ith(state,5);
    total += Ith(state,10)*p->internalFrac;
    total += Ith(state,11)*p->internalFrac;
    
    return total/fgMgConv;
}

struct rates Param(double *params) {
    struct rates out;
    
    for (size_t ii = 0; ii < 11; ii++) {
        if (params[ii] < 0) throw invalid_argument(string("Input model parameter outside the physical range."));
    }
    
    out.Binding1 = 1.2;
    out.Binding2 = 0.06;
    out.Unbinding1 = 0.042;
    out.Unbinding2 = params[0];
    out.xFwd1 = params[1];
    out.xRev1 = params[0];
    out.xFwd3 = params[1];
    out.xRev3 = params[2];
    out.internalize = params[3];
    out.pYinternalize = params[4];
    out.kRec = params[5];
    out.kDeg = params[6];
    out.fElse = params[7];
    out.expression = params[8];
    out.autocrine = params[9];
    out.fD2 = 1;
    out.internalFrac = 0.5;
    out.internalV = 623;
    out.pD1 = (int) params[10];
    
    out.xRev5 = out.xRev3*out.Unbinding1/out.Unbinding2;
    out.xRev4 = out.xRev3*out.Unbinding2/out.Unbinding1;
    out.xRev2 = out.xRev1*out.Unbinding1/out.Unbinding2;
    out.xFwd2 = out.xFwd1*out.Binding1/out.Binding2;
    out.xFwd4 = out.xFwd3*out.Binding2/out.Binding1;
    out.xFwd5 = out.xFwd3*out.Binding1/out.Binding2;
    out.xFwd6 = out.xFwd3*out.Binding2/out.xFwd1;
    out.xRev6 = out.xRev3*out.Unbinding2/out.xRev1;
    
    return out;
}


/// END REACTION CODE


#define NVp N_VGetArrayPointer_Serial

struct inData {
    N_Vector fitt;           ///< Calculated model values.
    const double *pYmeas;    ///< pY measurement.
    const double *errorMeas; ///< Error for pY measurement.
};







static void calcParallelFull (double *outData, double *totData, double *surfData, size_t stimuli, N_Vector init_state, int *retVal, struct rates *params) {
    
    struct rates paramTwo = *params;
    
    *retVal = 0;
    N_Vector state = N_VClone(init_state);
    void *cvode_mem = solver_setup(state, &paramTwo, AXL_react);
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
    
    paramTwo.gasCur = paramTwo.gasCur + GassDoseFull[stimuli];
    CVodeSetUserData(cvode_mem, &paramTwo);
    
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

static void calcKinetic (double *outData, double *totData, double *surfData, double *earlyPY, struct rates *params) {
    N_Vector init_state = N_VNew_Serial(Nspecies);
    double t;
    
    void *cvode_mem = initState(init_state, params);
    // Initialize state based on autocrine ligand
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(init_state);
        throw runtime_error(string("Error during solver threads."));
        return;
    }
    
    //
    // This part will calculate the dose response
    //
    struct rates paramTwo = *params;
    N_Vector state = N_VClone(init_state);
    // Initialize state based on autocrine ligand
    
    for (size_t stimuli = 0; stimuli < 6; stimuli++) {
        for (int xx = 0; xx < Nspecies; xx++) Ith(state,xx) = Ith(init_state,xx);
        
        paramTwo.gasCur = paramTwo.gasCur + Gass[stimuli];
        CVodeSetUserData(cvode_mem, &paramTwo);
        
        t = 0;
        
        solverReset(cvode_mem, state);
        
        /* In loop, call CVode, print results, and test for error.
         Break out of loop when NOUT preset output times have been reached.  */
        
        for (unsigned int ii = 0; ii < NELEMS(times); ii++) {
            int flag = CVode(cvode_mem, times[ii], state, &t, CV_NORMAL);
            
            if (flag < 0) {
                N_VDestroy_Serial(state);
                CVodeFree(&cvode_mem);
                N_VDestroy_Serial(init_state);
                throw runtime_error(string("Error during solver threads."));
                return;
            }
            
            outData[stimuli*NELEMS(times) + ii] = pYcalc(state,params);
            totData[stimuli*NELEMS(times) + ii] = totCalc(state,params);
            surfData[stimuli*NELEMS(times) + ii] = surfAXL(state);
        }
        
        stimuli++;
    }
    
    //
    // This part calculates the kinetic response
    //
    earlyPY[0] = pYcalc(init_state,params);
    // Initialize state based on autocrine ligand
    
    /* We've got the initial state, so now run through the kinetic data */
    for (int xx = 0; xx < Nspecies; xx++) Ith(state,xx) = Ith(init_state,xx);
    
    paramTwo.gasCur = paramTwo.gasCur + 1.25;
    CVodeSetUserData(cvode_mem, &paramTwo);
    
    t = 0;
    
    solverReset(cvode_mem, state);
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    
    for (unsigned int ii = 1; ii < NELEMS(kTPS); ii++) {
        int flag = CVode(cvode_mem, kTPS[ii], state, &t, CV_NORMAL);
        
        if (flag < 0) {
            N_VDestroy_Serial(state);
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(init_state);
            throw runtime_error(string("Error during solver threads."));
            return;
        }
        
        earlyPY[ii] = pYcalc(state,params);
    }

    N_VDestroy_Serial(state);
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(init_state);
}



static void calcKineticFull (double *outData, double *totData, double *surfData, double *, struct rates *params) {
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

static double errorFunc (double fitt, double pYmeas, double errorMeas) {
    return pow((((double) fitt) - pYmeas) / errorMeas, 2) / 2;
}

static double errorOpt(unsigned, const double *x, double *, void *data) {
    struct inData *dataS = (struct inData *) data;
    double xx = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(dataS->fitt); ii++)
        xx += errorFunc((double) Ith(dataS->fitt,ii) * x[0], dataS->pYmeas[ii], dataS->errorMeas[ii]);
    
    return xx;
}

static double initialCondition (struct inData *dataS) {
    double meas = 0;
    double fit = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(dataS->fitt); ii++) {
        meas += Ith(dataS->fitt,ii);
        fit += dataS->pYmeas[ii];
    }
    
    return fit / meas;
}

static double errorFuncOpt (N_Vector fitt, const double *pYmeas, const double *errorMeas) {
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
    stop.ftol_rel = 1E-8;
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


static double errorFuncFix (N_Vector fitt, const double *pYmeas, const double *errorMeas) {
    double xx = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(fitt); ii++)
        xx += errorFunc((double) Ith(fitt,ii), pYmeas[ii], errorMeas[ii]);
    
    return xx;
}

double calcError (struct rates inP) {
    N_Vector outData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector totData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector surfData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector earlyPY = N_VNew_Serial(NELEMS(kTPS));
    
    double error = 0;
    
    try {
        calcKinetic(NVp(outData), NVp(totData), NVp(surfData), NVp(earlyPY), &inP);
        
        //error += errorFuncOpt (outData, pY[0], pYerror[0]);
        error += errorFuncFix (totData, tot[0], totError[0]);
        //error += errorFuncOpt (surfData, surf[0], surfError[0]);
        error += errorFuncOpt (earlyPY, pYk, pYkErr);
        
    } catch (runtime_error &e) {
        errorLogger(&e);
        error = 1E8;
    }
    
    N_VDestroy_Serial(outData);
    N_VDestroy_Serial(totData);
    N_VDestroy_Serial(surfData);
    N_VDestroy_Serial(earlyPY);
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
    double t;
    
    for (int ii = 0; ii < Nspecies ; ii++) Ith(init,ii) = 0;
    
    params->gasCur = params->autocrine;
    
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
    struct rates paramTwo = *params;
    
    double t; ///< Time position of the solver.
    
    void *cvode_mem = NULL;
    int flag;
    
    // Initialize state based on autocrine ligand
    try {
        cvode_mem = initState(state, &paramTwo);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    
    /* We've got the initial state, so now run through the kinetic data */
    paramTwo.gasCur = paramTwo.gasCur + GasStim;
    CVodeSetUserData(cvode_mem, &paramTwo);
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




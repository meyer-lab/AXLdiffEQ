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



double surfCalc (N_Vector state) {
    return (Ith(state,0) + Ith(state,1) + Ith(state,2) + Ith(state,3) + 2*Ith(state,4) + 2*Ith(state,5))/fgMgConv;
}


// This takes the model state and calculates the amount of phosphorylated species
double pYcalc (N_Vector state, struct rates *p) {
    if (p->pD1 == 1) {
        return 2*Ith(state,5) + 2*Ith(state,4) + p->internalFrac*(2*Ith(state,11) + 2*Ith(state,10));
    } else {
        return 2*Ith(state,5) + p->internalFrac*(2*Ith(state,11));
    }
}

// This takes the model state and calculates the total amount of receptor in a cell
double totCalc (const N_Vector state, const struct rates * const p) {
    double total = 0;
    
    for (int ii = 0; ii < 6; ii++) total += Ith(state,ii);
    for (int ii = 6; ii < 12; ii++) total += Ith(state,ii)*p->internalFrac;
    
    total += Ith(state,4);
    total += Ith(state,5);
    total += Ith(state,10)*p->internalFrac;
    total += Ith(state,11)*p->internalFrac;
    
    return total/fgMgConv;
}

struct rates Param(const double * const params) {
    struct rates out;
    
    if (min_element(params,params+(Nparams - 1)) < 0)
        throw invalid_argument(string("Parameter outside the physical range."));
        
        out.Binding1 = 1.2;
        out.Binding2 = 0.06;
        out.Unbinding1 = 0.042;
        out.Unbinding2 = params[0];
        
        out.xFwd1 = params[1];
        out.xFwd2 = out.xFwd1;
        out.xFwd3 = out.xFwd1;
        out.xFwd4 = out.xFwd1;
        out.xFwd5 = out.xFwd1;
        out.xFwd6 = out.Binding1;
        
        out.xRev4 = params[2];
        
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
        
        const double KD1 = out.Unbinding1 / out.Binding1;
        const double KD2 = out.Unbinding2 / out.Binding2;
        
        out.xRev1 = out.Unbinding2; /// Checked
        out.xRev2 = KD1*out.xRev1/KD2; // Checked
        out.xRev5 = out.xRev4*KD1*KD1/KD2/KD2; // Checked
        out.xRev6 = KD1*out.xFwd6*out.xRev4/out.xRev1; // Checked
        out.xRev3 = out.xRev4*KD1/KD2; // Checked
        
        
        
        return out;
}



/// END REACTION CODE


#define NVp N_VGetArrayPointer_Serial

struct inData {
    size_t N;
    const double *fitt;           ///< Calculated model values.
    const double *pYmeas;    ///< pY measurement.
    const double *errorMeas; ///< Error for pY measurement.
};

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
        
        paramTwo.gasCur = paramTwo.autocrine + Gass[stimuli];
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
            surfData[stimuli*NELEMS(times) + ii] = surfCalc(state);
        }
    }
    
    //
    // This part calculates the kinetic response
    //
    earlyPY[0] = pYcalc(init_state,params);
    // Initialize state based on autocrine ligand
    
    /* We've got the initial state, so now run through the kinetic data */
    for (int xx = 0; xx < Nspecies; xx++) Ith(state,xx) = Ith(init_state,xx);
    
    paramTwo.gasCur = paramTwo.autocrine + 1.25;
    
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

static double errorFunc (double fitt, double pYmeas, double errorMeas) {
    return pow((((double) fitt) - pYmeas) / errorMeas, 2) / 2;
}

static double errorOpt(unsigned, const double *x, double *, void *data) {
    struct inData *dataS = (struct inData *) data;
    double xx = 0;
    
    for (int ii = 0; ii < dataS->N; ii++)
        xx += errorFunc(dataS->fitt[ii] * x[0], dataS->pYmeas[ii], dataS->errorMeas[ii]);
    
    return xx;
}

static double initialCondition (struct inData *dataS) {
    double meas = 0;
    double fit = 0;
    
    for (int ii = 0; ii < dataS->N; ii++) {
        meas += dataS->fitt[ii];
        fit += dataS->pYmeas[ii];
    }
    
    return fit / meas;
}

static double errorFuncOpt (const double *fitt, const double *pYmeas, const double *errorMeas, size_t inN, double *xx) {
    struct inData dataS;
    dataS.fitt = fitt;
    dataS.pYmeas = pYmeas;
    dataS.errorMeas = errorMeas;
    dataS.N = inN;
    
    double ff = 0;
    xx[0] = initialCondition(&dataS);
    const double lower = xx[0]/2;
    const double upper = xx[0]*2;
    double dx = 3*xx[0]/8;
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
    
    int flag = cobyla_minimize(1, errorOpt, &dataS, 0, NULL, 0, NULL, &lower, &upper, xx, &ff, &stop, &dx);
    
    if (flag < 0) throw runtime_error(string("Error during error optimization step."));
    
    return ff;
}

static double errorFuncFix (const double *fitt, const double *pYmeas, const double *errorMeas, size_t inN) {
    double xx = 0;
    
    for (int ii = 0; ii < inN; ii++)
        xx += errorFunc(fitt[ii], pYmeas[ii], errorMeas[ii]);
    
    return xx;
}



double U87calcError (struct rates inP, double *fitParam) {
    double outData[NELEMS(Gass)*NELEMS(times)];
    double totData[NELEMS(Gass)*NELEMS(times)];
    double surfData[NELEMS(Gass)*NELEMS(times)];
    double earlyPY[NELEMS(kTPS)];
    
    double error = 0;
    
    try {
        calcKinetic(outData, totData, surfData, earlyPY, &inP);
        
        error += errorFuncOpt (outData, U87pY[0], U87pYerror[0], NELEMS(outData), &fitParam[0]);
        error += errorFuncOpt (totData, tot[0], totError[0], NELEMS(totData), &fitParam[1]);
        error += errorFuncOpt (earlyPY, U87pYk, U87pYkErr, NELEMS(earlyPY), &fitParam[2]);
        
    } catch (runtime_error &e) {
        errorLogger(&e);
        error = 1E8;
    }
    
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
void calcProfileSet (double *pYData, double *totData, double *surfData, double *speciesData, double *tps, struct rates *params, unsigned int nTps, double GasStim, double *convFac) {
    calcError(*params, convFac);
    
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
        pYData[ii] = pYcalc(state,params);
        totData[ii] = totCalc(state,params);
        surfData[ii] = surfCalc(state);
        
        for (size_t ss = 0; ss < Nspecies; ss++) {
            speciesData[ii*Nspecies + ss] = Ith(state,ss);
        }
        
        
        ii = 1;
    }
    
    for (; ii < nTps; ii++) {
        flag = CVode(cvode_mem, tps[ii], state, &t, CV_NORMAL);
        if (flag < 0) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(state);
            throw runtime_error(string("Error at CVode Time Course."));
        }
        
        pYData[ii] = pYcalc(state,params);
        totData[ii] = totCalc(state,params);
        surfData[ii] = surfCalc(state);
        
        for (size_t ss = 0; ss < Nspecies; ss++) {
            speciesData[ii*Nspecies + ss] = Ith(state,ss);
        }
    }
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(state);
}






static int AXL_react_diff(const double t, N_Vector xx , N_Vector dxxdt, void *user_data) {
    struct rates *pInD = (struct rates *) user_data;
    double reactIn[Nspecies];
    double reactOut[Nspecies];
    
    double* const xx_d = NV_DATA_S(xx);
    double* const dxxdt_d = NV_DATA_S(dxxdt);
    size_t pos, ii;
    const size_t grid_size = (size_t) NV_LENGTH_S(xx)/Nspecies;
    const double dRdRMaxRMaxR = maxR*maxR*(1.0/grid_size)*(1.0/grid_size);
    
    for (size_t spec = 0; spec < Nspecies; spec++) {
        if (pInD->diffD[spec] == 0) {
            for (ii = 1; ii < (grid_size-1); ii++) dxxdt_d[spec*grid_size + ii] = 0;
            dxxdt_d[spec*grid_size] = 0;
            dxxdt_d[(spec+1)*grid_size - 1] = 0;
        } else {
            for (ii = 1; ii < (grid_size-1); ii++) {
                pos = spec*grid_size + ii;
                dxxdt_d[pos] = (-4.0*xx_d[pos] + (2.0-1.0/ii)*xx_d[pos - 1] + (2.0+1.0/ii)*xx_d[pos + 1])/2/dRdRMaxRMaxR;
            }
            
            // Take care of boundary conditions
            dxxdt_d[spec*grid_size] = 4*(xx_d[spec*grid_size + 1] - xx_d[spec*grid_size])/dRdRMaxRMaxR;
            dxxdt_d[(spec+1)*grid_size - 1] = -4*(xx_d[(spec+1)*grid_size - 1] - xx_d[(spec+1)*grid_size - 2])/dRdRMaxRMaxR;
        }
    }
    
    // Add in the reaction for each location
    for (size_t jj = 0; jj < grid_size; jj++) {
        for (ii = 0; ii < Nspecies; ii++) reactIn[ii] = xx_d[ii*grid_size + jj];
        
        pInD->gasCur = pInD->gasProfile[jj];
        
        AXL_react(t,reactIn,reactOut,user_data);
        
        // Convert by diffusion coefficient and add in reaction
        for (ii = 0; ii < Nspecies; ii++) {
            dxxdt_d[ii*grid_size + jj] = dxxdt_d[ii*grid_size + jj]*pInD->diffD[ii] + reactOut[ii];
        }
    }
    
    return 0;
}



void diffusionSolution(double *dataPtr, double *GasIn, unsigned int gridIn, double *params, double *tps, unsigned int nTps, double *dIn) {
    // Common
    double t = 0;
    void *cvode_mem = NULL;
    int flag;
    
    // Create the parameter structure
    struct rates pInS = Param(params);
    pInS.diffD = dIn;
    pInS.gasProfile = GasIn;
    
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
    
    
    // Initialize full diffusion model
    N_Vector state = N_VNew_Serial(Nspecies * gridIn);
    
    for (size_t ii = 0; ii < gridIn; ii++) {
        for (size_t spec = 0; spec < Nspecies; spec++) Ith(state,spec*((size_t) gridIn) + ii) = Ith(init_state,spec);
    }
    
    N_VDestroy_Serial(init_state);
    
    try {
        cvode_mem = solver_setup (state, &pInS, AXL_react_diff);
        
        if (cvode_mem == NULL) throw runtime_error(string("Uncaught NULL return from solver_setup in diffusion setup."));
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    size_t tIDX = 0;
    
    if (tps[0] == 0) {
        for (size_t jj = 0; jj < (size_t) NV_LENGTH_S(state); jj++) dataPtr[jj] = Ith(state,jj);
        tIDX = 1;
    }
    
    for (; tIDX < nTps; tIDX++) {
        flag = CVode(cvode_mem, tps[tIDX], state, &t, CV_NORMAL);
        
        if (flag < 0) {
            N_VDestroy_Serial(state);
            CVodeFree(&cvode_mem);
            throw runtime_error(string("CVode error on diffusion solution."));
        }
        
        for (size_t jj = 0; jj < (size_t) NV_LENGTH_S(state); jj++) {
            dataPtr[tIDX*((size_t) NV_LENGTH_S(state)) + jj] = Ith(state,jj);
        }
    }
    
    /* Free vectors */
    N_VDestroy_Serial(state);
    CVodeFree(&cvode_mem);
}



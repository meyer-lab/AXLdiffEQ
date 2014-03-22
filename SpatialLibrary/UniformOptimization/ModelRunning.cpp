//
//  ModelRunning.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include <sundials/sundials_nvector.h>
#include <cvode/cvode.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <nlopt.hpp>
#include "ModelRunning.h"
#include "CVodeHelpers.h"

using namespace std;

struct inData {
    N_Vector fitt;           ///< Calculated model values.
    const double *pYmeas;    ///< pY measurement.
    const double *errorMeas; ///< Error for pY measurement.
};

// Calculate phosphorylation at time points measured
void calcProfile (N_Vector outData, N_Vector outStim, N_Vector outStimTot, struct rates params, double autocrine, double expression) {
    params.expression = expression;
    
    N_Vector init_state = N_VNew_Serial(Nspecies);
    N_Vector state = N_VNew_Serial(Nspecies);
    
    realtype t;
    
    void *cvode_mem = NULL;
    int flag;
    
    // Initialize state based on autocrine ligand
    try {
        initState(init_state, params, autocrine);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        N_VDestroy_Serial(init_state);
        throw;
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    for (unsigned int stimuli = 0; stimuli < NELEMS(Gass); stimuli++) {
        
        for (int xx = 0; xx < Nspecies; xx++) {
            Ith(state,xx) = Ith(init_state,xx);
        }
        
        Ith(state,0) += Gass[stimuli];
        t = 0;
        
        try {
            cvode_mem = solver_setup (state, &params, AXL_react);
            if (cvode_mem == NULL) throw runtime_error(string("Uncaught NULL return from solver_setup."));
        } catch (exception &e) {
            N_VDestroy_Serial(state);
            N_VDestroy_Serial(init_state);
            throw;
        }
        
        /* In loop, call CVode, print results, and test for error.
         Break out of loop when NOUT preset output times have been reached.  */
        
        Ith(outData,stimuli*NELEMS(times)) = pYcalc(state,params);
        
        for (unsigned int ii = 1; ii < NELEMS(times); ii++) {
            flag = CVode(cvode_mem, times[ii], state, &t, CV_NORMAL);
            
            if (flag < 0) {
                CVodeFree(&cvode_mem);
                N_VDestroy_Serial(state);
                N_VDestroy_Serial(init_state);
                throw runtime_error(string("Error from CVode on time course solve."));
            }
            
            Ith(outData,stimuli*NELEMS(times) + ii) = pYcalc(state,params);
        }
        
        /* Free integrator memory */
        CVodeFree(&cvode_mem);
        delete params.JacP;
    }
    
    
    /* We've got the initial state, so now run through the dose data */
    for (unsigned int stimuli = 0; stimuli < NELEMS(GassDose); stimuli++) {
        // Load the initial state (t = 0)
        for (int xx = 0; xx < Nspecies; xx++) {
            Ith(state,xx) = Ith(init_state,xx);
        }
        
        Ith(state,0) += GassDose[stimuli];
        t = 0;
        
        try {
            cvode_mem = solver_setup (state, &params, AXL_react);
            if (cvode_mem == NULL) throw runtime_error(string("Uncaught NULL return from solver_setup."));
        } catch (exception &e) {
            N_VDestroy_Serial(state);
            N_VDestroy_Serial(init_state);
            throw;
        }
        
        flag = CVode(cvode_mem, DoseTime, state, &t, CV_NORMAL);
        if (flag < 0) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(state);
            N_VDestroy_Serial(init_state);
            throw runtime_error(string("Error in CVode integration within calcProfile."));
        }
        
        Ith(outStim,stimuli) = pYcalc(state,params)/totCalc(state);
        Ith(outStimTot,stimuli) = totCalc(state);
        
        /* Free integrator memory */
        CVodeFree(&cvode_mem);
        delete params.JacP;
    }
    
    /* Free y and abstol vectors */
    N_VDestroy_Serial(state);
    N_VDestroy_Serial(init_state);
}

double errorOpt(unsigned n, const double *x, double *grad, void *data) {
    struct inData dataS;
    double xx = 0;
    
    memcpy(&dataS, data, sizeof(dataS));
    
    for (int ii = 0; ii < NV_LENGTH_S(dataS.fitt); ii++) {
        xx += pow((((double) Ith(dataS.fitt,ii) * x[0]) - dataS.pYmeas[ii]) / dataS.errorMeas[ii], 2);
    }
    
    return xx;
}

double errorFuncOpt (N_Vector fitt, const double *pYmeas, const double *errorMeas) {
    struct inData dataS;
    dataS.fitt = fitt;
    dataS.pYmeas = pYmeas;
    dataS.errorMeas = errorMeas;
    
    double ff = 0;
    vector<double> xx = {0.0000001};
    
    nlopt::opt opter = nlopt::opt(nlopt::algorithm::LN_COBYLA, 1);
    opter.set_lower_bounds(0);
    opter.set_upper_bounds(1E9);
    opter.set_min_objective(errorOpt,&dataS);
    opter.set_xtol_rel(1E-5);
    nlopt::result flag = opter.optimize(xx, ff);
    
    
    if (flag < 0) throw runtime_error(string("Error during error optimization step."));
    
    return ff;
}

double errorFuncFix (N_Vector fitt, const double *pYmeas, const double *errorMeas) {
    double xx = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(fitt); ii++)
        xx += pow((((double) Ith(fitt,ii)) - pYmeas[ii]) / errorMeas[ii], 2);
    
    return xx;
}

double calcErrorOneLine (struct rates inP, size_t cellLine, double autocrine) {
    N_Vector outData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector outStim = N_VNew_Serial(NELEMS(GassDose));
    N_Vector outStimTot = N_VNew_Serial(NELEMS(GassDose));
    N_Vector outDataAll = N_VNew_Serial(NELEMS(Gass)*NELEMS(times)*NfitCells);
    
    double error = 0;
    
    try {
        calcProfile (outData, outStim, outStimTot, inP, autocrine, inP.expression);
        
        error += errorFuncOpt (outData, &pY[cellLine*NELEMS(Gass)*NELEMS(times)], &pYerror[cellLine*NELEMS(Gass)*NELEMS(times)]);
        error += errorFuncOpt (outStim, pYdose[cellLine], DoseError[cellLine]);
        error += errorFuncFix (outStimTot, DoseTot[cellLine], DoseTotErr[cellLine]);
    } catch (exception &e) {
        N_VDestroy_Serial(outData);
        N_VDestroy_Serial(outStim);
        N_VDestroy_Serial(outStimTot);
        N_VDestroy_Serial(outDataAll);
        errorLogger(&e);
        
        return 1E6;
    }
    
    N_VDestroy_Serial(outData);
    N_VDestroy_Serial(outStim);
    N_VDestroy_Serial(outStimTot);
    N_VDestroy_Serial(outDataAll);
    
    return error;
}

double calcError (param_type inP) {
    struct rates params = Param(inP);
    
    N_Vector outData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector outStim = N_VNew_Serial(NELEMS(GassDose));
    N_Vector outStimTot = N_VNew_Serial(NELEMS(GassDose));
    N_Vector outDataAll = N_VNew_Serial(NELEMS(Gass)*NELEMS(times)*NfitCells);
    
    double error = 0;
    
    for (unsigned short ii = 0; ii < NfitCells; ii++) {
        try {
            calcProfile (outData, outStim, outStimTot, params, inP[15+ii], inP[15+NfitCells+ii]);
            
            error += errorFuncOpt (outData, &pY[ii*NELEMS(Gass)*NELEMS(times)], &pYerror[ii*NELEMS(Gass)*NELEMS(times)]);
            error += errorFuncOpt (outStim, pYdose[ii], DoseError[ii]);
            error += errorFuncFix (outStimTot, DoseTot[ii], DoseTotErr[ii]);
        } catch (exception &e) {
            N_VDestroy_Serial(outData);
            N_VDestroy_Serial(outStim);
            N_VDestroy_Serial(outStimTot);
            N_VDestroy_Serial(outDataAll);
            errorLogger(&e);
            
            return 1E6;
        }
    }
    
    N_VDestroy_Serial(outData);
    N_VDestroy_Serial(outStim);
    N_VDestroy_Serial(outStimTot);
    N_VDestroy_Serial(outDataAll);
    
    return error;
}

double errorFunc (double fitt, double pYmeas, double errorMeas) {
    return pow((((double)fitt) - pYmeas) / errorMeas, 2);
}

void calcErrorRef (param_type params, double *out, atomic<bool> *done) {
    *out = calcError(params);
    *done = true;
}



// Calculate the initial state by waiting a long time with autocrine Gas
void initState( N_Vector init, struct rates params, double autocrine) {
    endoImpair = 1.0;
    degImpair = 1.0;
    realtype t;
    
    for (int ii = 0; ii < Nspecies ; ii++) Ith(init,ii) = 0;
    
    Ith(init,0) = autocrine;
    Ith(init,1) = 1.5e5;
    Ith(init,7) = 1.5e5;
    
    void *cvode_mem = solver_setup (init, &params, AXL_react);
    if (cvode_mem == NULL) throw runtime_error(string("Error with solver setup in initState."));
    
    int flag = CVode(cvode_mem, autocrineT, init, &t, CV_NORMAL);
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Integration failure at initial condition."));
    }
    
    long int nJac;
    CVDlsGetNumJacEvals(cvode_mem, &nJac);
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    delete params.JacP;
}



/// Calculate phosphorylation at time points measured
void calcProfileSet (double *outData, double *tps, struct rates params, int nTps, double autocrine, double AXL, double GasStim, int frac) {
    params.expression = AXL;
    
    N_Vector state = N_VNew_Serial(Nspecies);
    
    realtype t; ///< Time position of the solver.
    
    void *cvode_mem = NULL;
    int flag;
    
    // Initialize state based on autocrine ligand
    try {
        initState(state, params, autocrine);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    Ith(state,0) += GasStim;
    t = 0;
    
    try {
        cvode_mem = solver_setup (state, &params, AXL_react);
        if (cvode_mem == NULL) throw runtime_error(string("Uncaught NULL return from solver_setup in calcProfileSet."));
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    size_t ii = 0;
    
    if (tps[0] == 0) {
        if (frac == 0) {
            outData[0] = pYcalc(state,params);
        } else if (frac == 1) {
            outData[0] = pYcalc(state,params) / totCalc(state);
        } else {
            outData[0] = totCalc(state);
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
            outData[ii] = pYcalc(state,params) / totCalc(state);
        } else {
            outData[ii] = totCalc(state);
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



void diffusionSolution(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn) {
    
    // Common
    realtype t = 0;
    void *cvode_mem = NULL;
    int flag;
    
    for (int ii = 0; ii < Nspecies; ii++) diffD[ii] = dIn[ii];
    
    // Create the parameter structure
    param_type pIn;
    struct rates pInS;
    for (size_t ii = 0; ii < pIn.size(); ii++) {
        pIn[ii] = params[ii];
    }
    pInS = Param(pIn);
    pInS.expression = AXLin;
    // Done creating parameter structure
    
    // Get initial state
    N_Vector init_state = N_VNew_Serial(Nspecies);
    
    // Initialize state based on autocrine ligand
    try {
        initState(init_state, pInS, autocrine);
    } catch (exception &e) {
        N_VDestroy_Serial(init_state);
        throw;
    }
    
    endoImpair = endoImpairIn;
    degImpair = degImpairIn;
    
    // Initialize full diffusion model
    N_Vector state = N_VNew_Serial(Nspecies * gridIn);
    
    for (size_t ii = 0; ii < (size_t) abs(gridIn); ii++) {
        for (size_t spec = 0; spec < Nspecies; spec++) Ith(state,spec*((size_t) gridIn) + ii) = Ith(init_state,spec);
    }
    
    N_VDestroy_Serial(init_state);
    
    for (size_t ii = 0; ii < (size_t) abs(gridIn); ii++) Ith(state,ii) = GasIn[ii];
    // Done initializing diffusion model
    
    try {
        cvode_mem = solver_setup (state, &pInS, AXL_react_diff);
        
        if (cvode_mem == NULL) {
            throw runtime_error(string("Uncaught NULL return from solver_setup in diffusion setup."));
        }
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    size_t tIDX;
    
    if (tps[0] == 0) {
        for (size_t jj = 0; jj < (size_t) NV_LENGTH_S(state); jj++) {
            dataPtr[jj] = Ith(state,jj);
            tIDX = 1;
        }
    } else tIDX = 0;
    
    
    for (; tIDX < (size_t) abs(nTps); tIDX++) {
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


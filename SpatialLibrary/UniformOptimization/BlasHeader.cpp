/**
 *   \file BlasHeader.cpp
 *   FiniteDiffGSL
 *
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

#include <iostream>
#include <thread>
#include <queue>
#include <iostream>
#include <fstream>
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <cvode/cvode_impl.h>
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_lapack.h>
#include <cvode/cvode_diag.h>
#include <nlopt.hpp>
#include "BlasHeader.h"

using namespace std;

const double internalFrac = 0.5;
const double fgMgConv = 135.2;

double diffD[Nspecies];

double endoImpair = 1; ///< Extent by which to impair endocytosis of Gas6-bound species.
double degImpair = 1;


struct inData {
    N_Vector fitt;           ///< Calculated model values.
    const double *pYmeas;    ///< pY measurement.
    const double *errorMeas; ///< Error for pY measurement.
};

// This takes the model state and calculates the amount of phosphorylated species
double pYcalc (N_Vector state, struct rates p) {
    double pY = Ith(state,1) + Ith(state,2) + Ith(state,3) + Ith(state,4) + 2*Ith(state,5) +
    internalFrac*(Ith(state,7) + Ith(state,8) + Ith(state,9) + Ith(state,10) + 2*Ith(state,11));
    
    pY *= p.scaleA;
    
    pY += 2*Ith(state,6) + internalFrac*(2*Ith(state,12));
    
    return pY;
}

// This takes the model state and calculates the total amount of receptor in a cell
double totCalc (N_Vector state) {
    double total = 0;
    
    for (int ii = 1; ii < 7; ii++) total += Ith(state,ii);
    for (int ii = 7; ii < 13; ii++) total += Ith(state,ii)*internalFrac;
    
    total += Ith(state,5);
    total += Ith(state,6);
    total += Ith(state,11)*internalFrac;
    total += Ith(state,12)*internalFrac;
    
    return total/fgMgConv;
}

int AXL_react(realtype t, N_Vector xIn , N_Vector dxdtIn,void *user_data) {
    struct rates *r = (struct rates *) user_data;
    realtype* x_d = NV_DATA_S(xIn);
    realtype* dxdt_d = NV_DATA_S(dxdtIn);
    
    // 0 Gas6   // 1 AXL   // 2 A1    // 3 A2
    // 4 A12    // 5 D1    // 6 D2    // 7 AXLi
    // 8 A1i    // 9 A2i   // 10 A12i // 11 D1i    // 12 D2i   // 13 Gasi
    
    dxdt_d[0] = 0;
    
    realtype dR1 = r->Binding1 * x_d[1] * x_d[0] - r->Unbinding1 * x_d[2];
    realtype dR2 = r->Binding2 * x_d[1] * x_d[0] - r->Unbinding2 * x_d[3];
    realtype dR3 = r->Binding2 * x_d[2] * x_d[0] - r->Unbinding2 * x_d[4];
    realtype dR4 = r->Binding1 * x_d[3] * x_d[0] - r->Unbinding1 * x_d[4];
    realtype dR5 = r->xFwd1 * x_d[1] * x_d[2] - r->xRev1 * x_d[5];
    realtype dR6 = r->xFwd2 * x_d[1] * x_d[3] - r->xRev2 * x_d[5];
    realtype dR7 = r->xFwd3 * x_d[1] * x_d[4] - r->xRev3 * x_d[6];
    realtype dR8 = r->xFwd4 * x_d[2] * x_d[2] - r->xRev4 * x_d[6];
    realtype dR9 = r->xFwd5 * x_d[3] * x_d[3] - r->xRev5 * x_d[6];
    realtype dR11 = r->xFwd6 * x_d[0] * x_d[5] - r->xRev6 * x_d[6];
    
    realtype dR32 = r->Binding1 * x_d[7] * x_d[13] / 623 - r->Unbinding1 * x_d[8];
    realtype dR33 = r->Binding2 * x_d[7] * x_d[13] / 623 - r->Unbinding2 * x_d[9];
    realtype dR34 = r->Binding2 * x_d[8] * x_d[13] / 623 - r->Unbinding2 * x_d[10];
    realtype dR35 = r->Binding1 * x_d[9] * x_d[13] / 623 - r->Unbinding1 * x_d[10];
    realtype dR36 = r->xFwd1 * x_d[7] * x_d[8] - r->xRev1 * x_d[11];
    realtype dR37 = r->xFwd2 * x_d[7] * x_d[9] - r->xRev2 * x_d[11];
    realtype dR38 = r->xFwd3 * x_d[7] * x_d[10] - r->xRev3 * x_d[12];
    realtype dR39 = r->xFwd4 * x_d[8] * x_d[8] - r->xRev4 * x_d[12]; // Checked
    realtype dR40 = r->xFwd5 * x_d[9] * x_d[9] - r->xRev5 * x_d[12]; // Checked
    realtype dR41 = r->xFwd6 * x_d[13] * x_d[11] / 623 - r->xRev6 * x_d[12]; // Checked
    
    dxdt_d[1] = - dR7 - dR6 - dR5 - dR1 - dR2 + r->expression; // AXL
    dxdt_d[2] = -2*(dR8) - dR5 + dR1 - dR3                   ; // AXLgas1
    dxdt_d[3] = -2*(dR9) - dR6 + dR2 - dR4                   ; // AXLgas2
    dxdt_d[4] = -dR7 + dR3 + dR4                             ; // AXLgas12
    dxdt_d[5] = -dR11 + dR6 + dR5                            ; // AXLdimer1
    dxdt_d[6] = dR11 + dR9 + dR8 + dR7                       ; // AXLdimer2
    
    dxdt_d[7]  = - dR38 - dR37 - dR36 - dR32 - dR33          ; // AXLi
    dxdt_d[8]  = -2*(dR39) - dR36 + dR32 - dR34              ; // AXLgas1i
    dxdt_d[9]  = -2*(dR40) - dR37 + dR33 - dR35              ; // AXLgas2i
    dxdt_d[10] = -dR38 + dR34 + dR35                         ; // AXLgas12i
    dxdt_d[11] = -dR41 + dR37 + dR36                          ; // AXLdimer1i
    dxdt_d[12] = dR41 + dR40 + dR39 + dR38                   ; // AXLdimer2i
    
    dxdt_d[13] = -dR41 - dR32 - dR33 - dR34 - dR35 - r->kDeg*x_d[13];
    
    
    dxdt_d[1] += -x_d[1]*(r->internalize + r->pYinternalize*r->scaleA) + r->kRec*(1-r->fElse)*x_d[7]*internalFrac; // Endocytosis, recycling
    dxdt_d[7] += x_d[1]*(r->internalize + r->pYinternalize*r->scaleA)/internalFrac - r->kRec*(1-r->fElse)*x_d[7] - r->kDeg*r->fElse*x_d[7]; // Endocytosis, recycling, degradation
    
    for (int ii = 2; ii < 6; ii++) {
        dxdt_d[ii]  += -x_d[ii]*(r->internalize + r->pYinternalize*r->scaleA)*endoImpair + r->kRec*(1-r->fElse)*x_d[ii+6]*internalFrac; // Endocytosis, recycling
        dxdt_d[ii+6] += x_d[ii]*(r->internalize + r->pYinternalize*r->scaleA)/internalFrac*endoImpair - r->kRec*(1-r->fElse)*x_d[ii+6] // Endocytosis, recycling
        - r->kDeg*r->fElse*x_d[ii+6]*degImpair; // Degradation
    }
    
    dxdt_d[6]  += -x_d[6]*(r->internalize + r->pYinternalize)*endoImpair + r->kRec*(1-r->fD2)*x_d[12]*internalFrac; // Endocytosis, recycling
    dxdt_d[12] += x_d[6]*(r->internalize + r->pYinternalize)/internalFrac*endoImpair - r->kRec*(1-r->fD2)*x_d[12] - r->kDeg*r->fD2*x_d[12]*degImpair; // Endocytosis, recycling, degradation
    
    return 0;
}

struct rates Param( param_type params) {
    struct rates out;
    
    out.Binding1 = params[0];
    out.Binding2 = params[1];
    out.Unbinding1 = params[2];
    out.Unbinding2 = params[3];
    out.xFwd1 = params[4];
    out.xRev1 = params[5];
    out.xFwd3 = params[6];
    out.xRev3 = params[7];
    out.internalize = params[8];
    out.pYinternalize = params[9];
    out.scaleA = params[10];
    out.kRec = params[11];
    out.kDeg = params[12];
    out.fElse = params[13];
    out.fD2 = params[14];
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

// Calculate the initial state by waiting a long time with autocrine Gas
int initState( N_Vector init, struct rates params, double autocrine, N_Vector abstol) {
    endoImpair = 1.0;
    degImpair = 1.0;
    void *cvode_mem = NULL;
    realtype t;
    
    for (int ii = 0; ii < Nspecies ; ii++) {
        Ith(init,ii) = 0;
    }
    Ith(init,0) = autocrine;
    Ith(init,1) = 1.5e5;
    Ith(init,7) = 1.5e5;
    
    
    int flag = solver_setup (&cvode_mem, init, &params, abstol, SELECTED_MODEL);
    if (flag == 1) {
        cout << "Error at solver setup" << endl;
        return(1);
    }
    
    flag = CVode(cvode_mem, autocrineT, init, &t, CV_NORMAL);
    if (check_flag(&flag, (char*)"CVode Initial Value", 1)) {
        CVodeFree(&cvode_mem);
        return(1);
    }
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    
    return 0;
}

int AXL_react_diff(realtype t, N_Vector xx , N_Vector dxxdt, void *user_data) {
    realtype* xx_d = NV_DATA_S(xx);
    realtype* dxxdt_d = NV_DATA_S(dxxdt);
    size_t pos, spec;
    size_t grid_size = (size_t) NV_LENGTH_S(xx)/Nspecies;
    double dRdRMaxRMaxR = maxR*maxR*(1.0/grid_size)*(1.0/grid_size);
    
    for (spec = 0; spec < Nspecies; spec++) {
        if (diffD[spec] == 0) continue; // Skip calculating diffusion if we're just going to multiply it by zero.
        
        for (size_t ii = 0; ii < (grid_size-1); ii++) {
            pos = spec*grid_size + ii;
            dxxdt_d[pos] = (-4.0*xx_d[pos] + (2.0-1.0/ii)*xx_d[pos - 1] + (2.0+1.0/ii)*xx_d[pos + 1])/2/dRdRMaxRMaxR;
        }
    }
    
    for (size_t ii = 0; ii < Nspecies; ii++) {
        if (diffD[ii] == 0) continue;
        
        // Take care of the inner boundary condition
        dxxdt_d[ii*grid_size] = 4*(xx_d[ii*grid_size + 1] - xx_d[ii*grid_size])/dRdRMaxRMaxR;
        
        // Outer boundary condition
        dxxdt_d[(ii+1)*grid_size - 1] = -4*(xx_d[(ii+1)*grid_size - 1] - xx_d[(ii+1)*grid_size - 2])/dRdRMaxRMaxR;
    }
    
    N_Vector reactIn = N_VNew_Serial(Nspecies);
    //if (check_flag((void *)reactIn, (char*)"N_VNew_Serial", 0)) return(1);
    N_Vector reactOut = N_VNew_Serial(Nspecies);
    //if (check_flag((void *)reactOut, (char*)"N_VNew_Serial", 0)) return(1);
    
    // Add in the reaction for each location
    for (size_t jj = 0; jj < grid_size; jj++) {
        
        for (size_t ii = 0; ii < Nspecies; ii++)
            Ith(reactIn,ii) = xx_d[ii*grid_size + jj];
        
        AXL_react(t,reactIn,reactOut,user_data);
        
        for (size_t ii = 0; ii < Nspecies; ii++) {
            // Also convert by diffusion coefficient
            dxxdt_d[ii*grid_size + jj] *= diffD[ii];
            
            // Add in reaction difference
            dxxdt_d[ii*grid_size + jj] += Ith(reactOut,ii);
        }
    }
    
    N_VDestroy_Serial(reactIn);
    N_VDestroy_Serial(reactOut);
    
    return 0;
}


void errorPrint(string In) {
    ofstream errfile;
    
    errfile.open ("error.txt", ios::app);
    errfile << In << endl;
    errfile.close();
}



/// Check function return value...
///  opt == 0 means SUNDIALS function allocates memory so check if
///           returned NULL pointer
///  opt == 1 means SUNDIALS function returns a flag so check if
///           flag >= 0
///  opt == 2 means function allocates memory so check if returned
///           NULL pointer
int check_flag(void *flagvalue, char *funcname, int opt) {
    ofstream errfile;
    
    int *errflag;
    
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        
        if (print_CV_err == 1) {
            cout << "SUNDIALS_ERROR: " << funcname << "() failed - returned NULL pointer" << endl;
        } else if (print_CV_err == 2) {
            errfile.open ("error.txt", ios::app);
            errfile << "SUNDIALS_ERROR: " << funcname << "() failed - returned NULL pointer" << endl;
            errfile.close();
        }
        return(1); }
    
    /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            if (print_CV_err == 1) {
                cout << "SUNDIALS_ERROR: " << funcname << "() failed with flag = " << *errflag << endl;
            } else if (print_CV_err == 2) {
                errfile.open ("error.txt", ios::app);
                errfile << "SUNDIALS_ERROR: " << funcname << "() failed with flag = " << *errflag << endl;
                errfile.close();
            }
            return(1);
        }
    }
    
    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        if (print_CV_err == 1) {
            cout << "MEMORY_ERROR: " << funcname << "() failed - returned NULL pointer" << endl;
        } else if (print_CV_err == 2) {
            errfile.open ("error.txt", ios::app);
            errfile << "MEMORY_ERROR: " << funcname << "() failed - returned NULL pointer" << endl;
            errfile.close();
        }
        return(1); }
    
    return(0);
}

/// Get and print some final statistics
/// \param *cvode_mem solver memory
void PrintFinalStats(void *cvode_mem) {
    long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
    int flag;
    
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(&flag, (char*)"CVodeGetNumSteps", 1);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_flag(&flag, (char*)"CVodeGetNumRhsEvals", 1);
    flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_flag(&flag, (char*)"CVodeGetNumLinSolvSetups", 1);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_flag(&flag, (char*)"CVodeGetNumErrTestFails", 1);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_flag(&flag, (char*)"CVodeGetNumNonlinSolvIters", 1);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_flag(&flag, (char*)"CVodeGetNumNonlinSolvConvFails", 1);
    
    flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
    check_flag(&flag, (char*)"CVDlsGetNumJacEvals", 1);
    flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
    check_flag(&flag, (char*)"CVDlsGetNumRhsEvals", 1);
    
    flag = CVodeGetNumGEvals(cvode_mem, &nge);
    check_flag(&flag, (char*)"CVodeGetNumGEvals", 1);
    
    printf("\nFinal Statistics:\n");
    printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
           nst, nfe, nsetups, nfeLS, nje);
    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
           nni, ncfn, netf, nge);
}

/// Error handler
/// \param error_code an int
void errorHandler (int error_code, const char *module, const char *function, char *msg, void *eh_data) {
    ofstream errfile;
    
    if (print_CV_err == 1) {
        cout << msg << endl;
    } else if (print_CV_err == 2) {
        errfile.open ("error.txt", ios::app);
        errfile << msg << endl;
        errfile.close();
    }
}

int solver_setup (void **cvode_mem, N_Vector init, void *params, N_Vector abstol, CVRhsFn f) {
    int flag;
    
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, (char*)"CVodeCreate", 0)) {
        CVodeFree(cvode_mem);
        return(1);
    }
    
    CVodeSetErrHandlerFn(*cvode_mem, &errorHandler, nullptr);
    
    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    flag = CVodeInit(*cvode_mem, f, 0.0, init);
    if (check_flag(&flag, (char*)"CVodeInit", 1)) {
        CVodeFree(cvode_mem);
        return(1);
    }
    
    // Pass along the parameter structure to the differential equations
    flag = CVodeSetUserData(*cvode_mem, params);
    
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(*cvode_mem, rel_tol, abstol);
    if (check_flag(&flag, (char*)"CVodeSVtolerances", 1)) {
        CVodeFree(cvode_mem);
        return(1);
    }
    
    /* Call CVDense to specify the CVDENSE dense linear solver
     * If the problem is large enough, it's better to use LAPACK
     */
    if (NV_LENGTH_S(init) > 20) {
        flag = CVLapackDense(*cvode_mem, (int) NV_LENGTH_S(init));
    } else {
        flag = CVDense(*cvode_mem, NV_LENGTH_S(init));
    }
    if (check_flag(&flag, (char*)"CVDense", 1)) {
        CVodeFree(cvode_mem);
        return(1);
    }
    
    CVodeSetMaxNumSteps(*cvode_mem, 2E6);
    CVodeSetMinStep(*cvode_mem,0);
    CVodeSetMaxStep(*cvode_mem,1000);
    
    return 0;
}


// Calculate phosphorylation at time points measured
int calcProfile (N_Vector outData, N_Vector outStim, N_Vector outStimTot, struct rates params, double autocrine, double expression) {
    params.expression = expression;
    
    N_Vector init_state = N_VNew_Serial(Nspecies);
    if (check_flag((void *)init_state, (char*)"N_VNew_Serial", 0)) return(1);
    N_Vector state = N_VNew_Serial(Nspecies);
    if (check_flag((void *)state, (char*)"N_VNew_Serial", 0)) return(1);
    
    realtype t;
    
    N_Vector abstol = N_VNew_Serial(Nspecies);
    if (check_flag((void *)abstol, (char*)"N_VNew_Serial", 0)) return(1);
    
    
    void *cvode_mem;
    int flag;
    
    cvode_mem = NULL;
    
    for (int ii = 0; ii < Nspecies; ii++) {
        Ith(abstol,ii) = abs_tol;
    }
    
    // Initialize state based on autocrine ligand
    flag = initState(init_state, params, autocrine, abstol);
    if (check_flag(&flag, (char*)"Initializer", 1)) {
        N_VDestroy_Serial(abstol);
        N_VDestroy_Serial(state);
        N_VDestroy_Serial(init_state);
        return(1);
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    for (unsigned int stimuli = 0; stimuli < NELEMS(Gass); stimuli++) {
        
        for (int xx = 0; xx < Nspecies; xx++) {
            Ith(state,xx) = Ith(init_state,xx);
        }
        
        Ith(state,0) += Gass[stimuli];
        t = 0;
        
        flag = solver_setup (&cvode_mem, state, &params, abstol, SELECTED_MODEL);
        if (flag == 1) {
            cout << "Error at solver setup" << endl;
            N_VDestroy_Serial(abstol);
            N_VDestroy_Serial(state);
            N_VDestroy_Serial(init_state);
            return(1);
        }
        
        /* In loop, call CVode, print results, and test for error.
         Break out of loop when NOUT preset output times have been reached.  */
        
        Ith(outData,stimuli*NELEMS(times)) = pYcalc(state,params);
        
        
        for (unsigned int ii = 1; ii < NELEMS(times); ii++) {
            flag = CVode(cvode_mem, times[ii], state, &t, CV_NORMAL);
            
            Ith(outData,stimuli*NELEMS(times) + ii) = pYcalc(state,params);
            
            
            if (check_flag(&flag, (char*)"CVode Time Course", 1)) {
                CVodeFree(&cvode_mem);
                N_VDestroy_Serial(abstol);
                N_VDestroy_Serial(state);
                N_VDestroy_Serial(init_state);
                return(1);
            }
        }
        
        
        /* Free integrator memory */
        CVodeFree(&cvode_mem);
    }
    
    
    /* We've got the initial state, so now run through the dose data */
    for (unsigned int stimuli = 0; stimuli < NELEMS(GassDose); stimuli++) {
        // Load the initial state (t = 0)
        for (int xx = 0; xx < Nspecies; xx++) {
            Ith(state,xx) = Ith(init_state,xx);
        }
        
        Ith(state,0) += GassDose[stimuli];
        t = 0;
        
        flag = solver_setup (&cvode_mem, state, &params, abstol, SELECTED_MODEL);
        if (flag == 1) {
            cout << "Error at solver setup" << endl;
            N_VDestroy_Serial(abstol);
            N_VDestroy_Serial(state);
            N_VDestroy_Serial(init_state);
            return(1);
        }
        
        flag = CVode(cvode_mem, DoseTime, state, &t, CV_NORMAL);
        if (check_flag(&flag, (char*)"CVode Dose Response", 1)) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(abstol);
            N_VDestroy_Serial(state);
            N_VDestroy_Serial(init_state);
            return(1);
        }
        
        Ith(outStim,stimuli) = pYcalc(state,params)/totCalc(state);
        Ith(outStimTot,stimuli) = totCalc(state);
        
        
        /* Free integrator memory */
        CVodeFree(&cvode_mem);
    }
    
    /* Free y and abstol vectors */
    N_VDestroy_Serial(abstol);
    N_VDestroy_Serial(state);
    N_VDestroy_Serial(init_state);
    
    return 0;
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

void errorOptPrint(vector<double> x, void *data) {
    struct inData dataS;
    double xx = 0;
    
    memcpy(&dataS, data, sizeof(dataS));
    
    
    for (int ii = 0; ii < NV_LENGTH_S(dataS.fitt); ii++) {
        xx += pow((((double) Ith(dataS.fitt,ii) * x[0]) - dataS.pYmeas[ii]) / dataS.errorMeas[ii], 2);
        
        
        cout << pow((((double) Ith(dataS.fitt,ii) * x[0]) - dataS.pYmeas[ii]) / dataS.errorMeas[ii], 2) << "\t\t" <<
        ((double) Ith(dataS.fitt,ii) * x[0]) << "\t" << dataS.pYmeas[ii] << "\t" << dataS.errorMeas[ii] << endl;
    }
}

double errorFuncOpt (N_Vector fitt, const double *pYmeas, const double *errorMeas, int print) {
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
    
    if (flag < 0) return 1E6;
    
    
    if (print == 1) {
        errorOptPrint(xx, &dataS);
        cout << "Opt: " << xx[0] << endl;
    }
    
    return ff;
}

double errorFuncFix (N_Vector fitt, const double *pYmeas, const double *errorMeas, int print) {
    double xx = 0;
    
    for (int ii = 0; ii < NV_LENGTH_S(fitt); ii++)
        xx += pow((((double) Ith(fitt,ii)) - pYmeas[ii]) / errorMeas[ii], 2);
    
    return xx;
}

double calcError (param_type inP) {
    struct rates params = Param(inP);
    
    N_Vector outData = N_VNew_Serial(NELEMS(Gass)*NELEMS(times));
    N_Vector outStim = N_VNew_Serial(NELEMS(GassDose));
    N_Vector outStimTot = N_VNew_Serial(NELEMS(GassDose));
    N_Vector outDataAll = N_VNew_Serial(NELEMS(Gass)*NELEMS(times)*NfitCells);
    
    double error = 0;
    
    for (unsigned short ii = 0; ii < NfitCells; ii++) {
        
        int flag = calcProfile (outData, outStim, outStimTot, params, inP[15+ii], inP[15+NfitCells+ii]) ;
        
        if (flag != 0) {
            error = 1e10;
            break;
        }
        
        error += errorFuncOpt (outData, &pY[ii*NELEMS(Gass)*NELEMS(times)], &pYerror[ii*NELEMS(Gass)*NELEMS(times)], 0);
        error += errorFuncOpt (outStim, pYdose[ii], DoseError[ii], 0);
        error += errorFuncFix (outStimTot, DoseTot[ii], DoseTotErr[ii], 0);
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




/* MATLAB INTERFACE ELEMENTS */


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
    
    const int nThreads = 8; 		/**      */
    queue<queT> runThese;
    std::thread t[nThreads];
    atomic<bool> done[nThreads];
    queT in;
    
    for (size_t ii = 0; ii < nIn; ii++) {
        param_type pInSlice;
        
        for (size_t jj = 0; jj < pInSlice.size(); jj++) {
            pInSlice[jj] = pIn[(size_t) ii*pInSlice.size() + jj];
        }
        
        in.In = pInSlice;
        in.out = &dataPtr[ii];
        runThese.push(in);
     }

    
    for (int ii = 0; ii < nThreads; ii++) {
        if (runThese.size() == 0) {
            break;
        }
        t[ii] = std::thread(calcErrorRef,runThese.front().In,runThese.front().out, &done[ii]);
        runThese.pop();
    }
    
    while (runThese.size() > 0) {
        for (int ii = 0; ii < nThreads; ii++) {
            if (runThese.size() == 0) {
                break;
            }
            
            if (done[ii] == true) {
                t[ii].join();
                done[ii] = false;
                t[ii] = std::thread(calcErrorRef,runThese.front().In,runThese.front().out, &done[ii]);
                runThese.pop();
            }
        }
    }
    
    if (nIn > nThreads) {
        for (int ii = 0; ii < nThreads; ii++) {
            t[ii].join();
        }
    } else {
        for (int ii = 0; ii < nIn; ii++) {
            t[ii].join();
        }
    }
    
    return 0;
}

/// Function for fitting a parameter set from R
extern "C" void rEntry(double *dataPtr, const double *pIn) {
    atomic<bool> done;
    
    param_type pInSlice;
    for (size_t jj = 0; jj < pInSlice.size(); jj++) {
        pInSlice[jj] = pIn[jj];
    }
    
    calcErrorRef(pInSlice,dataPtr, &done);
}

/// Calculate phosphorylation at time points measured
int calcProfileSet (double *outData, double *tps, struct rates params, int nTps, double autocrine, double AXL, double GasStim, int frac) {
    params.expression = AXL;
    
    N_Vector state = N_VNew_Serial(Nspecies);
    if (check_flag((void *)state, (char*)"N_VNew_Serial", 0)) return(1);
    
    realtype t; ///< Time position of the solver.
    
    N_Vector abstol = N_VNew_Serial(Nspecies);
    if (check_flag((void *)abstol, (char*)"N_VNew_Serial", 0)) return(1);
    
    void *cvode_mem = NULL;
    int flag;
    
    for (int ii = 0; ii < Nspecies; ii++) {
        Ith(abstol,ii) = abs_tol;
    }
    
    // Initialize state based on autocrine ligand
    flag = initState(state, params, autocrine, abstol);
    
    if (check_flag(&flag, (char*)"Initializer", 1)) {
        N_VDestroy_Serial(abstol);
        N_VDestroy_Serial(state);
        return(1);
    }
    
    /* We've got the initial state, so now run through the kinetic data */
    Ith(state,0) += GasStim;
    t = 0;
    
    flag = solver_setup (&cvode_mem, state, &params, abstol, SELECTED_MODEL);
    if (flag == 1) {
        cout << "Error at solver setup" << endl;
        N_VDestroy_Serial(abstol);
        N_VDestroy_Serial(state);
        return(1);
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
    
    for (; ii < nTps; ii++) {
        flag = CVode(cvode_mem, tps[ii], state, &t, CV_NORMAL);
        
        if (frac == 0) {
            outData[ii] = pYcalc(state,params);
        } else if (frac == 1) {
            outData[ii] = pYcalc(state,params) / totCalc(state);
        } else {
            outData[ii] = totCalc(state);
        }
        
        if (check_flag(&flag, (char*)"CVode Time Course", 1)) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(abstol);
            N_VDestroy_Serial(state);
            return(1);
        }
    }
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    /* Free y and abstol vectors */
    N_VDestroy_Serial(abstol);
    N_VDestroy_Serial(state);
    
    return 0;
}

extern "C" int calcProfileMatlab(double *dataPtr, double *params, double *tps, int nTps, double autocrine, double AXL, double GasStim, int frac) {
    param_type pIn;
    
    for (size_t ii = 0; ii < pIn.size(); ii++)
        pIn[ii] = params[ii];

    int flag = calcProfileSet (dataPtr, tps, Param(pIn), nTps, autocrine, AXL, GasStim, frac);
    
    return flag;
}

extern "C" int matlabDiffTPS(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn) {
    // Common
    realtype t = 0;
    void *cvode_mem = NULL;
    int flag;
    
    for (int ii = 0; ii < Nspecies; ii++)
        diffD[ii] = dIn[ii];
    
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
    if (check_flag((void *)init_state, (char*)"N_VNew_Serial", 0)) return(1);
    N_Vector abstol = N_VNew_Serial(Nspecies);
    if (check_flag((void *)abstol, (char*)"N_VNew_Serial", 0)) return(1);
    
    
    for (int ii = 0; ii < Nspecies; ii++)
        Ith(abstol,ii) = abs_tol;
    
    // Initialize state based on autocrine ligand
    flag = initState(init_state, pInS, autocrine, abstol);
    if (check_flag(&flag, (char*)"Initializer", 1)) {
        N_VDestroy_Serial(abstol);
        N_VDestroy_Serial(init_state);
        return(1);
    }
    // Done getting initial state
    
    endoImpair = endoImpairIn;
    degImpair = degImpairIn;
    
    // Initialize full diffusion model
    N_Vector state = N_VNew_Serial(Nspecies * gridIn);
    if (check_flag((void *)state, (char*)"N_VNew_Serial", 0)) return(1);
    N_VDestroy_Serial(abstol);
    abstol = N_VNew_Serial(Nspecies * gridIn);
    if (check_flag((void *)abstol, (char*)"N_VNew_Serial", 0)) return(1);
    for (int ii = 0; ii < Nspecies*gridIn; ii++) {
        Ith(abstol,ii) = abs_tol;
    }
    
    for (size_t ii = 0; ii < gridIn; ii++) {
        for (size_t spec = 0; spec < Nspecies; spec++) {
            Ith(state,spec*((size_t) gridIn) + ii) = Ith(init_state,spec);
        }
    }
    
    for (size_t ii = 0; ii < gridIn; ii++) {
        Ith(state,ii) = GasIn[ii];
    }
    // Done initializing diffusion model
    
    
    flag = solver_setup (&cvode_mem, state, &pInS, abstol, AXL_react_diff);
    if (flag == 1) {
        cout << "Error at solver setup" << endl;
        N_VDestroy_Serial(abstol);
        N_VDestroy_Serial(state);
        return(1);
    }
    
    size_t tIDX;
    
    if (tps[0] == 0) {
        for (size_t jj = 0; jj < NV_LENGTH_S(state); jj++) {
            dataPtr[jj] = Ith(state,jj);
            tIDX = 1;
        }
    } else {
        tIDX = 0;
    }
    
    for (; tIDX < nTps; tIDX++) {
        flag = CVode(cvode_mem, tps[tIDX], state, &t, CV_NORMAL);
        
        if (check_flag(&flag, (char*)"CVode", 1)) return(1);
        
        for (size_t jj = 0; jj < NV_LENGTH_S(state); jj++) {
            dataPtr[tIDX*((size_t) NV_LENGTH_S(state)) + jj] = Ith(state,jj);
        }
    }
    
    /* Free y and abstol vectors */
    N_VDestroy_Serial(state);
    N_VDestroy_Serial(abstol);
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    
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
    if (check_flag((void *)state, (char*)"N_VNew_Serial", 0)) return(1);
    
    for (size_t time = 0; time < nTps; time++) {
        for (size_t gridP = 0; gridP < gridIn; gridP++) {
            for (size_t spec = 0; spec < Nspecies; spec++)
                Ith(state,spec) = dataPtrTemp[time*(gridIn*Nspecies) + spec*gridIn + gridP];
            
            if (frac == 0) {
                dataPtr[time*gridIn + gridP] = pYcalc(state, pInS);
            } else if (frac == 1) {
                dataPtr[time*gridIn + gridP] = pYcalc(state, pInS) / totCalc(state);
            } else {
                dataPtr[time*gridIn + gridP] = totCalc(state);
            }
        }
    }
    
    
    return 0;
}

extern "C" int matlabDiffTPS_pYavg(double *dataPtr, double AXLin, double *GasIn, int gridIn, double autocrine, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn, int frac) {
    
    double dataPtrTemp[gridIn*nTps];
    
    int flag = matlabDiffTPS_pY(dataPtrTemp, AXLin, GasIn, gridIn, autocrine, params, tps, nTps, dIn, endoImpairIn, degImpairIn, frac);
    if (flag == 1) return(1);
    
    double summ;
    
    for (size_t time = 0; time < nTps; time++) {
        summ = 0;
        
        for (size_t gridP = 0; gridP < gridIn; gridP++) {
            summ += dataPtrTemp[time*gridIn + gridP]/(gridIn)*(gridP);
        }
        
        dataPtr[time] = 2*summ/gridIn;
    }
    
    return 0;
}
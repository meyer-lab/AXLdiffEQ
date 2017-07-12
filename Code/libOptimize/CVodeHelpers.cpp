//
//  CVodeHelpers.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/11/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <string>
#include <cvode/cvode_dense.h>     /* prototype for CVDense */
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "CVodeHelpers.h"
#include "ModelRunning.h"
#include <sundials/sundials_dense.h>

#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

using namespace std;

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

static void errorHandler(int error_code, const char *module, const char *function, char *msg, void *) {
    if (error_code == CV_WARNING) return;
    
    stringstream OutMesg;

    OutMesg << "Internal CVode error in " << function << endl;
    OutMesg << msg << endl;
    OutMesg << "In module: " << module << endl;
    OutMesg << "Error code: " << error_code << endl;
    
    errorLogger(OutMesg);
}


void* solver_setup (N_Vector init, void *params, double abstolIn, double reltolIn, CVRhsFn f) {
    void *cvode_mem = NULL;
    
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeCreate in solver_setup."));
    }
    
    CVodeSetErrHandlerFn(cvode_mem, &errorHandler, NULL);
    

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    if (CVodeInit(cvode_mem, f, 0.0, init) < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeInit in solver_setup."));
    }
    
    N_Vector abbstol = N_VNew_Serial(NV_LENGTH_S(init));
    for (int ii = 0; ii < NV_LENGTH_S(init); ii++) {
        Ith(abbstol,ii) = abstolIn;
    }
    
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    if (CVodeSVtolerances(cvode_mem, reltolIn, abbstol) < 0) {
        N_VDestroy_Serial(abbstol);
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeSVtolerances in solver_setup."));
    }
    N_VDestroy_Serial(abbstol);
    
    // Call CVDense to specify the CVDENSE dense linear solver
    if (CVDense(cvode_mem, (int) NV_LENGTH_S(init)) < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVDense in solver_setup."));
    }
    
    // Pass along the parameter structure to the differential equations
    if (CVodeSetUserData(cvode_mem, params) < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeSetUserData in solver_setup."));
    }

    CVodeSetMaxNumSteps(cvode_mem, 5000);
    
    return cvode_mem;
}

void solverReset (void *cvode_mem, N_Vector init) {
    if (CVodeReInit(cvode_mem, 0.0, init) < 0) throw runtime_error(string("Error at CVode reinit."));
}

void* solver_setup (N_Vector init, void *params, CVRhsFn f) {
    return solver_setup (init, params, 1E-3, 1E-6, f);
}


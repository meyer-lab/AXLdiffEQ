//
//  CVodeHelpers.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/11/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include "nvector_serial.h"  /* serial N_Vector types, fcts., macros */
#include "cvode.h"             /* prototypes for CVODE fcts., consts. */
#include <string>
#include "cvode_dense.h"       /* prototype for CVDense */
#include "cvode_band.h"
#include <sstream>
#include "CVodeHelpers.h"
#include "ReactionCode.h"
#include "sundials_dense.h"

#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

using namespace std;

void errorHandler(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
    //stringstream OutMesg;
    
    //OutMesg << "Internal CVode error in " << function << endl;
    //OutMesg << msg << endl;
    //OutMesg << "In module: " << module << endl;
    //OutMesg << "Error code: " << error_code << endl;
    
    //throw runtime_error(OutMesg.str());
}

void* solver_setup (N_Vector init, void *params, double abstolIn, double reltolIn, CVRhsFn f) {
    int flag;
    void *cvode_mem = NULL;
    
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeCreate in solver_setup."));
    }
    
    CVodeSetErrHandlerFn(cvode_mem, &errorHandler, nullptr);
    
    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, f, 0.0, init);
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeInit in solver_setup."));
    }
    
    N_Vector abbstol = N_VNew_Serial(NV_LENGTH_S(init));
    for (int ii = 0; ii < NV_LENGTH_S(init); ii++) {
        Ith(abbstol,ii) = abstolIn;
    }
    
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltolIn, abbstol);
    N_VDestroy_Serial(abbstol);
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeSVtolerances in solver_setup."));
    }
    
    // Call CVDense to specify the CVDENSE dense linear solver
    if (NV_LENGTH_S(init) > 20) {
        flag = CVBand(cvode_mem, (int) NV_LENGTH_S(init), 2, 2);
    } else {
        flag = CVDense(cvode_mem, NV_LENGTH_S(init));
    }
    
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVDense in solver_setup."));
    }
    
    // Pass along the parameter structure to the differential equations
    flag = CVodeSetUserData(cvode_mem, params);
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        throw runtime_error(string("Error calling CVodeSetUserData in solver_setup."));
    }
    
    CVodeSetMaxConvFails(cvode_mem, 50);
    CVodeSetMaxNumSteps(cvode_mem, 2E6);
    
    return cvode_mem;
}

void solverReset (void *cvode_mem, N_Vector init) {
    int flag = CVodeReInit(cvode_mem, 0.0, init);
    if (flag < 0) throw runtime_error(string("Error at CVode reinit."));
}

void* solver_setup (N_Vector init, void *params, CVRhsFn f) {
    return solver_setup (init, params, 1E-4, 1E-6, f);
}

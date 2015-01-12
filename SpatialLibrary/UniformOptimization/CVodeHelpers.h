//
//  CVodeHelpers.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/11/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__CVodeHelpers__
#define __UniformOptimization__CVodeHelpers__

#include "nvector/nvector_serial.h"  /* serial N_Vector types, fcts., macros */
#include "cvode/cvode.h"
#include "cvode/cvode_direct.h"


#define print_CV_err 0

#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
void* solver_setup (N_Vector, void *, double, double, CVRhsFn);
void* solver_setup (N_Vector, void *, CVRhsFn);
void solverReset (void *, N_Vector);
void PrintFinalStats(void *);

void errorLogger (std::exception *);
void errorLogger (std::stringstream &);

#endif /* defined(__UniformOptimization__CVodeHelpers__) */

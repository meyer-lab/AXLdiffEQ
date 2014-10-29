/**
 *   \file BlasHeader.h
 *   FiniteDiffGSL
 *
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

#ifndef __UniformOptimization__HelperFunctions__
#define __UniformOptimization__HelperFunctions__

#ifdef __cplusplus
extern "C" {
#endif

int matlabEntry(double * , double * , int);
int calcProfileMatlab(double *, double *, double *, int, double, int);
//int matlabDiffTPS(double *dataPtr, double AXLin, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn);
//int matlabDiffTPS_pY(double *, double , double *, int , double , double *, double *, int, double *, double, double, int);
//int matlabDiffTPS_pYavg(double *, double, double *, int, double , double *, double *, int, double *, double, double, int);


    
#ifdef __cplusplus
}
#endif

#endif /* defined(__UniformOptimization__HelperFunctions__) */

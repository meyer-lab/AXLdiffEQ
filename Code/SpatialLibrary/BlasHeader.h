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

    double pyEntry(double *pIn);
    int calcProfileMatlab(double *, double *, double *, double *, double *, double *, unsigned int, double, double *);
    void pyEntryVec(double *pIn, double *pOut, int n);
    
    int diffCalc(double *, double *, double *, double *, double *, unsigned int, double *, double *, unsigned int, double *);
    int getNspecies ();
    void U87pyEntryVec(double *pIn, double *pOut, int n);
    
#ifdef __cplusplus
}
#endif

#endif /* defined(__UniformOptimization__HelperFunctions__) */

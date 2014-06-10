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
int calcProfileMatlab(double *, double *, double *, int, double, double, double, int);
int matlabDiffTPS(double *, double , double *, int , double , double *, double *, int, double *, double, double);
int matlabDiffTPS_pY(double *, double , double *, int , double , double *, double *, int, double *, double, double, int);
int matlabDiffTPS_pYavg(double *, double, double *, int, double , double *, double *, int, double *, double, double, int);
void rEntry(double *, const double *);
double calcErrorOneCellLine (int, const double *);
int matlabEntryWithSi(double *, double *, int);
int matlabEntryA549(double *, double *, int);
int matlabEntryA549VaryEndo(double *dataPtr, double *pIn, int nIn);
    
#ifdef __cplusplus
}
#endif

#endif /* defined(__UniformOptimization__HelperFunctions__) */

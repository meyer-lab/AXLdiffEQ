/**
 *   \file BlasHeader.h
 *   FiniteDiffGSL
 *
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

#ifndef __UniformOptimization__HelperFunctions__
#define __UniformOptimization__HelperFunctions__

extern "C" int matlabEntry(double * , double * , int);
extern "C" int calcProfileMatlab(double *, double *, double *, int, double, double, double, int);
extern "C" int matlabDiffTPS(double *, double , double *, int , double , double *, double *, int, double *, double, double);
extern "C" int matlabDiffTPS_pY(double *, double , double *, int , double , double *, double *, int, double *, double, double, int);
extern "C" int matlabDiffTPS_pYavg(double *, double, double *, int, double , double *, double *, int, double *, double, double, int);
extern "C" void rEntry(double *, const double *);
extern "C" double calcErrorOneCellLine (int, const double *);
extern "C" int matlabEntryWithSi(double *, double *, int);
extern "C" int matlabEntryA549(double *, double *, int);

#endif /* defined(__UniformOptimization__HelperFunctions__) */

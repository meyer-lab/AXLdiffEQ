//
//  ModelRunning.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ModelRunning__
#define __UniformOptimization__ModelRunning__

#include "ReactionCode.h"
#include <atomic>

#define autocrineT 10000
#define print_CV_err 1
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))


static const double times[] = {60, 240}; ///< Times of kinetic measurements.
static const double Gass[] = {64, 16, 4, 1, 0.25, 0}; ///< Kinetic Gas6 doses.

// Wrapping is outermost cell line, then Gas, then time
static const double pY[6][2] = { ///< pY measurements on short time scales.
    {10.75952427, 8.305264139},
    {7.390399159, 7.056438019},
    {7.144036441, 7.680851079},
    {4.570826833, 8.184089069},
    {6.107714557, 7.204021903},
    {7.535575387, 7.535575387}};

static const double pYerror[6][2] = { ///< Error for short time scale pY measurements.
    {1.74431775, 2.100723242},
    {1.267611, 1.260108508},
    {0.898008437, 1.680415875},
    {1.521844479, 0.871927763},
    {0.932623012, 0.563182873},
    {0.812417951, 0.812417951}};

static const double tot[6][2] = {
    {242.9578093, 227.1923827},
    {221.8099578, 236.6570748},
    {213.0225396, 254.8630804},
    {184.0914948, 243.3181911},
    {189.8325909, 223.5550426},
    {188.5457582, 188.5457582}}; // A549

static const double totError[6][2] = {
    {12.30, 9.32},
    {13.34, 9.17},
    {17.34, 15.91},
    {10.93, 14.38},
    {9.08, 13.22},
    {5.83, 5.83}}; // A549

double calcError (param_type);
void errorLogger (std::exception *);
void*initState(N_Vector, struct rates *);
void diffusionSolution(double *dataPtr, double AXLin, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn);
void calcErrorRef (param_type, double *, std::atomic<bool> *);
double calcErrorOneLine (struct rates, size_t, double);
void errorLogger (std::stringstream &);
void calcProfile (N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, struct rates *, double, double);
void calcProfileSet (double *outData, double *tps, struct rates *params, int nTps, double GasStim, int frac);

#endif /* defined(__UniformOptimization__ModelRunning__) */

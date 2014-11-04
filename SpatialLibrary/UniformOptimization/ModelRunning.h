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

#define autocrineT 10000
#define print_CV_err 0
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
    {3443.11, 3219.69},
    {3143.41, 3353.82},
    {3018.88, 3611.82},
    {2608.88, 3448.21},
    {2690.24, 3168.14},
    {2672.00, 2672.00}}; // A549

static const double totError[6][2] = {
    {174.38, 132.10},
    {189.03, 129.93},
    {245.75, 225.42},
    {154.89, 203.72},
    {128.72, 187.34},
    {82.62, 82.62}}; // A549

double calcError (struct rates);
void errorLogger (std::exception *);
void*initState(N_Vector, struct rates *);
void diffusionSolution(double *dataPtr, double *GasIn, int gridIn, double *params, double *tps, int nTps, double *dIn, double endoImpairIn, double degImpairIn);
double calcErrorOneLine (struct rates, size_t, double);
void errorLogger (std::stringstream &);
void calcProfile (N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, struct rates *, double, double);
void calcProfileSet (double *outData, double *tps, struct rates *params, int nTps, double GasStim, int frac);


#endif /* defined(__UniformOptimization__ModelRunning__) */

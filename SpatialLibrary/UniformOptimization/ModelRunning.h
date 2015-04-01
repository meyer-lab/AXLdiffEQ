//
//  ModelRunning.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ModelRunning__
#define __UniformOptimization__ModelRunning__

#include "sundials_nvector.h"
#include "cvode_impl.h"
#include "cvode.h"
#include "sundials_dense.h"

#define autocrineT 10000
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))

#define Nspecies 13
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */

#define Nparams 11

//#define numParams 14

struct rates {
    double Binding1;   ///< Forward binding rate for Ig1
    double Binding2;   ///< Forward binding rate for Ig2
    double Unbinding1; ///< Reverse binding rate for Ig1
    double Unbinding2; ///< Reverse binding rate for Ig2
    double xFwd1;      ///< Reaction 1 forward rate.
    double xRev1;      ///< Reaction 1 reverse rate.
    double xFwd2;      ///< Reaction 2 forward rate.
    double xRev2;      ///< Reaction 2 reverse rate.
    double xFwd3;      ///< Reaction 3 forward rate.
    double xRev3;      ///< Reaction 3 reverse rate.
    double xFwd4;      ///< Reaction 4 forward rate.
    double xRev4;      ///< Reaction 4 reverse rate.
    double xFwd5;      ///< Reaction 5 forward rate.
    double xRev5;      ///< Reaction 5 reverse rate.
    double xFwd6;      ///< Reaction 6 forward rate.
    double xRev6;      ///< Reaction 6 reverse rate.
    double expression; ///< AXL expression rate.
    double internalize;///< Non-pY species internalization rate.
    double pYinternalize;///< pY species internalization rate.
    double kRec;       ///< Recycling rate.
    double kDeg;       ///< Degradation rate.
    double fElse;      ///< Recycling fraction for non-D2 species.
    double fD2;        ///< Recycling fraction for D2.
    int pD1;
    double internalFrac;
    double internalV;
    double autocrine;
    double gasCur;
};

static const double times[2] = {60, 240}; ///< Times of kinetic measurements.
static const double Gass[6] = {64, 16, 4, 1, 0.25, 0}; ///< Kinetic Gas6 doses.
static const double kTPS[5] = {0, 0.5, 1, 5, 10};


static const double pYk[5] = {4.1, 3.6, 5.3, 11.3, 11.6};
static const double pYkErr[5] = {1.0, 1.3, 1.1, 0.7, 0.4};

// Wrapping is outermost cell line, then Gas, then time
static const double pY[6][2] = { ///< pY measurements on short time scales.
    {10.75952427, 8.305264139},
    {7.390399159, 7.056438019},
    {7.144036441, 7.680851079},
    {4.570826833, 8.184089069},
    {6.107714557, 7.204021903},
    {7.535575387, 7.535575387}};

static const double pYerror[6][2] = { ///< Error for short time scale pY measurements.
    {1.74431775,  2.100723242},
    {1.267611,    1.260108508},
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

static const double surf[6][2] = {
    {0.206, 0.239},
    {0.274, 0.316},
    {0.281, 0.251},
    {0.220, 0.302},
    {0.256, 0.281},
    {0.257, 0.337}}; // A549

static const double surfError[6][2] = {
    {0.043, 0.015},
    {0.047, 0.037},
    {0.032, 0.030},
    {0.025, 0.036},
    {0.044, 0.035},
    {0.030, 0.023}}; // A549


double calcError (struct rates);
void*initState(N_Vector, struct rates *);
void calcProfile (N_Vector, N_Vector, N_Vector, N_Vector, N_Vector, struct rates *, double, double);
void calcProfileSet (double *outData, double *tps, struct rates *params, int nTps, double GasStim, int frac);
int AXL_react(double, N_Vector, N_Vector, void *);
struct rates Param(const double *);
int AXL_jac (long int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
double calcErrorMer (struct rates inP);

#endif /* defined(__UniformOptimization__ModelRunning__) */

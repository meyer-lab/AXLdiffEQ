//
//  ReactionCode.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/11/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ReactionCode__
#define __UniformOptimization__ReactionCode__

#include "nvector/nvector_serial.h"  /* serial N_Vector types, fcts., macros */
#include <array>
#include <vector>

#define Nspecies 14
#define maxR 1.0
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */

extern double diffD[Nspecies];
extern double endoImpair; ///< Extent by which to impair endocytosis of Gas6-bound species.
extern double degImpair;

#define numParams 13

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
    double internalFrac;
    double internalV;
    double autocrine;
};

struct diffRates {
    struct rates params;
    N_Vector reactIn;
    N_Vector reactOut;
};

int AXL_react(double, N_Vector, N_Vector, void *);
int AXL_react_diff(double, N_Vector, N_Vector, void *);
double pYcalc (N_Vector, struct rates *);
double totCalc (N_Vector, struct rates *);
struct rates Param(double*);
double surfAXL (N_Vector);

#endif /* defined(__UniformOptimization__ReactionCode__) */

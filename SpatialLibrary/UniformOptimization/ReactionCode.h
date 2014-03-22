//
//  ReactionCode.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/11/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ReactionCode__
#define __UniformOptimization__ReactionCode__

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <array>
#include <vector>
#include <cvode/cvode_dense.h>
#include <cppad/cppad.hpp>


#define Nspecies 14
#define maxR 1.0
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */

extern double diffD[Nspecies];
extern double endoImpair; ///< Extent by which to impair endocytosis of Gas6-bound species.
extern double degImpair;
extern const double internalFrac;

typedef std::array< double , 19 > param_type;

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
    double scaleA;     ///< Signaling capacity of non-D2 species.
    double expression; ///< AXL expression rate.
    double internalize;///< Non-pY species internalization rate.
    double pYinternalize;///< pY species internalization rate.
    double kRec;       ///< Recycling rate.
    double kDeg;       ///< Degradation rate.
    double fElse;      ///< Recycling fraction for non-D2 species.
    double fD2;        ///< Recycling fraction for D2.
    CppAD::ADFun<double> *JacP;
};

int AXL_react(realtype, N_Vector, N_Vector, void *);
int AXL_react_diff(realtype, N_Vector, N_Vector, void *);
double pYcalc (N_Vector, struct rates);
double totCalc (N_Vector);
struct rates Param(param_type);

template <class Type>
inline void AXL_react_new(std::vector<Type> &x_d, std::vector<Type> &dxdt_d, void *user_data) {
    struct rates *r = (struct rates *) user_data;
    
    // 0 Gas6   // 1 AXL   // 2 A1    // 3 A2
    // 4 A12    // 5 D1    // 6 D2    // 7 AXLi
    // 8 A1i    // 9 A2i   // 10 A12i // 11 D1i    // 12 D2i   // 13 Gasi
    
    dxdt_d[0] = 0;
    
    Type dR1 = r->Binding1 * x_d[1] * x_d[0] - r->Unbinding1 * x_d[2];
    Type dR2 = r->Binding2 * x_d[1] * x_d[0] - r->Unbinding2 * x_d[3];
    Type dR3 = r->Binding2 * x_d[2] * x_d[0] - r->Unbinding2 * x_d[4];
    Type dR4 = r->Binding1 * x_d[3] * x_d[0] - r->Unbinding1 * x_d[4];
    Type dR5 = r->xFwd1 * x_d[1] * x_d[2] - r->xRev1 * x_d[5];
    Type dR6 = r->xFwd2 * x_d[1] * x_d[3] - r->xRev2 * x_d[5];
    Type dR7 = r->xFwd3 * x_d[1] * x_d[4] - r->xRev3 * x_d[6];
    Type dR8 = r->xFwd4 * x_d[2] * x_d[2] - r->xRev4 * x_d[6];
    Type dR9 = r->xFwd5 * x_d[3] * x_d[3] - r->xRev5 * x_d[6];
    Type dR11 = r->xFwd6 * x_d[0] * x_d[5] - r->xRev6 * x_d[6];
    
    Type dR32 = r->Binding1 * x_d[7] * x_d[13] / 623 - r->Unbinding1 * x_d[8];
    Type dR33 = r->Binding2 * x_d[7] * x_d[13] / 623 - r->Unbinding2 * x_d[9];
    Type dR34 = r->Binding2 * x_d[8] * x_d[13] / 623 - r->Unbinding2 * x_d[10];
    Type dR35 = r->Binding1 * x_d[9] * x_d[13] / 623 - r->Unbinding1 * x_d[10];
    Type dR36 = r->xFwd1 * x_d[7] * x_d[8] - r->xRev1 * x_d[11];
    Type dR37 = r->xFwd2 * x_d[7] * x_d[9] - r->xRev2 * x_d[11];
    Type dR38 = r->xFwd3 * x_d[7] * x_d[10] - r->xRev3 * x_d[12];
    Type dR39 = r->xFwd4 * x_d[8] * x_d[8] - r->xRev4 * x_d[12]; // Checked
    Type dR40 = r->xFwd5 * x_d[9] * x_d[9] - r->xRev5 * x_d[12]; // Checked
    Type dR41 = r->xFwd6 * x_d[13] * x_d[11] / 623 - r->xRev6 * x_d[12]; // Checked
    
    dxdt_d[1] = - dR7 - dR6 - dR5 - dR1 - dR2 + r->expression; // AXL
    dxdt_d[2] = -2*(dR8) - dR5 + dR1 - dR3                   ; // AXLgas1
    dxdt_d[3] = -2*(dR9) - dR6 + dR2 - dR4                   ; // AXLgas2
    dxdt_d[4] = -dR7 + dR3 + dR4                             ; // AXLgas12
    dxdt_d[5] = -dR11 + dR6 + dR5                            ; // AXLdimer1
    dxdt_d[6] = dR11 + dR9 + dR8 + dR7                       ; // AXLdimer2
    
    dxdt_d[7]  = - dR38 - dR37 - dR36 - dR32 - dR33          ; // AXLi
    dxdt_d[8]  = -2*(dR39) - dR36 + dR32 - dR34              ; // AXLgas1i
    dxdt_d[9]  = -2*(dR40) - dR37 + dR33 - dR35              ; // AXLgas2i
    dxdt_d[10] = -dR38 + dR34 + dR35                         ; // AXLgas12i
    dxdt_d[11] = -dR41 + dR37 + dR36                          ; // AXLdimer1i
    dxdt_d[12] = dR41 + dR40 + dR39 + dR38                   ; // AXLdimer2i
    
    dxdt_d[13] = -dR41 - dR32 - dR33 - dR34 - dR35 - r->kDeg*x_d[13];
    
    
    dxdt_d[1] += -x_d[1]*(r->internalize + r->pYinternalize*r->scaleA) + r->kRec*(1-r->fElse)*x_d[7]*internalFrac; // Endocytosis, recycling
    dxdt_d[7] += x_d[1]*(r->internalize + r->pYinternalize*r->scaleA)/internalFrac - r->kRec*(1-r->fElse)*x_d[7] - r->kDeg*r->fElse*x_d[7]; // Endocytosis, recycling, degradation
    
    for (int ii = 2; ii < 6; ii++) {
        dxdt_d[ii]  += -x_d[ii]*(r->internalize + r->pYinternalize*r->scaleA)*endoImpair + r->kRec*(1-r->fElse)*x_d[ii+6]*internalFrac; // Endocytosis, recycling
        dxdt_d[ii+6] += x_d[ii]*(r->internalize + r->pYinternalize*r->scaleA)/internalFrac*endoImpair - r->kRec*(1-r->fElse)*x_d[ii+6] // Endocytosis, recycling
        - r->kDeg*r->fElse*x_d[ii+6]*degImpair; // Degradation
    }
    
    dxdt_d[6]  += -x_d[6]*(r->internalize + r->pYinternalize)*endoImpair + r->kRec*(1-r->fD2)*x_d[12]*internalFrac; // Endocytosis, recycling
    dxdt_d[12] += x_d[6]*(r->internalize + r->pYinternalize)/internalFrac*endoImpair - r->kRec*(1-r->fD2)*x_d[12] - r->kDeg*r->fD2*x_d[12]*degImpair; // Endocytosis, recycling, degradation
}



template <class Type>
inline int AXL_react_diff_new(std::vector<Type> &xx_d, std::vector<Type> &dxxdt_d, void *user_data) {
    size_t pos, spec;
    size_t grid_size = (size_t) xx_d.size()/Nspecies;
    double dRdRMaxRMaxR = maxR*maxR*(1.0/grid_size)*(1.0/grid_size);
    
    for (spec = 0; spec < Nspecies; spec++) {
        if (diffD[spec] == 0) {
            for (size_t ii = 1; ii < (grid_size-1); ii++) dxxdt_d[spec*grid_size + ii] = 0;
            dxxdt_d[spec*grid_size] = 0;
            dxxdt_d[(spec+1)*grid_size - 1] = 0;
        } else {
            for (size_t ii = 1; ii < (grid_size-1); ii++) {
                pos = spec*grid_size + ii;
                dxxdt_d[pos] = (-4.0*xx_d[pos] + (2.0-1.0/ii)*xx_d[pos - 1] + (2.0+1.0/ii)*xx_d[pos + 1])/2/dRdRMaxRMaxR;
            }
            
            // Take care of boundary conditions
            dxxdt_d[spec*grid_size] = 4*(xx_d[spec*grid_size + 1] - xx_d[spec*grid_size])/dRdRMaxRMaxR;
            dxxdt_d[(spec+1)*grid_size - 1] = -4*(xx_d[(spec+1)*grid_size - 1] - xx_d[(spec+1)*grid_size - 2])/dRdRMaxRMaxR;
        }
    }
    
    std::vector<Type> reactIn(Nspecies);
    std::vector<Type> reactOut(Nspecies);
    
    // Add in the reaction for each location
    for (size_t jj = 0; jj < grid_size; jj++) {
        for (size_t ii = 0; ii < Nspecies; ii++) reactIn[ii] = xx_d[ii*grid_size + jj];
        
        AXL_react_new(reactIn,reactOut,user_data);
        
        // Convert by diffusion coefficient and add in reaction
        for (size_t ii = 0; ii < Nspecies; ii++) {
            dxxdt_d[ii*grid_size + jj] *= diffD[ii];
            dxxdt_d[ii*grid_size + jj] += reactOut[ii];
        }
    }
    
    return 0;
}





















#endif /* defined(__UniformOptimization__ReactionCode__) */

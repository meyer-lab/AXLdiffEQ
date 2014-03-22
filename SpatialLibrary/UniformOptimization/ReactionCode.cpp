//
//  ReactionCode.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/11/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include "ReactionCode.h"
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_dense.h>
#include <string>
#include <array>
#include <vector>
#include <cppad/cppad.hpp>
#include <string>
#include <stdexcept>


using namespace std;

double endoImpair = 1; ///< Extent by which to impair endocytosis of Gas6-bound species.
double degImpair = 1;
const double internalFrac = 0.5;
double diffD[Nspecies];
const double fgMgConv = 135.2;


int AXL_react(realtype t, N_Vector xIn, N_Vector dxdtIn, void *user_data) {
    struct rates *r = (struct rates *) user_data;
    realtype* x_d = NV_DATA_S(xIn);
    realtype* dxdt_d = NV_DATA_S(dxdtIn);
    
    // 0 Gas6   // 1 AXL   // 2 A1    // 3 A2
    // 4 A12    // 5 D1    // 6 D2    // 7 AXLi
    // 8 A1i    // 9 A2i   // 10 A12i // 11 D1i    // 12 D2i   // 13 Gasi
    
    dxdt_d[0] = 0;
    
    realtype dR1 = r->Binding1 * x_d[1] * x_d[0] - r->Unbinding1 * x_d[2];
    realtype dR2 = r->Binding2 * x_d[1] * x_d[0] - r->Unbinding2 * x_d[3];
    realtype dR3 = r->Binding2 * x_d[2] * x_d[0] - r->Unbinding2 * x_d[4];
    realtype dR4 = r->Binding1 * x_d[3] * x_d[0] - r->Unbinding1 * x_d[4];
    realtype dR5 = r->xFwd1 * x_d[1] * x_d[2] - r->xRev1 * x_d[5];
    realtype dR6 = r->xFwd2 * x_d[1] * x_d[3] - r->xRev2 * x_d[5];
    realtype dR7 = r->xFwd3 * x_d[1] * x_d[4] - r->xRev3 * x_d[6];
    realtype dR8 = r->xFwd4 * x_d[2] * x_d[2] - r->xRev4 * x_d[6];
    realtype dR9 = r->xFwd5 * x_d[3] * x_d[3] - r->xRev5 * x_d[6];
    realtype dR11 = r->xFwd6 * x_d[0] * x_d[5] - r->xRev6 * x_d[6];
    
    realtype dR32 = r->Binding1 * x_d[7] * x_d[13] / 623 - r->Unbinding1 * x_d[8];
    realtype dR33 = r->Binding2 * x_d[7] * x_d[13] / 623 - r->Unbinding2 * x_d[9];
    realtype dR34 = r->Binding2 * x_d[8] * x_d[13] / 623 - r->Unbinding2 * x_d[10];
    realtype dR35 = r->Binding1 * x_d[9] * x_d[13] / 623 - r->Unbinding1 * x_d[10];
    realtype dR36 = r->xFwd1 * x_d[7] * x_d[8] - r->xRev1 * x_d[11];
    realtype dR37 = r->xFwd2 * x_d[7] * x_d[9] - r->xRev2 * x_d[11];
    realtype dR38 = r->xFwd3 * x_d[7] * x_d[10] - r->xRev3 * x_d[12];
    realtype dR39 = r->xFwd4 * x_d[8] * x_d[8] - r->xRev4 * x_d[12]; // Checked
    realtype dR40 = r->xFwd5 * x_d[9] * x_d[9] - r->xRev5 * x_d[12]; // Checked
    realtype dR41 = r->xFwd6 * x_d[13] * x_d[11] / 623 - r->xRev6 * x_d[12]; // Checked
    
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
    
    return 0;
}


int AXL_react_diff(realtype t, N_Vector xx , N_Vector dxxdt, void *user_data) {
    realtype* xx_d = NV_DATA_S(xx);
    realtype* dxxdt_d = NV_DATA_S(dxxdt);
    size_t pos, spec;
    size_t grid_size = (size_t) NV_LENGTH_S(xx)/Nspecies;
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
    
    N_Vector reactIn = N_VNew_Serial(Nspecies);
    N_Vector reactOut = N_VNew_Serial(Nspecies);
    
    // Add in the reaction for each location
    for (size_t jj = 0; jj < grid_size; jj++) {
        for (size_t ii = 0; ii < Nspecies; ii++) Ith(reactIn,ii) = xx_d[ii*grid_size + jj];
        
        AXL_react(t,reactIn,reactOut,user_data);
        
        // Convert by diffusion coefficient and add in reaction
        for (size_t ii = 0; ii < Nspecies; ii++) {
            dxxdt_d[ii*grid_size + jj] *= diffD[ii];
            dxxdt_d[ii*grid_size + jj] += Ith(reactOut,ii);
        }
    }
    
    N_VDestroy_Serial(reactIn);
    N_VDestroy_Serial(reactOut);
    
    return 0;
}

// This takes the model state and calculates the amount of phosphorylated species
double pYcalc (N_Vector state, struct rates p) {
    double pYa = Ith(state,1) + Ith(state,2) + Ith(state,3) + Ith(state,4) + 2*Ith(state,5) +
    internalFrac*(Ith(state,7) + Ith(state,8) + Ith(state,9) + Ith(state,10) + 2*Ith(state,11));
    
    pYa *= p.scaleA;
    
    pYa += 2*Ith(state,6) + internalFrac*(2*Ith(state,12));
    
    return pYa;
}

// This takes the model state and calculates the total amount of receptor in a cell
double totCalc (N_Vector state) {
    double total = 0;
    
    for (int ii = 1; ii < 7; ii++) total += Ith(state,ii);
    for (int ii = 7; ii < 13; ii++) total += Ith(state,ii)*internalFrac;
    
    total += Ith(state,5);
    total += Ith(state,6);
    total += Ith(state,11)*internalFrac;
    total += Ith(state,12)*internalFrac;
    
    return total/fgMgConv;
}


struct rates Param(param_type params) {
    struct rates out;
    
    for (size_t ii = 0; ii < 15; ii++) {
        if (params[ii] < 0) throw invalid_argument(string("An input model parameter is outside the physical range."));
    }
    
    out.Binding1 = params[0];
    out.Binding2 = params[1];
    out.Unbinding1 = params[2];
    out.Unbinding2 = params[3];
    out.xFwd1 = params[4];
    out.xRev1 = params[5];
    out.xFwd3 = params[6];
    out.xRev3 = params[7];
    out.internalize = params[8];
    out.pYinternalize = params[9];
    out.scaleA = params[10];
    out.kRec = params[11];
    out.kDeg = params[12];
    out.fElse = params[13];
    out.fD2 = params[14];
    out.xRev5 = out.xRev3*out.Unbinding1/out.Unbinding2;
    out.xRev4 = out.xRev3*out.Unbinding2/out.Unbinding1;
    out.xRev2 = out.xRev1*out.Unbinding1/out.Unbinding2;
    out.xFwd2 = out.xFwd1*out.Binding1/out.Binding2;
    out.xFwd4 = out.xFwd3*out.Binding2/out.Binding1;
    out.xFwd5 = out.xFwd3*out.Binding1/out.Binding2;
    out.xFwd6 = out.xFwd3*out.Binding2/out.xFwd1;
    out.xRev6 = out.xRev3*out.Unbinding2/out.xRev1;
    out.JacP = nullptr;
    
    
    return out;
}


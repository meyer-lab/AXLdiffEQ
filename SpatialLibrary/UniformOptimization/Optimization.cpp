//
//  Optimization.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/23/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include "Optimization.h"
#include "ModelRunning.h"
#include "BlasHeader.h"
#include <nlopt.hpp>
#include <random>
#include <vector>

using namespace nlopt;
using namespace std;

void getLimits (vector<double> &minn, vector<double> &maxx, int nCells) {
    minn.clear();
    maxx.clear();
    
    // Ig1 bind
    minn.push_back(0.07918124604);
    maxx.push_back(0.07918124604);
    
    // Ig2 bind
    minn.push_back(-5);
    maxx.push_back(1);
    
    minn.push_back(-1.3767507096);
    maxx.push_back(-1.3767507096);
    
    minn.push_back(-5);
    maxx.push_back(5);
    
    for (size_t ii = 0; ii < 4; ii++) {
        minn.push_back(-5);
        maxx.push_back(5);
    }
    
    minn.push_back(-6);
    maxx.push_back(0);
    
    minn.push_back(-3);
    maxx.push_back(2);
    
    minn.push_back(-6);
    maxx.push_back(-1);
    
    minn.push_back(-3);
    maxx.push_back(-1);
    
    minn.push_back(-4);
    maxx.push_back(-1);
    
    minn.push_back(-5);
    maxx.push_back(0);
    
    minn.push_back(-1);
    maxx.push_back(0);
    
    // Gas6 autocrine
    for (size_t ii = 0; ii < nCells; ii++) {
        minn.push_back(-5);
        maxx.push_back(2);
    }
    
    // AXL expression
    for (size_t ii = 0; ii < nCells; ii++) {
        minn.push_back(0);
        maxx.push_back(5);
    }
}

double calcErrorOptLog (unsigned n, const double *x, double *grad, void *data) {
    param_type inner;
    for (size_t i = 0; i < n; i++) inner[i] = pow(10,x[i]);
    
    return calcError(inner);
}

double calcErrorOptOneLog (unsigned n, const double *x, double *grad, void *data) {
    int *line = (int *) data;
    
    double xIn[n];
    
    for (size_t ii = 0; ii < n; ii++) xIn[ii] = pow(10,x[ii]);
    
    return calcErrorOneCellLine (*line, xIn);
}

void bumpOptim(vector<double> minn, vector<double> maxx, vector<double> xx, double *ff, double strength,
                 unsigned int seed, nlopt_func minFun, void *data) {
    
    double outVal;
    default_random_engine generator;
    generator.seed(seed);
    normal_distribution<double> normRnd(0,strength);
    uniform_real_distribution<double> uniRnd(0,1);
    
    for (int ii = 0; ii < xx.size(); ii++) {
        xx[ii] = (xx[ii] + normRnd(generator)*(maxx[ii] - minn[ii]));
        
        if ((xx[ii] < minn[ii]) || (xx[ii] > maxx[ii]))
            xx[ii] = minn[ii] + (maxx[ii] - minn[ii]) * uniRnd(generator);
    }
    
    nlopt::opt opter = nlopt::opt(nlopt::algorithm::LN_COBYLA, (unsigned int) xx.size());
    opter.set_lower_bounds(minn);
    opter.set_upper_bounds(maxx);
    opter.set_xtol_rel(1E-6);
    opter.set_ftol_rel(1E-5);
    opter.set_min_objective(minFun, data);
    opter.optimize(xx, outVal);
    
    *ff = outVal;
}



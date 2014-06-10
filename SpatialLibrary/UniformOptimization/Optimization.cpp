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
#include "ReactionCode.h"
#include <random>
#include <vector>
#include <algorithm>
#include <mutex>

using namespace nlopt;
using namespace std;


vector<double> xBest;
double fbest = 1e10;
mutex mtx;

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}

void globalBest (size_t n, const double *x, double f) {
    mtx.lock();
    
    if (f < fbest) {
        fbest = f;
        
        xBest.clear();
        
        cout << "Best! " << f << endl;
        
        for (size_t i = 0; i < n; i++) {
            xBest.push_back(x[i]);
            cout << x[i] << " ";
        }
        
        cout << endl;
    }
    
    mtx.unlock();
}


void getLimits (vector<double> &minn, vector<double> &maxx, int nCells) {
    minn.clear();
    maxx.clear();
    
    // Ig1 bind
    minn.push_back(-4);
    maxx.push_back(1);
    
    // Ig2 bind
    minn.push_back(-4);
    maxx.push_back(1);
    
    minn.push_back(-5);
    maxx.push_back(5);
    
    minn.push_back(-5);
    maxx.push_back(5);
    
    for (size_t ii = 0; ii < 4; ii++) {
        minn.push_back(-5);
        maxx.push_back(5);
    }
    
    minn.push_back(-6);
    maxx.push_back(0);
    
    minn.push_back(-3);
    maxx.push_back(1);
    
    minn.push_back(-6);
    maxx.push_back(0);
    
    minn.push_back(-3);
    maxx.push_back(-1);
    
    minn.push_back(-4);
    maxx.push_back(-1);
    
    minn.push_back(-5);
    maxx.push_back(0);
    
    minn.push_back(-1);
    maxx.push_back(0);
    
    // Gas6 autocrine
    for (size_t ii = 0; ii < (size_t) nCells; ii++) {
        minn.push_back(-6);
        maxx.push_back(0);
    }
    
    // AXL expression
    for (size_t ii = 0; ii < (size_t) nCells; ii++) {
        minn.push_back(0);
        maxx.push_back(5);
    }
}

void getRandLimits (vector<double> &minn, vector<double> &maxx, int nCells) {
    minn.clear();
    maxx.clear();
    
    // Ig1 bind
    minn.push_back(-2.3829);
    maxx.push_back(-2.3829);
    
    // Ig2 bind
    minn.push_back(0.77811	);
    maxx.push_back(0.77811	);
    
    minn.push_back(-5);
    maxx.push_back(0);
    
    minn.push_back(-5);
    maxx.push_back(5);
    
    // Receptor reactions
    minn.push_back(-5);
    maxx.push_back(0);
    
    minn.push_back(-5);
    maxx.push_back(5);
    
    minn.push_back(-5);
    maxx.push_back(0);
    
    minn.push_back(-5);
    maxx.push_back(5);
    
    // AXLint1
    minn.push_back(-2.2774);
    maxx.push_back(-2.2774);
    
    // AXLint2
    minn.push_back(-3);
    maxx.push_back(-3);
    
    // ScaleA
    minn.push_back(-8);
    maxx.push_back(0);
    
    // kRec
    minn.push_back(-2.5648);
    maxx.push_back(-2.5648);
    
    // kDeg
    minn.push_back(-1);
    maxx.push_back(-1);
    
    // fElse
    minn.push_back(-1.2476);
    maxx.push_back(-1.2476);
    
    // fD2
    minn.push_back(0);
    maxx.push_back(0);
    
    // Gas6 autocrine
    minn.push_back(-2);
    maxx.push_back(-2);
    
    // AXL expression
    minn.push_back(2);
    maxx.push_back(2);
}

double calcErrorOptLog (unsigned n, const double *x, double *grad, void *data) {
    param_type inner;
    for (size_t i = 0; i < n; i++) inner[i] = pow(10,x[i]);
    
    double error = calcError(inner);
    
    globalBest(n, x, error);
    
    return error;
}

double calcErrorOptOneLog (unsigned n, const double *x, double *grad, void *data) {
    int *line = (int *) data;
    
    double xIn[n];
    
    for (size_t ii = 0; ii < n; ii++) xIn[ii] = pow(10,x[ii]);
    
    double error = calcErrorOneCellLine (*line, xIn);
    
    globalBest(n, x, error);
    
    return error;
}

double calcErrorOptA549Full (unsigned n, const double *x, double *grad, void *data) {
    
    param_type xIn;
    
    for (size_t ii = 0; ii < n; ii++) xIn[ii] = pow(10,x[ii]);
    
    double error = calcErrorA549Full (Param(xIn), xIn[15]);
    
    globalBest(n, x, error);
    
    return error;
}


double calcErrorOptAllLog (unsigned n, const double *x, double *grad, void *data) {
    param_type xIn;
    const int nCellLines = 4;

    double autocrine[nCellLines];
    double expression[nCellLines];

    for (size_t ii = 0; ii < nCellLines; ii++) {
    	autocrine[ii] = pow(10,x[ii+15]);
    	expression[ii] = pow(10,x[ii+15+nCellLines]);
    }

    for (size_t ii = 0; ii < NELEMS(xIn); ii++) xIn[ii] = pow(10,x[ii]);

    double error = calcErrorAll(Param(xIn), expression, autocrine);
    
    globalBest(n, x, error);
    
    return error;
}

double calcErrorOptPaperSiLog (unsigned n, const double *x, double *grad, void *data) {
    param_type xIn;

    for (size_t ii = 0; ii < xIn.size(); ii++) xIn[ii] = pow(10,x[ii]);

    double error = calcError(xIn) + calcErrorSi(xIn);
    
    globalBest(n, x, error);
    
    return error;
}

void bumpOptimGlobal(vector<double> minn, vector<double> maxx, nlopt_func minFun, void *data, int method) {
    using nlopt::algorithm;
    
    vector<double> xx;
    double outVal;
    nlopt::opt opterG;
    default_random_engine generator;
    generator.seed(rdtsc());
    uniform_real_distribution<double> uniRnd(0,1);
    uniform_int_distribution<int> intRand(0,3);
    
    nlopt_srand(rdtsc());
    
    for (size_t i = 0; i < minn.size(); i++) {
        xx.push_back(minn[i] + (maxx[i] - minn[i]) * uniRnd(generator));
    }
    
    if (method == -1) method = intRand(generator);
    
    if (method == 0) opterG = nlopt::opt(G_MLSL_LDS, (unsigned int) xx.size());
    if (method == 1) opterG = nlopt::opt(GN_ISRES, (unsigned int) xx.size());
    if (method == 2) opterG = nlopt::opt(GN_MLSL_LDS, (unsigned int) xx.size());
    if (method == 3) opterG = nlopt::opt(GN_CRS2_LM, (unsigned int) xx.size());
    //if (method == 4) opterG = nlopt::opt(GN_DIRECT_L_RAND, (unsigned int) xx.size());
    
    opterG.set_lower_bounds(minn);
    opterG.set_upper_bounds(maxx);
    opterG.set_min_objective(minFun, data);
    
    if (method == 0) {
        nlopt::opt opter = nlopt::opt(nlopt::algorithm::LN_BOBYQA, (unsigned int) xx.size());
        opter.set_lower_bounds(minn);
        opter.set_upper_bounds(maxx);
        opter.set_xtol_rel(1E-6);
        opter.set_ftol_rel(1E-5);
        opter.set_min_objective(minFun, data);
        opterG.set_local_optimizer(opter);
    }
    
    //opterG.set_maxtime(30);
    opterG.optimize(xx, outVal);
}

void randLargeResponse(vector<double> minn, vector<double> maxx) {
    double xx[18];
    double dataPtr[2];
    double tps[2] = {0, 0.5};
    const double threshold = 20;
    
    vector<double> foundMin = maxx;
    vector<double> foundMax = minn;
    
    size_t N = 0;
    
    default_random_engine generator;
    generator.seed(rdtsc());
    uniform_real_distribution<double> uniRnd(0,1);
    nlopt_srand(rdtsc());
    
    while (1) {
        N++;
        
        for (size_t i = 0; i < minn.size(); i++) {
            xx[i] = pow(10,(minn[i] + (maxx[i] - minn[i]) * uniRnd(generator)));
        }
        
        if ((xx[3] - xx[1]) < (xx[2] - xx[0])) continue;
        if ((xx[2] - xx[0]) > 0) continue;
        if ((xx[3] - xx[1]) < 1) continue;
        if ((xx[3] - xx[1]) > 3) continue;
        
        
        if (calcProfileMatlab(dataPtr, xx, tps, 2, xx[15], xx[16], 1.25, 5)) continue;

        if (dataPtr[1] < 0.3) continue;
        if ((dataPtr[1] / dataPtr[0]) < threshold) continue;
        
        
        for (size_t i = 0; i < minn.size(); i++) {
            if (log10(xx[i]) < foundMin[i]) foundMin[i] = log10(xx[i]);
            if (log10(xx[i]) > foundMax[i]) foundMax[i] = log10(xx[i]);
        }
        
        mtx.lock();
        
        for (size_t i = 0; i < minn.size(); i++) {
            cout << foundMin[i] << " ";
        }
        //cout << endl;
//        for (size_t i = 0; i < minn.size(); i++) {
//            cout << foundMax[i] << " ";
//        }
        
        cout << N << endl;
        
        mtx.unlock();
    }
}

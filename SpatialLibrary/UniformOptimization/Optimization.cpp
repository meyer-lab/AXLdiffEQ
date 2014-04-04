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


#define noiseOut 0

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
    for (size_t ii = 0; ii < (size_t) nCells; ii++) {
        minn.push_back(-5);
        maxx.push_back(2);
    }
    
    // AXL expression
    for (size_t ii = 0; ii < (size_t) nCells; ii++) {
        minn.push_back(0);
        maxx.push_back(5);
    }
}

void getLimits_sepA (vector<double> &minn, vector<double> &maxx, int nCells) {
    minn.clear();
    maxx.clear();
    
    // Ig1 bind
    minn.push_back(0.07918124604);
    maxx.push_back(0.07918124604);
    
    // Ig2 bind
    minn.push_back(-5);
    maxx.push_back(0.08);
    
    minn.push_back(-1.3767507096);
    maxx.push_back(-1.3767507096);
    
    minn.push_back(-1);
    maxx.push_back(2);
    
    // Receptor dimer parameters
    for (size_t ii = 0; ii < 4; ii++) {
        minn.push_back(-5);
        maxx.push_back(2);
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
    
    minn.push_back(-5);
    maxx.push_back(0);
    
    // Gas6 autocrine
    for (size_t ii = 0; ii < (size_t) nCells; ii++) {
        minn.push_back(-5);
        maxx.push_back(2);
    }
    
    // AXL expression
    for (size_t ii = 0; ii < (size_t) nCells; ii++) {
        minn.push_back(0);
        maxx.push_back(5);
    }
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



double calcErrorOptPaperSiOneLog_sepA (unsigned n, const double *x, double *grad, void *data) {
    int *line = (int *) data;
    
    param_type xIn;
    
    for (size_t ii = 0; ii < n; ii++) xIn[ii] = pow(10,x[ii]);
    
    struct rates_sepA rates = Param_sepA(xIn);
    
    rates.expression = xIn[n-1];
    
    double error = calcErrorOneLine_sepA (rates, *line, xIn[n-2]) + calcErrorSiOneLine_sepA (rates, *line, xIn[n-2]);
    
    //cout << error << endl;
    
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

double calcErrorOptPaperSiLog_sepA (unsigned n, const double *x, double *grad, void *data) {
    vector<double> xIn(n);
    
    for (size_t ii = 0; ii < n; ii++) xIn[ii] = pow(10,x[ii]);
    
    double error = calcError_sepA(xIn) + calcErrorSi_sepA(xIn);
    
    if (noiseOut) cout << error << endl;
    
    globalBest(n, x, error);
    
    return error;
}

double calcErrorSiLog_sepA (unsigned n, const double *x, double *grad, void *data) {
    vector<double> xIn;
    
    for (size_t ii = 0; ii < n; ii++) xIn.push_back(pow(10,x[ii]));
    
    double error = calcErrorSi_sepA (xIn);
    
    globalBest(n, x, error);
    
    return error;
}

double calcErrorOptAllSiLog_sepA (unsigned n, const double *x, double *grad, void *data) {
    param_type xIn;
    const size_t nCellLines = (n-15)/2;
    
    
    vector<double> siIn;
    for (size_t ii = 0; ii < n; ii++) siIn.push_back(pow(10,x[ii]));
    
    
    double autocrine[nCellLines];
    double expression[nCellLines];
    
    for (size_t ii = 0; ii < nCellLines; ii++) {
    	autocrine[ii] = pow(10,x[ii+16]);
    	expression[ii] = pow(10,x[ii+16+nCellLines]);
    }
    
    for (size_t ii = 0; ii < NELEMS(xIn); ii++) xIn[ii] = pow(10,x[ii]);
    
    double error = calcErrorAll_sepA(Param_sepA(xIn), expression, autocrine) + calcErrorSi_sepA (siIn);
    
    globalBest(n, x, error);
    
    return error;
}

double calcErrorOptPaperSiAllLog_sepA (unsigned n, const double *x, double *grad, void *data) {
    clock_t bbegin, endd;
    
    if (noiseOut) bbegin = clock();
    
    
    vector<double> xIn;
    vector<double> siIn;
    for (size_t ii = 0; ii < n; ii++) siIn.push_back(pow(10,x[ii]));
    
    xIn = siIn;
    
    xIn.erase(xIn.end()-3, xIn.end());
    xIn.erase(xIn.end()-5, xIn.end()-2);
    
    double error = calcError_sepA(xIn) + calcErrorSi_sepA (siIn);
    
    if (noiseOut) {
        endd = clock();
    
        cout << "Time: " << ((double) (endd - bbegin))/((double) CLOCKS_PER_SEC) << endl;
        cout << "Err: " << error << endl << endl;
    }
    
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
    
    nlopt_srand(rdtsc());
    
    for (size_t i = 0; i < minn.size(); i++) {
        xx.push_back(minn[i] + (maxx[i] - minn[i]) * uniRnd(generator));
    }
    
    if (method == 0) opterG = nlopt::opt(G_MLSL_LDS, (unsigned int) xx.size());
    if (method == 1) opterG = nlopt::opt(GN_ISRES, (unsigned int) xx.size());
    if (method == 2) opterG = nlopt::opt(GN_MLSL_LDS, (unsigned int) xx.size());
    if (method == 3) opterG = nlopt::opt(GN_CRS2_LM, (unsigned int) xx.size());
    if (method == 4) opterG = nlopt::opt(GN_DIRECT_L_RAND, (unsigned int) xx.size());
    
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




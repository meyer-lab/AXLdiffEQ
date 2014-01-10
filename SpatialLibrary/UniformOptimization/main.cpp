/**
 *   \file main.cpp
 *   FiniteDiffGSL
 *
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */
#include <nlopt.hpp>
#include <random>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <fstream>

#include "BlasHeader.h"
#include "Fitting.h"

using namespace nlopt;
using namespace std;

mutex mtx;
atomic<double> fBest;

double calcErrorOpt (unsigned n, const double *x, double *grad, void *data) {
    param_type inner;
    for (size_t ii = 0; ii < inner.size(); ii++) inner[ii] = x[ii];
    
    return calcError(inner);
}

void printFitState(param_type in) {
    cout << endl << endl;
}


double calcErrorOptLog (unsigned n, const double *x, double *grad, void *data) {
    ofstream errfile;
    param_type inner;
    for (size_t ii = 0; ii < inner.size(); ii++) inner[ii] = pow(10,x[ii]);
    
    double error = calcError(inner);

    if (error < fBest) {
        fBest = error;
        
        mtx.lock();
        cout << fBest << endl;
        
        printFitState(inner);
        
        //errfile.open ("output.log", ios::app);
        //errfile << fBest << endl;
        //errfile.close();
        mtx.unlock();
    }
    
    return error;
}

double calcErrorOpt (vector<double> x) {
    param_type inner;
    for (size_t ii = 0; ii < x.size(); ii++) inner[ii] = x[ii];
    return calcError(inner);
}

double calcErrorOptLog (vector<double> x) {
    param_type inner;
    for (size_t ii = 0; ii < x.size(); ii++) inner[ii] = pow(10,x[ii]);
    return calcError(inner);
}

void calcOptimizer(double seed) {
    double ff = 0;
    fBest = 1E10;
    vector<double> xx;
    
    vector<double> minn = {-1,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-3,-3, -5, -5, 0, 0};
    vector<double> maxx = {1 ,5 ,1 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 , 5, 5, 0, 0, 1, 1, 5, 5};
    
    cout << minn.size() << endl;
    cout << maxx.size() << endl;
    
    default_random_engine generator;
    generator.seed(seed);
    uniform_real_distribution<double> uniRnd(0,1);
    
    for (int ii = 0; ii < minn.size(); ii++) {
        xx.push_back(minn[ii] + (maxx[ii] - minn[ii]) * uniRnd(generator));
    }
    
    opt opter = opt(algorithm::GN_CRS2_LM, (unsigned int) minn.size());
    opter.set_lower_bounds(minn);
    opter.set_upper_bounds(maxx);
    opter.set_population(1000);
    opter.set_min_objective(calcErrorOptLog, nullptr);
    
    int flag = opter.optimize(xx, ff);
    
    mtx.lock();
    cout << flag << ", " << ff << endl;
    mtx.unlock();
}

int main()
{
//    std::thread threads[8];
//    
//    for (int i=0; i < NELEMS(threads); ++i)
//        threads[i] = std::thread(calcOptimizer,random());
//    
//    for (auto& th : threads) th.join();
    
    param_type in = {0.6,
        2.2822811E-05,
        1.0904953E-05,
        1.0436166E-03, // Binding params
        
        1.6057974E-02,
        1.4127670E+00,
        3.0836805E+03,
        1.6612537E-02, // xFwd Rev
        
        1.2994370E+00, // AXLint1
        6.0380584E+00, // AXLint2
        
        1.2561171E-01, // scaleA
        
        1.3818834E-01,
        1.2482924E-01,
        2.8639311E-01,
        4.4657377E-01,
        
        2.1667962E-02,
        1.9510010E-01,
        7.7814510E+03,
        1.1154666E+04};
    
    cout << calcError(in) << endl;
    
    return 0;
}

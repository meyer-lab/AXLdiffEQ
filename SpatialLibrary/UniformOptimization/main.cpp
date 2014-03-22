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
#include <atomic>
#include <fstream>
#include "BlasHeader.h"
#include "ReactionCode.h"
#include "ModelRunning.h"


//using namespace nlopt;
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

void calcOptimizer(unsigned int seed) {
    double ff = 0;
    fBest = 1E10;
    vector<double> xx;
    
    vector<double> minn = {-1,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-3,-3, -5, -5, 0, 0};
    vector<double> maxx = {1 ,5 ,1 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 , 5, 5, 0, 0, 1, 1, 5, 5};
    
    default_random_engine generator;
    generator.seed(seed);
    uniform_real_distribution<double> uniRnd(0,1);
    
    for (int ii = 0; ii < minn.size(); ii++) {
        xx.push_back(minn[ii] + (maxx[ii] - minn[ii]) * uniRnd(generator));
    }
    
    nlopt::opt opter = nlopt::opt(nlopt::algorithm::GD_STOGO_RAND, (unsigned int) minn.size());
    opter.set_lower_bounds(minn);
    opter.set_upper_bounds(maxx);
    opter.set_population(1000);
    opter.set_maxtime(60);
    opter.set_min_objective(calcErrorOptLog, nullptr);
    
    int flag = opter.optimize(xx, ff);
    
    mtx.lock();
    cout << flag << ", " << ff << endl;
    mtx.unlock();
}


void testDiffModel () {
    const size_t gridSize = 100;
    
    double tps[] = {30};
    double data[1];
    double GasIn[gridSize];
    
    double params[] = {1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571, 0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1, 0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46};
    
    for (int ii = 0; ii < gridSize; ii++) {
        GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
    }
    
    double dIn[] = {0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    cout << matlabDiffTPS_pYavg(data, 359.46, GasIn, gridSize, 0.001, params, tps, 1, dIn, 1, 1, 0) << endl;
}



int main()
{
//    std::thread threads[1];
//    
//    for (int i=0; i < NELEMS(threads); ++i)
//        threads[i] = std::thread(calcOptimizer,random());
//    
//    for (auto& th : threads) th.join();
//
//    param_type in = {0.6, 2.2822811E-05, 1.0904953E-05, 1.0436166E-03, // Binding params
//        1.6057974E-02,
//        1.4127670E+00,
//        3.0836805E+03,
//        1.6612537E-02, // xFwd Rev
//        1.2994370E+00, // AXLint1
//        6.0380584E+00, // AXLint2
//        1.2561171E-01, // scaleA
//        1.3818834E-01,
//        1.2482924E-01,
//        2.8639311E-01,
//        4.4657377E-01,
//        2.1667962E-02,
//        1.9510010E-01,
//        7.7814510E+03,
//        1.1154666E+04};
//    
//    for (int ii = 0; ii < 100; ii++) calcError(in);
//    cout << calcError(in) << endl;
    
    for (int ii = 1; ii <= 100; ii++) {
        cout << "Test " << ii << " of 100." << endl;
        testDiffModel();
    }
    
//    double in[] = {0.6,
//        2.2822811E-05,
//        1.0904953E-05,
//        1.0436166E-03, // Binding params
//        
//        1.6057974E-02,
//        1.4127670E+00,
//        3.0836805E+03,
//        1.6612537E-02, // xFwd Rev
//        
//        1.2994370E+00, // AXLint1
//        6.0380584E+00, // AXLint2
//        
//        1.2561171E-01, // scaleA
//        
//        1.3818834E-01,
//        1.2482924E-01,
//        2.8639311E-01,
//        4.4657377E-01,
//        7.7814510E+03,
//        2.1667962E-02};
//    
//    
//    for (int ii = 0; ii < 1; ii++) cout << calcErrorOneCellLine (0, in) << endl;
    
    //calcOptimizer(0);
    
    
    return 0;
}

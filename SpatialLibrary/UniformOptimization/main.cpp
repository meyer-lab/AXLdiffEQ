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
#include "Optimization.h"
#include "BlasHeader.h"
#include "ReactionCode.h"
#include "ModelRunning.h"

using namespace std;

void testDiffModel () {
    const size_t gridSize = 100;
    
    double tps[] = {30};
    double data[1];
    double GasIn[gridSize];
    
    double params[] = {1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571, 0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1, 0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46};
    
    for (int ii = 0; ii < gridSize; ii++) GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
    
    double dIn[] = {0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    cout << matlabDiffTPS_pYavg(data, 359.46, GasIn, gridSize, 0.001, params, tps, 1, dIn, 1, 1, 0) << endl;
}

int main() {
//    std::thread threads[4];
//    
//    for (int i=0; i < NELEMS(threads); ++i)
//        threads[i] = std::thread(calcOptimizer,time(NULL) + i);
//    
//    for (auto& th : threads) th.join();
    
    vector<double> minn;
    vector<double> maxx;
    vector<double> in;
    
    getLimits(minn, maxx, 2);
    
    default_random_engine generator;
    generator.seed(time(NULL));
    uniform_real_distribution<double> uniRnd(0,1);
    
    for (int ii = 0; ii < minn.size(); ii++) {
        in.push_back((minn[ii] + (maxx[ii] - minn[ii]) * uniRnd(generator)));
    }
    
    double fBest = 1e6;
    double out;
    int cellLine = 1;
    
    vector<double> best(in);
    
    while (fBest == 1e6) {
        out = bumpOptim(minn, maxx, in, 1, time(NULL), calcErrorOptLog, &cellLine);
        
        if (out < fBest) {
            fBest = out;
            best = in;
        } else in = best;
    }
    
    while (1) {
        for (double strength = 1; strength > 0.001; strength /= 5) {
            out = bumpOptim(minn, maxx, in, strength, time(NULL), calcErrorOptLog, &cellLine);
            
            if (out < fBest) {
                fBest = out;
                best = in;
                cout << "Best! " << fBest << endl;
                
                for(auto& i : best) cout << i << ' ';
                
                cout << endl;
                
            } else in = best;
        }
    }
    
    return 0;
}

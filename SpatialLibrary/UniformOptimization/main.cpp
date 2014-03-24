/**
 *   \file main.cpp
 *   FiniteDiffGSL
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */
#include <random>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include "Optimization.h"

#define powRND pow(10,-4.0*uniRnd(generator))

using namespace std;

int main() {
    const int nthreads = 4;
    chrono::milliseconds dura( 200 );
    vector<double> minn, maxx, best;
    vector<double> inn[nthreads];
    thread threads[nthreads];
    size_t ii;
    double fBest = 1e6;
    double out[nthreads];
    int cellLine[] = {1};
    
    getLimits(minn, maxx, 1);
    
    default_random_engine generator;
    generator.seed(time(NULL));
    uniform_real_distribution<double> uniRnd(0,1);
    
    for (ii = 0; ii < minn.size(); ii++) best.push_back(minn[ii] + (maxx[ii] - minn[ii]) * uniRnd(generator));
    
    for (ii = 0; ii < nthreads; ii++) {
        inn[ii] = best;
        out[ii] = -1;
        
        threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,time(NULL),calcErrorOptOneLog,cellLine);
    }
    
    while (true) {
        for (ii = 0; ii < nthreads; ii++) {
            if (out[ii] == -1) continue;
            
            if (out[ii] < fBest) {
                fBest = out[ii];
                best = inn[ii];
                cout << "Best! " << fBest << endl;
                
                for(auto& i : best) cout << i << ' ';
                
                cout << endl;
            }
            
            inn[ii] = best;
            out[ii] = -1;
            
            threads[ii].join();
            
            threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,time(NULL),calcErrorOptOneLog,cellLine);
        }
        this_thread::sleep_for(dura);
    }
    
    return 0;
}

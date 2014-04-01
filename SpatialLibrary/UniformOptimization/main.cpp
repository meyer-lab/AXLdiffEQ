/**
 *   \file main.cpp
 *   FiniteDiffGSL
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */
#include <random>
#include <iostream>
#include <vector>
#include <boost/thread/thread.hpp>
#include <climits>
#include <chrono>
#include <sstream>
#include "Optimization.h"
#include "BlasHeader.h"

#define powRND pow(10,-4.0*uniRnd(generator))
#define seedRND intRnd(generator)

using namespace std;
using namespace boost;

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}


int main( int argc, char *argv[] ) {

	if (argc < 2) {
		cout << "You need to tell me the cell line you want." << endl;

		return 1;
	}

	istringstream ss(argv[1]);
	int x;
	if (!(ss >> x) || (x > 3) || (x < -6)) {
	    cerr << "Invalid cell line number " << argv[1] << '\n';
	    return 1;
	}

	cout << "Running optimization for cell line " << x << endl;

    const int nthreads = 8;
    boost::chrono::milliseconds dura( 200 );
    vector<double> minn, maxx, best;
    vector<double> inn[nthreads];
    thread threads[nthreads];
    size_t ii;
    double fBest = 1e6;
    double out[nthreads];
    int cellLine[] = {x};
    
    if (cellLine[0] == -1) getLimits(minn, maxx, 4);
    else if (cellLine[0] == -2) getLimits(minn, maxx, 2);
    else if (cellLine[0] == -3) getLimits_sepA(minn, maxx, 5);
    else if (cellLine[0] == -4) getLimits_sepA(minn, maxx, 2);
    else if (cellLine[0] == -5) getLimits_sepA(minn, maxx, 5);
    else if (cellLine[0] == -6) getLimits_sepA(minn, maxx, 5);
    else getLimits(minn, maxx, 1);
    
    default_random_engine generator;

    generator.seed(rdtsc());
    uniform_real_distribution<double> uniRnd(0,1);
    uniform_int_distribution<unsigned int> intRnd(0, UINT_MAX);
    
    for (ii = 0; ii < minn.size(); ii++) best.push_back(minn[ii] + (maxx[ii] - minn[ii]) * uniRnd(generator));
    
    for (ii = 0; ii < nthreads; ii++) {
        inn[ii] = best;
        out[ii] = -1;
        
        if (cellLine[0] == -1) {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptAllLog,nullptr);
        } else if (cellLine[0] == -2) {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptPaperSiLog,nullptr);
        } else if (cellLine[0] == -3) {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptAllSiLog_sepA,nullptr);
        } else if (cellLine[0] == -4) {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptPaperSiLog_sepA,nullptr);
        } else if (cellLine[0] == -5) {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorSiLog_sepA,nullptr);
        } else if (cellLine[0] == -6) {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptPaperSiAllLog_sepA,nullptr);
        } else {
        	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptOneLog,cellLine);
        }
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
            
            if (cellLine[0] == -1) {
            	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptAllLog,nullptr);
            } else if (cellLine[0] == -2) {
            	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptPaperSiLog,nullptr);
            } else if (cellLine[0] == -3) {
                threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptAllSiLog_sepA,nullptr);
            } else if (cellLine[0] == -4) {
                threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptPaperSiLog_sepA,nullptr);
            } else if (cellLine[0] == -5) {
                threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorSiLog_sepA,nullptr);
            } else if (cellLine[0] == -6) {
                threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptPaperSiAllLog_sepA,nullptr);
            } else {
            	threads[ii] = thread(bumpOptim,minn,maxx,inn[ii],&out[ii],powRND,seedRND,calcErrorOptOneLog,cellLine);
            }
        }
        this_thread::sleep_for(dura);
    }
    
    return 0;
}

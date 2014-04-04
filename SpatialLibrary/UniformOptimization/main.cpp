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
#include <sstream>
#include "Optimization.h"
#include "BlasHeader.h"

using namespace std;
using namespace boost;

int main( int argc, char *argv[] ) {
   
	if (argc < 2) {
		cout << "You need to tell me the cell line you want." << endl;

		return 1;
	}

	istringstream ss(argv[1]);
	int x;
	if (!(ss >> x) || (x > 3) || (x < -7)) {
	    cerr << "Invalid cell line number " << argv[1] << '\n';
	    return 1;
	}

	cout << "Running optimization for cell line " << x << endl;

    const int nthreads = 4;
    vector<double> minn, maxx, best;
    vector<double> inn[nthreads];
    thread threads[nthreads];
    size_t ii;
    double out[nthreads];
    int cellLine[] = {x};
    
    int Acell[] = {1};
    
    if (cellLine[0] == -1) getLimits(minn, maxx, 4);
    else if (cellLine[0] == -2) getLimits(minn, maxx, 2);
    else if (cellLine[0] == -3) getLimits_sepA(minn, maxx, 5);
    else if (cellLine[0] == -4) getLimits_sepA(minn, maxx, 2);
    else if (cellLine[0] == -5) getLimits_sepA(minn, maxx, 5);
    else if (cellLine[0] == -6) getLimits_sepA(minn, maxx, 5);
    else if (cellLine[0] == -7) getLimits_sepA(minn, maxx, 1);
    else getLimits(minn, maxx, 1);
    
    int method = 4;
    
    for (ii = 0; ii < nthreads; ii++) {
        inn[ii] = best;
        out[ii] = -1;
        
        if (cellLine[0] == -1) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptAllLog,nullptr, method);
        } else if (cellLine[0] == -2) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptPaperSiLog,nullptr, method);
        } else if (cellLine[0] == -3) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptAllSiLog_sepA,nullptr, method);
        } else if (cellLine[0] == -4) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptPaperSiLog_sepA,nullptr, method);
        } else if (cellLine[0] == -5) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorSiLog_sepA,nullptr, method);
        } else if (cellLine[0] == -6) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptPaperSiAllLog_sepA,nullptr, method);
        } else if (cellLine[0] == -7) {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptPaperSiOneLog_sepA,Acell, method);
        } else {
        	threads[ii] = thread(bumpOptimGlobal,minn,maxx,calcErrorOptOneLog,cellLine, method);
        }
    }
    
    for (ii = 0; ii < nthreads; ii++) threads[ii].join();
    
    return 0;
}

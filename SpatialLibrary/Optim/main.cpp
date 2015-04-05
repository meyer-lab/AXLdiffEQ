//
//  main.cpp
//  Optim
//
//  Created by Aaron Meyer on 1/14/15.
//  Copyright (c) 2015 Aaron Meyer. All rights reserved.
//

#include <nlopt.hpp>
#include <iostream>
#include "BlasHeader.h"
#include <cmath>
#include <random>

double optimFunc (unsigned, const double *, double *, void *);

using namespace std;


double optimFunc (const unsigned n, const double *x, double *, void *) {
    double in[n];
    
    for (int ii = 0; ii < (n-1); ii++) {
        in[ii] = pow(10,x[ii]);
    }
    
    if (x[10] < 0.5) {
        in[10] = 0;
    } else {
        in[10] = 1;
    }
    
    double out = pyEntry(in);
    
    return out;
}

static void runOptim() {
    default_random_engine gen((unsigned int) chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> randN(0, 1);
    
    const int len = 11;
    double minn[] = {0.6 ,1E-15,1E-5,1E-3,1E-3,1E-4,1E-4,1E-2,100,1E-3, 1};
    double maxx[] = {600 ,1E2  , 600,   1,   1, 0.1,   1,   1,1E5,   1, 10};
    double init[len];
    double out;
    
    for (int ii = 0; ii < len; ii++) {
        minn[ii] = log10(minn[ii]);
        maxx[ii] = log10(maxx[ii]);
    }
    
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, len); /* algorithm and dim */
    nlopt_set_lower_bounds(opt, minn);
    nlopt_set_upper_bounds(opt, maxx);
    nlopt_set_min_objective(opt, optimFunc, NULL);
    
    nlopt_result flag;
    
    for (int jj = 0; jj < 1E6; jj++) {
        for (int ii = 0; ii < len; ii++) {
            init[ii] = minn[ii] + (maxx[ii] - minn[ii])*randN(gen);
        }
        
        flag = nlopt_optimize(opt, init, &out);
        
        
        if (out < 100) {
            cout << jj << "   ";
            
            for (int xx = 0; xx < len; xx++)
                cout << init[xx] << " ";
            
            cout << "   " << out << endl;
        }
    }
    
    nlopt_destroy(opt);
}

int main() {
    
    
    double params[] = {27.003, 5.0203e-06, 4.4569, 0.0018128, 4.467, 0.07214, 0.0010026, 0.0017783, 0.19216, 0.13615, 10190, 0.038456, 1};
    
    double GasIn[100];
    double tps[1] = {10};
    double dIn[20];
    double dataVoid[400];
    
    dIn[0] = 10;
    for (size_t ii = 1; ii < 20; ii++) dIn[ii] = 0;
    
    for (size_t ii = 0; ii < 100; ii++) {
        GasIn[ii] = 10 / ((double) ii + 1.0);
    }
    
    
    diffCalc(dataVoid, dataVoid, dataVoid, dataVoid, GasIn, 100, params, tps, 1, dIn);
    
    
    
    return 0;
}

//
//  Optimization.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/23/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__Optimization__
#define __UniformOptimization__Optimization__

#include <iostream>
#include <vector>

typedef double (*nlopt_func)(unsigned n, const double *x, double *gradient, /* NULL if not needed */ void *func_data);

void getLimits (std::vector<double> &, std::vector<double> &, int);
void bumpOptim(std::vector<double>, std::vector<double>, std::vector<double>, double *, double, unsigned int, nlopt_func, void *);
double calcErrorOptOneLog (unsigned, const double *, double *, void *);
double calcErrorOptLog (unsigned, const double *, double *, void *);
double calcErrorOptAllLog (unsigned, const double *, double *, void *);
double calcErrorOptPaperSiLog (unsigned, const double *, double *, void *);
double calcErrorSiLog_sepA (unsigned, const double *, double *, void *);
double calcErrorOptAllSiLog_sepA (unsigned, const double *, double *, void *);
void getLimits_sepA (std::vector<double> &, std::vector<double> &, int);
double calcErrorOptPaperSiLog_sepA (unsigned, const double *, double *, void *);
double calcErrorSiLog_sepA (unsigned, const double *, double *, void *);
double calcErrorOptPaperSiAllLog_sepA (unsigned n, const double *x, double *grad, void *data);

#endif /* defined(__UniformOptimization__Optimization__) */

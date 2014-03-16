//
//  ModelRunning.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ModelRunning__
#define __UniformOptimization__ModelRunning__

#include "ReactionCode.h"

#define autocrineT 10000
#define print_CV_err 2
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))
#define DoseTime 240
#define NfitCells 2


static const double times[] = {0, 0.5, 1, 5, 10}; ///< Times of kinetic measurements.
static const double Gass[] = {1.25, 0.25, 0.05, 0.01}; ///< Kinetic Gas6 doses.
static const double GassDose[] = {2.50, 1.25, 0.625, 0.3125, 0.15625, 0.078125, 0.00};

// Wrapping is outermost cell line, then Gas, then time
static const double pY[] = { ///< pY measurements on short time scales.
    1.000, 1.204, 1.140, 1.259, 1.108, 1.000, 1.699, 1.108, 1.180, 1.028,
    1.000, 1.120, 1.004, 0.932, 0.852, 1.000, 1.004, 0.916, 0.992, 0.936, // A172
    1.000, 1.211, 1.222, 2.037, 2.942, 1.000, 1.400, 1.211, 1.332, 1.118,
    1.000, 1.111, 0.962, 1.395, 1.184, 1.000, 1.116, 0.700, 1.043, 1.006, // A549
    1.000, 1.147, 1.383, 1.501, 2.149, 1.000, 1.306, 0.919, 1.468, 1.133,
    1.000, 0.952, 1.258, 1.169, 0.904, 1.000, 0.882, 0.919, 0.860, 0.930, // U87
    1.000, 1.298, 1.560, 1.610, 1.204, 1.000, 2.165, 1.739, 1.605, 1.271,
    1.000, 1.850, 1.730, 1.354, 1.026, 1.000, 1.591, 1.371, 1.235, 0.897};// BT549

static const double pYerror[] = { ///< Error for short time scale pY measurements.
    0.102, 0.084, 0.045, 0.366, 0.138, 0.102, 0.246, 0.109, 0.166, 0.080,
    0.102, 0.161, 0.054, 0.008, 0.080, 0.102, 0.146, 0.084, 0.137, 0.042, // A172
    0.087, 0.034, 0.058, 0.008, 0.035, 0.087, 0.239, 0.097, 0.096, 0.133,
    0.087, 0.183, 0.155, 0.335, 0.067, 0.087, 0.094, 0.138, 0.058, 0.197, // A549
    0.041, 0.191, 0.288, 0.119, 0.068, 0.041, 0.309, 0.060, 0.210, 0.048,
    0.041, 0.091, 0.218, 0.029, 0.075, 0.041, 0.139, 0.076, 0.135, 0.183, // U87
    0.294, 0.116, 0.143, 0.034, 0.125, 0.294, 0.194, 0.077, 0.054, 0.217,
    0.294, 0.050, 0.194, 0.373, 0.135, 0.294, 0.144, 0.035, 0.057, 0.069};// BT549

static const double pYdose[][7] = {
    {0.097755178, 0.104591926, 0.058509424, 0.035308947, 0.020346485, 0.019266966, 0.017639228}, // A172
    {0.034696084, 0.031362807, 0.025562428, 0.020713842, 0.013218411, 0.011107913, 0.012675294}, // A549
    {0.137176522, 0.129695581, 0.079213029, 0.072137190, 0.033326061, 0.025154533, 0.014065899}, // U87
    {0.087530778, 0.107545744, 0.121762241, 0.122572025, 0.105076825, 0.101415904, 0.108748274}};// BT549

static const double DoseError[][7] = {
    {0.006607048, 0.009925505, 0.002981499, 0.003968828, 0.003270587, 0.00318326,  0.001418565}, // A172
    {0.003295993, 0.003288019, 0.00346509,  0.003814638, 0.002328042, 0.002823496, 0.001398771}, // A549
    {0.026841367, 0.012427671, 0.003871628, 0.004112299, 0.007232776, 0.006543342, 0.003681512}, // U87
    {0.013008556, 0.028788423, 0.019185641, 0.028070405, 0.024978345, 0.016382209, 0.015430528}};// BT549

static const double DoseTot[][7] = {
    {1398, 1149, 1454, 1561, 2132, 1424, 1717}, // A172
    {1734, 1960, 1687, 1429, 2432, 2181, 2672}, // A549
    {425, 394, 470, 322, 452, 475, 544},        // U87
    {306, 294, 199, 247, 247, 278, 207}};       // BT549

static const double DoseTotErr[][7] = {
    {63 , 105, 62, 170, 267, 91 , 98}, // A172
    {224, 242, 98, 186, 198, 113, 78}, // A549
    {28, 37, 49, 22, 33, 24, 20},      // U87
    {20, 50, 13, 30, 13, 16, 21}};     // BT549



double calcError (param_type);
void errorLogger (std::exception *);
void initState(N_Vector, struct rates, double);
void diffusionSolution(double *, double, double *, int, double, double *, double *, int, double *, double, double);
void calcErrorRef (param_type, double *, std::atomic<bool> *);
void calcProfileSet (double *, double *, struct rates, int, double, double, double, int);
double calcErrorOneLine (struct rates, size_t, double);

#endif /* defined(__UniformOptimization__ModelRunning__) */

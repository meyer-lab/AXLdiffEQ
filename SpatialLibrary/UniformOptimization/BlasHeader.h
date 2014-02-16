/**
 *   \file BlasHeader.h
 *   FiniteDiffGSL
 *
 *
 *   Copyright (c) 2013 Aaron Meyer. All rights reserved.
 */

#ifndef __UniformOptimization__HelperFunctions__
#define __UniformOptimization__HelperFunctions__

#include <array>
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_impl.h>
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_lapack.h>
#include <cvode/cvode_diag.h>

using namespace std;

#define Nspecies 14
#define autocrineT 10000
#define maxR 1.0
#define print_CV_err 2
#define rel_tol 1E-4 // Relative tolerance of the solver
#define abs_tol 1E-4
#define SELECTED_MODEL AXL_react
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))
#define DoseTime 240
#define NfitCells 2

extern double diffD[Nspecies];

typedef array< double , 19 > param_type;

struct rates {
    double Binding1;   ///< Forward binding rate for Ig1
    double Binding2;   ///< Forward binding rate for Ig2
    double Unbinding1; ///< Reverse binding rate for Ig1
    double Unbinding2; ///< Reverse binding rate for Ig2
    double xFwd1;      ///< Reaction 1 forward rate.
    double xRev1;      ///< Reaction 1 reverse rate.
    double xFwd2;      ///< Reaction 2 forward rate.
    double xRev2;      ///< Reaction 2 reverse rate.
    double xFwd3;      ///< Reaction 3 forward rate.
    double xRev3;      ///< Reaction 3 reverse rate.
    double xFwd4;      ///< Reaction 4 forward rate.
    double xRev4;      ///< Reaction 4 reverse rate.
    double xFwd5;      ///< Reaction 5 forward rate.
    double xRev5;      ///< Reaction 5 reverse rate.
    double xFwd6;      ///< Reaction 6 forward rate.
    double xRev6;      ///< Reaction 6 reverse rate.
    double scaleA;     ///< Signaling capacity of non-D2 species.
    double expression; ///< AXL expression rate.
    double internalize;///< Non-pY species internalization rate.
    double pYinternalize;///< pY species internalization rate.
    double kRec;       ///< Recycling rate.
    double kDeg;       ///< Degradation rate.
    double fElse;      ///< Recycling fraction for non-D2 species.
    double fD2;        ///< Recycling fraction for D2.
};



static const double times[] = {0, 0.5, 1, 5, 10}; ///< Times of kinetic measurements.
static const double Gass[] = {1.25, 0.25, 0.05, 0.01}; ///< Kinetic Gas6 doses.
static const double GassDose[] = {2.50, 1.25, 0.625, 0.3125, 0.15625, 0.078125, 0.00};

// Wrapping is outermost cell line, then Gas, then time
static const double pY[] = { ///< pY measurements on short time scales.
    1.000, 1.204, 1.140, 1.259, 1.108, 1.000, 1.699, 1.108, 1.180, 1.028,
    1.000, 1.120, 1.004, 0.932, 0.852, 1.000, 1.004, 0.916, 0.992, 0.936, // A172
    1.000, 1.211, 1.222, 2.037, 2.942, 1.000, 1.400, 1.211, 1.332, 1.118,
    1.000, 1.111, 0.962, 1.395, 1.184, 1.000, 1.116, 0.700, 1.043, 1.006}; // A549

static const double pYerror[] = { ///< Error for short time scale pY measurements.
    0.102, 0.084, 0.045, 0.366, 0.138, 0.102, 0.246, 0.109, 0.166, 0.080,
    0.102, 0.161, 0.054, 0.008, 0.080, 0.102, 0.146, 0.084, 0.137, 0.042, // A172
    0.087, 0.034, 0.058, 0.008, 0.035, 0.087, 0.239, 0.097, 0.096, 0.133,
    0.087, 0.183, 0.155, 0.335, 0.067, 0.087, 0.094, 0.138, 0.058, 0.197}; // A549

static const double pYerrorAvg[] = {
    0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177,
    0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, 0.1177, // A172
    0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138,
    0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138, 0.1138}; // A549

static const double pYdose[][7] = {
    {0.097755178, 0.104591926, 0.058509424, 0.035308947, 0.020346485, 0.019266966, 0.017639228}, // A172
    {0.034696084, 0.031362807, 0.025562428, 0.020713842, 0.013218411, 0.011107913, 0.012675294}}; // A549

static const double DoseError[][7] = {
    {0.006607048, 0.009925505, 0.002981499, 0.003968828, 0.003270587, 0.00318326,  0.001418565}, // A172
    {0.003295993, 0.003288019, 0.00346509,  0.003814638, 0.002328042, 0.002823496, 0.001398771}}; // A549

static const double DoseErrorAvg[][7] = {
    {0.004479, 0.004479, 0.004479, 0.004479, 0.004479, 0.004479, 0.004479}, // A172
    {0.002916, 0.002916, 0.002916, 0.002916, 0.002916, 0.002916, 0.002916}}; // A549

static const double DoseTot[][7] = {
    {1398, 1149, 1454, 1561, 2132, 1424, 1717}, // A172
    {1734, 1960, 1687, 1429, 2432, 2181, 2672}}; // A549

static const double DoseTotErr[][7] = {
    {63 , 105, 62, 170, 267, 91 , 98}, // A172
    {224, 242, 98, 186, 198, 113, 78}}; // A549

int solver_setup (void **, N_Vector , void *, N_Vector , CVRhsFn);
int check_flag(void *, char *, int );

extern "C" int matlabEntry(double * , double * , int);
extern "C" int calcProfileMatlab(double *, double *, double *, int, double, double, double, int);
extern "C" int matlabDiffTPS(double *, double , double *, int , double , double *, double *, int, double *, double, double);
extern "C" int matlabDiffTPS_pY(double *, double , double *, int , double , double *, double *, int, double *, double, double, int);
extern "C" int matlabDiffTPS_pYavg(double *, double, double *, int, double , double *, double *, int, double *, double, double, int);
extern "C" void rEntry(double *, const double *);
double calcError (param_type);

#endif /* defined(__UniformOptimization__HelperFunctions__) */
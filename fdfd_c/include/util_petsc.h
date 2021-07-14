#ifndef _UTILS
#define _UTILS

// ----- Defines utility functions for CSparse vectors and matrices ------
// -----------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>

#include <petscksp.h>

#include "domain.h"

using namespace std;

// --- Define fundamental constants ---
#define PI 3.14159
#define C0 3e8
#define EPS0 8.85418781762039e-12
#define MU0 1.25663706212e-6
#define ETA0 376.7303
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef std::chrono::steady_clock Clock;

typedef struct {
    double k0;
    double kx;
    double ky;
    double Wy;
} FDparams;

// --- Declare global variables ---
extern Vec b_SFTF, fsrc; 
extern Vec b_NL1, b_NL2;
extern Vec ez_1, ez_1cc, ez_2;
extern Vec epszz, muxx, muyy; 
extern Mat Ae1, Ae2, Q; 
extern const int Nc;

// --- Define source parameters ---
static const double A0 = 4.0e6; // Incident FW amplitude [V/m]
static const double omeg = 1.7715733e15; // source frequency [1.064 microns]
static const double thetai = 0.0; // source angle

// --- Declare utility functions ---
int fsource(double A0, double kx, double ky);
int epszzfunc(double epsr0, double epsr1);
int mufunc();
int finitediffs(double k0, double ky, double Wy);
int defineQ(int SFx);
int defineAe(Mat *Ap, double k0, double ky, double Wy);


#endif 

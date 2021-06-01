#ifndef _HELPERS
#define _HELPERS

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include <complex.h>
#include <omp.h>

#include "parameters.h"

// --- Constants and simple macros
#define PI 3.14159
#define C0 3e8                      // Speed of light [m / s]
#define EPS0 8.85418781762039e-12   // Vacuum permittivity [C / V m]
#define MU0 1.25663706212e-6        // Vacuum permeability [N / A^2]
#define ETA0 376.730313668          // Vacuum impedance [Ohms]
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
typedef std::chrono::steady_clock Clock;

// Define global variables
extern double gi2[NX], gi3[NX];                     // PML parameters (see Sullivan, Ch. 3)
extern double fi1[NX], fi2[NX], fi3[NX];
extern double ga[NX], gb[NX];                           // Dz to Ez  
extern double Ix[NX];
extern double Dx[NX], Hy[NX];   // Field arrays
extern double Ex[NX];                           // Actual field in z-direction   
extern double chi3;
extern double Ew_re[NF][NX], Ew_im[NF][NX];     // Real and imaginary frequency domain output.
extern double Fsrc_re[NF], Fsrc_im[NF]; // Real and imaginary frequency domain source.
extern double omeg[NF];
extern double amp_in[NF];

// --- Function declarations ---
double src(double t);
void pmldef();      // Defines parameters for PML regions
void savegrid();    // Save XY-grid for plotting later
void initArrays(double dt);
double nfdtdsteps(int N, double T, double dt);    // Save XY-grid for plotting later

#endif

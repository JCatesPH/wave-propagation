#ifndef _HELPERS
#define _HELPERS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "parameters.h"

// --- Constants and simple macros
#define PI 3.14159
#define C0 3e8                      // Speed of light [m / s]
#define EPS0 8.85418781762039e-12   // Vacuum permittivity [C / V m]
#define MU0 1.25663706212e-6        // Vacuum permeability [N / A^2]
#define ETA0 376.730313668          // Vacuum impedance [Ohms]
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

// Define global variables
extern float gi2[NX], gi3[NX];                     // PML parameters (see Sullivan, Ch. 3)
extern float gj2[NY], gj3[NY];
extern float fi1[NX], fi2[NX], fi3[NX];
extern float fj1[NY], fj2[NY], fj3[NY];
extern float ga[NX][NY];                           // Dz to Ez  
extern float Ihx[NX][NY], Ihy[NX][NY];
extern float Dz[NX][NY], Hx[NX][NY], Hy[NX][NY];   // Field arrays
extern float Ez[NX][NY];                           // Actual field in z-direction   

// --- Function declarations ---
float src(float t);
void pmldef();      // Defines parameters for PML regions
void savegrid();    // Save XY-grid for plotting later
float nfdtdsteps(int nsteps, float T, float dt);    // Save XY-grid for plotting later

#endif

#ifndef _DOMAIN
#define _DOMAIN

// --- Domain size ---
#define NX 2400     // Number of x grid points
#define NY 101      // Number of y grid points
#define LX 40       // Num cells in PML
#define BX 120      // Num of cells in "buffer" (SF but not PML)

#define DX 6.65e-9     // x grid size
#define DY 6.65e-9     // y grid size

// --- PML Parameters ---
#define SX_P 3      // Polynomial order
#define SX_R 1.e-12  // Desired reflectivity coefficient
#define SX_M 3.0    // Max value of s_x

// --- Material Parameters
static const double d33 = 13.6e-12; // Susceptibility tensor element [m/V]; from OSA B paper
static const double epsr1 = 4.64834; // Extraordinary index of refraction for first harmonic
static const double epsr2 = 4.99165; // "..." for second harmonic
static const int idx1 = 320;
static const int idx2 = 2080;

#endif

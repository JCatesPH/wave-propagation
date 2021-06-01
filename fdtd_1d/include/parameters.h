#ifndef _DOMAIN
#define _DOMAIN

// --- Domain size ---
#define NX 20000      // Number of x grid points
#define NOUT 25     // Number of time steps where output is stored
#define DOUT 500    // Spacing between time steps saved
#define LX 10        // Num cells in PML
#define BX 50      // Num of cells in "buffer" (SF but not PML)
#define DX 5e-9     // x grid size
#define NF 2000        // Number of frequencies

// --- PML Parameters ---
#define SX_P 3          // Polynomial order
#define SX_R 1e-10      // Desired reflectivity coefficient
#define SX_M 4.0        // Max value of s_x


/* - Source parameters - */
//  Currently point source at (ic,jc)
static const int isrc     = NX/5;    // Index of source location
static const double I0    = 1.0;  //50e16;       // Mean pump intensity [W / m^2]
static const double t0    = 1e-13;      // Time of pulse peak [s]
static const double taup  = 50e-15;    // Pulse duration (FWHM) for fund. pulse, second harmonic has half duration [s]
static const double omeg0 = 1.884e15; // Angular freq. of fund. {1.884e15}
static const double R     = 0.0;         // Ratio of fund. and 2nd harmonic
static const double phi   = 0.0;         // Phase offset of 2nd harmonic

// --- Material Parameters ---
static const int    idie = NX/2;    // Index where dielectric starts
static const double epsz = 4.0;     // Relative permittivity
static const double sigma = 0.0;     // Conductivity
static const double n2 = 1e-23;      // Nonlinear index [m^2 / W]
static const double alpha = 3/4;     // Ratio of Kerr and Raman effect

#endif

#ifndef _DOMAIN
#define _DOMAIN

// --- Domain size ---
#define NX 500      // Number of x grid points
#define NY 500      // Number of y grid points
#define NOUT 75     // Number of time steps where output is stored
#define DOUT 5     // Spacing between time steps saved
#define LX 10        // Num cells in PML
#define LY 10        // Num cells in PML
#define BX 120      // Num of cells in "buffer" (SF but not PML)

#define DX 0.01     // x grid size
#define DY 0.01     // y grid size

// --- PML Parameters ---
#define SX_P 3          // Polynomial order
#define SX_R 1e-10      // Desired reflectivity coefficient
#define SX_M 4.0        // Max value of s_x


/* - Source parameters - */
//  Currently point source at (ic,jc)
static const int ic = NX/2 - 5;    // First index of source location
static const int jc = NY/2 - 5;    // Second index of source location
static const float I0 = 50e16;       // Mean pump intensity [W / m^2]
static const float t0 = 4e-10;      // Time of pulse peak [s]
static const float taup = 1e-10;    // Pulse duration (FWHM) for fund. pulse, second harmonic has half duration [s]
static const float omeg0 = 0.0; // Angular freq. of fund. {1.884e15}
static const float R = 0.0;         // Ratio of fund. and 2nd harmonic
static const float phi = 0.0;         // Phase offset of 2nd harmonic

// --- Material Parameters ---
static const float epsz = 1.0;     // Relative permittivity

#endif

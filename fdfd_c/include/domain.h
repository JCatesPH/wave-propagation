#ifndef _DOMAIN
#define _DOMAIN

// --- Domain size ---
#define NX 1000     // Number of x grid points
#define NY 200      // Number of y grid points
#define LX 40       // Num cells in PML
#define BX 120      // Num of cells in "buffer" (SF but not PML)

#define DX 5e-8     // x grid size
#define DY 5e-8     // y grid size

// --- PML Parameters ---
#define SX_P 3      // Polynomial order
#define SX_R 1e-10  // Desired reflectivity coefficient
#define SX_M 4.0    // Max value of s_x

// --- Material Parameters
const double chi3 = 5e-14;
const double epsr = 1.0;

#endif

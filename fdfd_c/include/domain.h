#ifndef DOMAIN
#define DOMAIN

// --- Domain size ---
#define NX 1000     // Number of x grid points
#define NY 200      // Number of y grid points
#define LX 40       // Num cells in PML
#define BX 120      // Num of cells in "buffer" (SF but not PML)

const float dx = 5e-8; // x grid size
const float dy = 5e-8; // y grid size

// --- PML Parameters ---
const int sx_p = 3; // Polynomial order
const float sx_R = 1e-10; // Desired reflectivity coefficient
const float sx_max = 4.0; // Max value of s_x

#endif
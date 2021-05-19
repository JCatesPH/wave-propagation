#ifndef _UTILS
#define _UTILS

// ----- Defines utility functions for CSparse vectors and matrices ------
// -----------------------------------------------------------------------

#include "cs.h"

// --- Define fundamental constants ---
#define PI 3.14159
#define C0 3e8
#define EPS0 8.85418781762039e-12
#define MU0 1.25663706212e-6
#define ETA0 376.7303
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

// --- Declare utility functions ---
cs_cl *cs_cl_scale (cs_cl *A, double alpha);
void cvecprint(char* path, cs_complex_t *vec, int length);
cs_complex_t *cvecdiff(cs_complex_t *x, cs_complex_t *y);
cs_complex_t *cvecadd(cs_complex_t *x, cs_complex_t *y);

#endif 

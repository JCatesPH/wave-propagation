#include "include/util.h"
#include "include/domain.h"

// --- Sparse matrix scaling function ---
cs_cl *cs_cl_scale (cs_cl *A, double alpha)
{
    cs_cl *C;
    cs_long_t n, triplet, nn, p, nz, *Ap, *Ai, *Cp, *Ci ;
    cs_complex_t *Ax, *Cx ;
    if (!A || !A->x) return (NULL) ;    /* return if A NULL or pattern-only */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    triplet = (A->nz >= 0) ;            /* true if A is a triplet matrix */
    nz = triplet ? A->nz : Ap [n] ;
    C = cs_cl_spalloc (A->m, n, A->nzmax, 1, triplet) ;
    if (!C) return (NULL) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    nn = triplet ? nz : (n+1) ;
    for (p = 0 ; p < nz ; p++) Ci [p] = Ai [p] ;
    for (p = 0 ; p < nn ; p++) Cp [p] = Ap [p] ;
    for (p = 0 ; p < nz ; p++) Cx [p] = alpha * Ax [p];
    if (triplet) C->nz = nz ;
    return (C) ;
}


// --- Utility function to print complex vector ---
void cvecprint(char* path, cs_complex_t *vec, int length){
    FILE *fp;
    fp = fopen(path, "w");
    for (int i=0; i<length; i++) {
        fprintf(fp, "%12.9e, %12.9e\n", creal(vec[i]), cimag(vec[i]));
    }
}


// --- Define difference between complex vectors ---
cs_complex_t *cvecdiff(cs_complex_t *x, cs_complex_t *y){
    // Compute: z = x - y
    static cs_complex_t z[NX*NY];
    for (int i=0; i<NX*NY; i++){
        z[i] = x[i] - y[i];
    }
    return z;
}

// --- Define addition of complex vectors ---
cs_complex_t *cvecadd(cs_complex_t *x, cs_complex_t *y){
    // Compute: z = x - y
    static cs_complex_t z[NX*NY];
    for (int i=0; i<NX*NY; i++){
        z[i] = x[i] + y[i];
    }
    return z;
}

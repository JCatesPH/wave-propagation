#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "include/cs.h"
#include "include/domain.h"

// --- Define global variables ---
#define PI 3.14159
#define C0 3e8
#define EPS0 8.85418781762039e-12
#define MU0 1.25663706212e-6
#define ETA0 376.7303

//#define VERBOSE // Uncomment to print additional debugging messages.

const int Nc = NX*NY;
cs_complex_t b[NX*NY], fsrc[NX*NY]; // Source vector
double X[NX]; // x grid points
double Y[NY]; // y grid points
double _Complex sx[LX];
cs_complex_t epszz[NX*NY], muxx[NX*NY], muyy[NX*NY]; 
cs_cl *Dxe, *Dye, *Dxh, *Dyh, *Ae, *Q; 


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
        fprintf(fp, "%12.9e, %+12.9e\n", creal(vec[i]), cimag(vec[i]));
    }
}

// --- Define source function ---
void fsource(float kx, float ky){
    printf("Defining source.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            fsrc[i*NY+j] = cexpf(I*(kx*X[i] + ky*Y[j]));
            b[i*NY+j] = 0.0;
        }
    }
    printf("Finished defining source.\n");
}

// --- Define permittivity function ---
void epszzfunc(){
    double sigmax = -(sx_p + 1) * log(sx_R) / (2 * ETA0 * LX);
    double frac, s0x, sigx;

    // Define PML region and change eps accordingly.
    printf("Defining permittivity tensor.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            //printf("%d, %d\n", i, j);
            if (i < LX) {
                frac = (float) (LX-i)/LX;
                s0x = 1.0 + sx_max * powf(frac, sx_p);
                sigx = sigmax * powf(sin(PI*frac/2), 2);
                sx[i] = s0x * (1.0 - I * ETA0 * sigx);
                epszz[i*NY+j] = sx[i];
            }
            else if (i > NX-LX) {
                epszz[i*NY+j] = sx[NX-i];
            }
            else{
                epszz[i*NY+j] = 1.0 + I*0.0;
            }
        }
    }
    printf("Finished defining permittivity tensor.\n");
}

// --- Define permeability tensors ---
void mufunc() {

    // Define PML region and change mu accordingly.
    printf("Defining permeability tensors.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            //printf("%d, %d\n", i, j);
            if (i < LX) {
                muxx[i*NY+j] = 1/sx[i];
                muyy[i*NY+j] = sx[i];
            }
            else if (i > NX-LX) {
                muxx[i*NY+j] = 1/sx[NX-i];
                muyy[i*NY+j] = sx[NX-i];
            }
            else {
                muxx[i*NY+j] = 1.0 + I*0.0;
                muyy[i*NY+j] = 1.0 + I*0.0;
            }
            if (isnan(creal(muxx[i*NY+j])) ||  isnan(cimag(muxx[i*NY+j]))) {
                fprintf(stderr, "Error: NaN encountered in muxx.\n"); exit(EXIT_FAILURE);
            }
            if (isnan(creal(muyy[i*NY+j])) ||  isnan(cimag(muyy[i*NY+j]))) {
                fprintf(stderr, "Error: NaN encountered in muyy.\n"); exit(EXIT_FAILURE);
            }
        }
    }
    printf("Finished defining permeability tensors.\n");
}

// --- Define finite-difference matrices ---
void finitediffs(float k0, float ky, float Wy) {
    int cserr;

    // Need to set the diagonals for finite differences.
    printf("Defining finite-difference matrices.\n");
    for (int i=0; i<Nc; i++){
        for (int j=0; j<Nc; j++) {
            // Set diagonal elements
            if (i==j){
                cserr = cs_cl_entry(Dxe, i, j, -1/(k0*dx));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dxe.\n"); exit(EXIT_FAILURE);}
                cserr = cs_cl_entry(Dye, i, j, -1/(k0*dy));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
            }
            // Set first superdiagonal
            else if (j==i+1){
                cserr = cs_cl_entry(Dxe, i, j, 1/(k0*dx));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dxe.\n"); exit(EXIT_FAILURE);}
            }
            // Set second band for Dye
            else if (j==i+NX){
                cserr = cs_cl_entry(Dye, i, j, 1/(k0*dx));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
            }
            // Impose Floquet boundary condition
            else if (j==i-NX*(NY-1)) {
                cserr = cs_cl_entry(Dye, i, j, cexpf(I*ky*Wy));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
            }
        }
    }
    //cs_cl_print(Dxe, 0);
    Dxe = cs_cl_compress(Dxe);
    Dye = cs_cl_compress(Dye);
    printf("Finished defining finite-difference matrices.\n");
}

// --- Define full linear operator ---
void defineAe(float k0, float ky, float Wy) {
    int cserr;

    // Create local matrices for pieces of product
    cs_cl *muxxmat = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);
    cs_cl *muyymat = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);
    cs_cl *epszzmat = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);
    // Define diagonal entries with mu vectors
    for (int i=0; i<Nc; i++){
        cserr = cs_cl_entry(muxxmat, i, i, 1/muxx[i]);
        if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for muxx.\n"); exit(EXIT_FAILURE);}
        cserr = cs_cl_entry(muyymat, i, i, 1/muyy[i]);
        if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for muyy.\n"); exit(EXIT_FAILURE);}
        cserr = cs_cl_entry(epszzmat, i, i, epszz[i]);
        if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for epszz.\n"); exit(EXIT_FAILURE);}
    }
    // Compress format for multiplication
    muxxmat = cs_cl_compress(muxxmat);
    muyymat = cs_cl_compress(muyymat);
    epszzmat = cs_cl_compress(epszzmat);

    // Set vector to be all ones
    Ae = cs_cl_multiply(Dxh, cs_cl_multiply(muyymat, Dxe));
    Ae = cs_cl_add(Ae, cs_cl_multiply(Dyh, cs_cl_multiply(muxxmat, Dye)), 1.0, 1.0);
    Ae = cs_cl_add(Ae, epszzmat, 1.0, 1.0);

    // Print brief versions of matrices to bug check
    #ifdef VERBOSE
        printf("===============\nDxe\n");
        cs_cl_print(Dxe, 1);
        printf("\n===============\nDye\n");
        cs_cl_print(Dye, 1);
        printf("\n===============\nDxh\n");
        cs_cl_print(Dxh, 1);
        printf("\n===============\nDyh\n");
        cs_cl_print(Dyh, 1);
        printf("===============\nmuxx\n");
        cs_cl_print(muxxmat, 1);
        printf("===============\nmuyy\n");
        cs_cl_print(muyymat, 1);
        printf("===============\nepszz\n");
        cs_cl_print(epszzmat, 1);
        printf("\n===============\nAe\n");
        cs_cl_print(Ae, 1);
    #endif

    // Clean up everything
    cs_cl_spfree(Dxh);
    cs_cl_spfree(Dxe);
    cs_cl_spfree(Dyh);
    cs_cl_spfree(Dye);
    cs_cl_spfree(epszzmat);
    cs_cl_spfree(muxxmat);
    cs_cl_spfree(muyymat);
    printf("Finished defining Ae and freed memory.\n");
}


// --- Define masking matrix for SF/TF ---
void defineQ(int SFx){
    int cserr;
    for (int i=0; i<Nc; i++) {
        if ((i % NY) < SFx) {
            cserr = cs_cl_entry(Q, i, i, 1.0+I*0.0);
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Q.\n"); exit(EXIT_FAILURE);}
        }
    }
    Q = cs_cl_compress(Q);
    //cs_cl_print(Q, 0);
}

// ------- MAIN -------
int main ()
{
    // --- Allocate memory for variables ---
    Dxe = cs_cl_spalloc(Nc, Nc, 3*Nc, 1, 1);
    Dye = cs_cl_spalloc(Nc, Nc, 3*Nc, 1, 1);
    Ae = cs_cl_spalloc(Nc, Nc, 5*Nc, 1, 0);
    Q = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);

    // --- Define source parameters ---
    float omeg = 9.4248e+14; // source frequency
    float thetai = 0.0; // source angle
    float k0 = omeg / C0; // source wave number
    float kx = k0 * cos(thetai); 
    float ky = k0 * sin(thetai);

    // --- Define derived domain variables ---
    float Wy = NY * dy; // width of y domain (for Floquet)
    int SFx = NX - LX - BX;

    FILE *fp;
    char fname[40];

    // --- Set the grid point values ---
    sprintf(fname, "data/X.csv");
    fp = fopen(fname, "w");
    for (int i=0; i<NX; i++){
        X[i] = i*dx;
        fprintf(fp, "%12.9e\n", X[i]);
    }
    fclose(fp);

    sprintf(fname, "data/Y.csv");
    fp = fopen(fname, "w");
    for (int j=0; j<NY; j++){
        Y[j] = j*dy;
        fprintf(fp, "%12.9e\n", Y[j]);
    }
    fclose(fp);


    // --- Define and save the permittivity and permeability tensors ---
    epszzfunc();
    mufunc();
    cvecprint("data/epszz.csv", epszz, Nc);


    // --- Generate and check the source field ---
    fsource(kx, ky);
    cvecprint("data/src.csv", fsrc, Nc);
   
    // --- Generate the finite-difference operators ---
    
    finitediffs(k0, ky, Wy);
    Dxh = cs_cl_transpose(Dxe, 1);
    Dyh = cs_cl_transpose(Dye, 1);
    Dxh = cs_cl_scale(Dxh, -1.0);
    Dyh = cs_cl_scale(Dyh, -1.0);

    // --- Define total linear operator ---
    printf("Defining linear operator Ae.\n");
    defineAe(k0, ky, Wy);

    // --- Define masking matrix ---
    printf("Defining masking matrix Q.\n");
    defineQ(SFx);

    // --- Use masking matrix to find source vector ---
    cs_cl_gaxpy(cs_cl_add(cs_cl_multiply(Q, Ae), cs_cl_multiply(Ae, Q), 1.0, -1.0), fsrc, b);
    printf("Source vector found. Printing.. \n");
    cvecprint("data/b.csv", b, Nc);
    printf("Finished printing. \n");

    // --- Set tolerance and solve system with LU factorization ---
    double tol = 1;
    cs_cl_lusol(0, Ae, b, tol);
    printf("System solved.\n");

    // --- Print output for further plotting/analysis ---
    cvecprint("data/Ez.csv", b, Nc);

    // --- Free remaining memory ---
    cs_cl_spfree(Ae);
    cs_cl_spfree(Q);
    //return 0;

    printf("Exiting.\n\n");
    exit(EXIT_SUCCESS);
}

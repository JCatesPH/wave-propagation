#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

#include "domain.h"

// --- Define global variables ---
#define PI 3.14159
#define C0 3e8
#define EPS0 8.85418781762039e-12
#define MU0 1.25663706212e-6
#define ETA0 376.7303

int Nc = NX*NY;
gsl_vector_complex *b; // Source vector
float X[NX]; // x grid points
float Y[NY]; // y grid points
float _Complex sx[LX];
gsl_vector_complex *epszz;
gsl_vector_complex *muxx;
gsl_vector_complex *muyy;
gsl_spmatrix_complex *Dxe;
gsl_spmatrix_complex *Dye;
gsl_spmatrix_complex *Dxh;
gsl_spmatrix_complex *Dyh;
gsl_spmatrix_complex *Ae_COO; // COO indicates storage format for sparse matrix
gsl_spmatrix_complex *Q_COO;
gsl_spmatrix_complex *Ae_CSC; // CSC : Compressed Sparse Column format
gsl_spmatrix_complex *Q_CSC;

// --- Define source function ---
void fsource(float kx, float ky){
    printf("Defining source.\n");
    gsl_complex tmp;
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            GSL_SET_COMPLEX(&tmp, crealf(cexpf(I*(kx*X[i] + ky*Y[j]))), cimag(cexpf(I*(kx*X[i] + ky*Y[j]))));
            gsl_vector_complex_set(b, i*NY+j, tmp);
        }
    }
    printf("Finished defining source.\n");
}

// --- Define permittivity function ---
void epszzfunc(){
    float sigmax = -(sx_p + 1) * log(sx_R) / (2 * ETA0 * LX);
    float frac, s0x, sigx;

    // Set all elements to 1 to start.
    gsl_complex tmp;
    GSL_SET_COMPLEX(&tmp, 1.0, 0.0);
    gsl_vector_complex_set_all(epszz, tmp);

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
                GSL_SET_COMPLEX(&tmp, creal(sx[i]), cimag(sx[i]));
                gsl_vector_complex_set(epszz, i*NY+j, tmp);
            }
            else if (i > NX-LX) {
                GSL_SET_COMPLEX(&tmp, creal(sx[NX-LX-i]), cimag(sx[NX-LX-i]));
                gsl_vector_complex_set(epszz, i*NY+j, tmp);
            }
        }
    }
    printf("Finished defining permittivity tensor.\n");
}

// --- Define permeability tensors ---
void mufunc() {
    // Set all elements to 1 to start.
    gsl_complex tmp;
    GSL_SET_COMPLEX(&tmp, 1.0, 0.0);
    gsl_vector_complex_set_all(muxx, tmp);
    gsl_vector_complex_set_all(muyy, tmp);

    // Define PML region and change eps accordingly.
    printf("Defining permeability tensors.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            //printf("%d, %d\n", i, j);
            if (i < LX) {
                GSL_SET_COMPLEX(&tmp, creal(1/sx[i]), cimag(1/sx[i]));
                gsl_vector_complex_set(muxx, i*NY+j, tmp);
                GSL_SET_COMPLEX(&tmp, creal(sx[i]), cimag(sx[i]));
                gsl_vector_complex_set(muyy, i*NY+j, tmp);
            }
            else if (i > NX-LX) {
                GSL_SET_COMPLEX(&tmp, creal(1/sx[NX-LX-i]), cimag(1/sx[NX-LX-i]));
                gsl_vector_complex_set(muxx, i*NY+j, tmp);
                GSL_SET_COMPLEX(&tmp, creal(sx[NX-LX-i]), cimag(sx[NX-LX-i]));
                gsl_vector_complex_set(muyy, i*NY+j, tmp);
            }
        }
    }
    printf("Finished defining permeability tensors.\n");
}

// --- Define finite-difference matrices ---
void finitediffs(float k0, float ky, float Wy) {
    // Set all elements to 1 to start.
    gsl_complex tmp; float _Complex z;

    // Define PML region and change eps accordingly.
    printf("Defining finite-difference matrices.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            // Set diagonal elements
            if (i==j){
                GSL_SET_COMPLEX(&tmp, -1/(k0*dx), 0.0);
                gsl_spmatrix_complex_set(Dxe, i, j, tmp);
                gsl_spmatrix_complex_set(Dye, i, j, tmp);
            }
            // Set first superdiagonal
            else if (j==i+1){
                GSL_SET_COMPLEX(&tmp, 1/(k0*dx), 0.0);
                gsl_spmatrix_complex_set(Dxe, i, j, tmp);
            }
            // Set second band for Dye
            else if (j==i+NX){
                GSL_SET_COMPLEX(&tmp, 1/(k0*dx), 0.0);
                gsl_spmatrix_complex_set(Dye, i, j, tmp);
            }
            // Impose Floquet boundary condition
            else if (j==i-NX*(NY-1)) {
                z = cexpf(I*ky*Wy);
                GSL_SET_COMPLEX(&tmp, creal(z)/(k0*dx), cimag(z)/(k0*dx));
                gsl_spmatrix_complex_set(Dye, i, j, tmp);
            }
        }
    }
    printf("Finished defining finite-difference matrices.\n");
}

// --- Define full linear operator ---
void defineAe() {
    // Create temporary vector for divisions
    gsl_complex z; GSL_SET_COMPLEX(&z, 1.0, 0.0);
    gsl_vector_complex *tmp = gsl_vector_complex_alloc(Nc);
    // Create local matrices for pieces of product
    gsl_spmatrix_complex *Ap1 = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    gsl_spmatrix_complex *Ap2 = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    // Set vector to be all ones
    gsl_vector_complex_set_all(tmp, z);
    // Compute first division and multiplication
    gsl_vector_complex_div(tmp, muyy);
    gsl_spmatrix_complex_scale_columns(Dxh, muyy);
    gsl_spblas_cgemm(1.0, Dxh, Dxe, Ap1);
    // Compute second division and multiplication
    gsl_vector_complex_set_all(tmp, z);
    gsl_vector_complex_div(tmp, muxx);
    gsl_spmatrix_complex_scale_columns(Dyh, muxx);
    gsl_spblas_cgemm(1.0, Dyh, Dye, Ap2);
    gsl_spmatrix_complex_add(Ae_COO, Ap1, Ap2);
    for (int i=0; i<Nc; i++){
        z = gsl_spmatrix_complex_get(Ae_COO, i, i);
        z = gsl_complex_add(z, gsl_vector_complex_get(epszz, i));
        gsl_spmatrix_complex_set(Ae_COO, i, i, z);
    }

    // Clean up everything
    gsl_vector_complex_free(tmp);
    gsl_spmatrix_complex_free(Ap1);
    gsl_spmatrix_complex_free(Ap2);
    gsl_spmatrix_complex_free(Dxh);
    gsl_spmatrix_complex_free(Dxe);
    gsl_spmatrix_complex_free(Dyh);
    gsl_spmatrix_complex_free(Dye);
    gsl_vector_complex_free(epszz);
    gsl_vector_complex_free(muxx);
    gsl_vector_complex_free(muyy);
}


// --- Define masking matrix for SF/TF ---
void defineQ(int SFx){
    gsl_complex z; 
    for (int i=0; i<Nc; i++) {
        if ((i % NY) < SFx) {
            GSL_SET_COMPLEX(&z, 1.0, 0.0);
            gsl_spmatrix_complex_set(Q_COO, i, i, z);
        }
        else {
            GSL_SET_COMPLEX(&z, 0.0, 0.0);
            gsl_spmatrix_complex_set(Q_COO, i, i, z);
        }
    }
}

// ------- MAIN -------
int main ()
{
    // --- Allocate memory for variables ---
    b = gsl_vector_complex_alloc(Nc);
    epszz = gsl_vector_complex_alloc(Nc);
    muxx = gsl_vector_complex_alloc(Nc);
    muyy = gsl_vector_complex_alloc(Nc);
    Dxe = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    Dye = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    Dxh = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    Dyh = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    Ae_COO = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    Q_COO = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_COO);
    Ae_CSC = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_CSC);
    Q_CSC = gsl_spmatrix_complex_alloc_nzmax(Nc, Nc, 3*Nc, GSL_SPMATRIX_CSC);

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


    // --- Define and save the permittivity tensor ---
    epszzfunc();
    sprintf(fname, "data/epszz.csv");
    fp = fopen(fname, "w");
    gsl_vector_complex_fprintf(fp, epszz, "%.5g");


    // --- Generate and check the source field ---
    fsource(kx, ky);
    sprintf(fname, "data/src.csv");
    fp = fopen(fname, "w");
    gsl_vector_complex_fprintf(fp, b, "%.5g");


    // --- Generate the finite-difference operators ---
    finitediffs(k0, ky, Wy);
    gsl_spmatrix_complex_transpose_memcpy(Dxh, Dxe);
    gsl_spmatrix_complex_transpose_memcpy(Dyh, Dye);


    // --- Define total linear operator ---
    defineAe();

    // --- Define masking matrix ---
    defineQ(SFx);

    // --- Change sparse matrix format to enable mult ---
    gsl_spmatrix_complex_csc(Ae_CSC, Ae_COO);
    gsl_spmatrix_complex_csc(Q_CSC, Q_COO);

    // --- Free old sparse matrix memory ---
    gsl_spmatrix_complex_free(Ae_COO);
    gsl_spmatrix_complex_free(Q_COO);


    // --- Use masking matrix to find source vector ---


    return 0;
}
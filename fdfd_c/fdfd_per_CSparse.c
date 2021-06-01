#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "include/cs.h"
#include "include/domain.h"
#include "include/util.h"

//#define VERBOSE // Uncomment to print additional debugging messages.
#define INCHECK // Uncomment to exit before solver is invoked. For checking inputs.

const int Nc = NX*NY;
cs_complex_t b[NX*NY], fsrc[NX*NY]; // Source vector
double X[NX]; // x grid points
double Y[NY]; // y grid points
double _Complex sx[LX];
cs_complex_t epszz[NX*NY], muxx[NX*NY], muyy[NX*NY]; 
cs_cl *Dxe, *Dye, *Dxh, *Dyh, *Ae, *Q; 

// --- Define source function ---
void fsource(double kx, double ky){
    printf("Defining source.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            fsrc[i+j*NX] = cexp(I*(kx*X[i] + ky*Y[j]));
            b[i+j*NX] = 0.0;
        }
    }
    printf("Finished defining source.\n");
}

// --- Define permittivity function ---
void epszzfunc(){
    double sigmax = -(SX_P + 1) * log(SX_R) / (2 * ETA0 * LX);
    #ifdef VERBOSE
        printf("sigmax = %f\n", sigmax);
    #endif
    double frac, s0x, sigx;

    // Define PML region and change eps accordingly.
    printf("Defining permittivity tensor.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            //printf("%d, %d\n", i, j);
            if (i < LX) {
                frac = (double) (LX-i)/LX;
                s0x = 1.0 + SX_M * powf(frac, SX_P);
                sigx = sigmax * powf(sin(PI*frac/2.0), 2.0);
                sx[i] = s0x * (1.0 - I * ETA0 * sigx);
                epszz[i+j*NX] = sx[i];
            }
            else if (i > NX-LX) {
                epszz[i+j*NX] = sx[NX-i-1];
            }
            else{
                epszz[i+j*NX] = 1.0 + I*0.0;
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
                muxx[i+j*NX] = 1/sx[i];
                muyy[i+j*NX] = sx[i];
            }
            else if (i > NX-LX) {
                muxx[i+j*NX] = 1/sx[NX-i-1];
                muyy[i+j*NX] = sx[NX-i-1];
            }
            else {
                muxx[i+j*NX] = 1.0 + I*0.0;
                muyy[i+j*NX] = 1.0 + I*0.0;
            }
            if (isnan(creal(muxx[i+j*NX])) ||  isnan(cimag(muxx[i+j*NX]))) {
                fprintf(stderr, "Error: NaN encountered in muxx.\n"); exit(EXIT_FAILURE);
            }
            if (isnan(creal(muyy[i+j*NX])) ||  isnan(cimag(muyy[i+j*NX]))) {
                fprintf(stderr, "Error: NaN encountered in muyy.\n"); exit(EXIT_FAILURE);
            }
        }
    }
    printf("Finished defining permeability tensors.\n");
}

// --- Define finite-difference matrices ---
void finitediffs(double k0, double ky, double Wy) {
    int cserr;

    // Need to set the diagonals for finite differences.
    printf("Defining finite-difference matrices.\n");
    /* OLD SLOW METHOD
    for (int i=0; i<Nc; i++){
        for (int j=0; j<Nc; j++) {
            // Set diagonal elements
            if (i==j){
                cserr = cs_cl_entry(Dxe, i, j, -1/(k0*DX));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dxe.\n"); exit(EXIT_FAILURE);}
                cserr = cs_cl_entry(Dye, i, j, -1/(k0*DY));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
            }
            // Set first superdiagonal
            else if (j==i+1){
                cserr = cs_cl_entry(Dxe, i, j, 1/(k0*DX));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dxe.\n"); exit(EXIT_FAILURE);}
            }
            // Set second band for Dye
            else if (j==i+NX){
                cserr = cs_cl_entry(Dye, i, j, 1/(k0*DY));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
            }
            // Impose Floquet boundary condition
            else if (j==i-NX*(NY-1)) {
                cserr = cs_cl_entry(Dye, i, j, cexpf(I*ky*Wy)/(k0*DY));
                if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
            }
        }
    } 
    */
    for (int i=0; i<Nc; i++){
        // Set diagonal elements
        cserr = cs_cl_entry(Dxe, i, i, -1/(k0*DX));
        if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dxe.\n"); exit(EXIT_FAILURE);}
        cserr = cs_cl_entry(Dye, i, i, -1/(k0*DY));
        if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}

        // Set first superdiagonal
        if (i < Nc-1){
            cserr = cs_cl_entry(Dxe, i, i+1, 1/(k0*DX));
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dxe.\n"); exit(EXIT_FAILURE);}
        }
        // Set second band for Dye
        if (i < Nc-NX){
            cserr = cs_cl_entry(Dye, i, i+NX, 1/(k0*DY));
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
        }
        // Impose Floquet boundary condition
        if (i >= NX*(NY-1)) {
            cserr = cs_cl_entry(Dye, i, i-NX*(NY-1), cexpf(I*ky*Wy)/(k0*DY));
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Dye.\n"); exit(EXIT_FAILURE);}
        }
    }
    Dxe = cs_cl_compress(Dxe);
    Dye = cs_cl_compress(Dye);
    printf("Finished defining finite-difference matrices.\n");
}

// --- Define full linear operator ---
void defineAe(double k0, double ky, double Wy) {
    int cserr;

    // Create local matrices for pieces of product
    cs_cl *muxxmat = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);
    cs_cl *muyymat = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);
    cs_cl *epszzmat = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);
    // Define diagonal entries with mu vectors
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            cserr = cs_cl_entry(muxxmat, i+j*NX, i+j*NX, 1/muxx[i+j*NX]);
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for muxx.\n"); exit(EXIT_FAILURE);}
            cserr = cs_cl_entry(muyymat, i+j*NX, i+j*NX, 1/muyy[i+j*NX]);
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for muyy.\n"); exit(EXIT_FAILURE);}
            cserr = cs_cl_entry(epszzmat, i+j*NX, i+j*NX, epszz[i+j*NX]);
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for epszz.\n"); exit(EXIT_FAILURE);}
        }
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
    cs_complex_t q[Nc];
    int cserr;
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            if (i<SFx) {
                q[i+j*NX] = 1.0;
            }
            else {
                q[i+j*NX] = 0.0;
            }
            cserr = cs_cl_entry(Q, i+j*NX, i+j*NX, q[i+j*NX]);
            if (cserr != 1) {fprintf(stderr, "Error encountered in cs_entry for Q.\n"); exit(EXIT_FAILURE);}
        }
    }

    cvecprint("data/q.csv", q, Nc);
    Q = cs_cl_compress(Q);

    #ifdef VERBOSE
        printf("===============\nQ\n");
        cs_cl_print(Q, 1);
    #endif
}


// --- Define nonlinear operator ---
cs_complex_t *nloperator(cs_complex_t *ez0){
    int SFx = LX + BX;
    static cs_complex_t ez1[NX*NY];
    for (int i=0; i<Nc; i++){
        ez1[i] = 0.0;
    }
    cs_cl_gaxpy(Ae, ez0, ez1);
    for (int i=0; i<Nc; i++){
        ez1[i] = ez1[i] - b[i];
        if ((i % NX) > SFx) {
            ez1[i] += 3/4 * chi3 * pow(cabs(ez0[i]),2)*ez0[i];
        }
    }
    return ez1;
}

// --- Define Anderson acceleration algorithm ---
/*
cs_complex_t *anderson(int m, cs_complex_t *ez0, float tol){
    // Adapted from pseudocode in https://doi.org/10.1016/j.cpc.2018.09.013
    int n = 0, k; // Iteration number
    float res = 1; // Residual
    cs_complex_t *ez1t, *ez1, *Gn1; // Updated guess

    ez1 = nloperator(ez0);
    Gn1 = cvecdiff(ez1, ez0, Nc);

    while (res >= tol){
        n += n;
        k = MIN(m, n);
        ez1 = nloperator(ez1);
    }
}
*/


// --- Define Newton algorithm ---
cs_complex_t *newtonsol(cs_complex_t *ez0, double tol){
    int n = 0; // Iteration number
    double res = 1; // Residual
    cs_complex_t *ez1; // Updated guess
    
    ez1 = nloperator(ez0);

    while (res >= tol){
        n += n;
        ez1 = nloperator(ez1);
    }
    return ez1;
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
    double omeg = 9.4248e+14; // source frequency
    double thetai = 0.0; // source angle
    double k0 = omeg / C0; // source wave number
    double kx = k0 * cos(thetai); 
    double ky = k0 * sin(thetai);

    // --- Define derived domain variables ---
    double Wy = NY * DY; // width of y domain (for Floquet)
    int SFx = LX + BX;

    FILE *fp;
    char fname[40];

    // --- Set the grid point values ---
    sprintf(fname, "data/X.csv");
    fp = fopen(fname, "w");
    for (int i=0; i<NX; i++){
        X[i] = i*DX;
        fprintf(fp, "%12.9e\n", X[i]);
    }
    fclose(fp);

    sprintf(fname, "data/Y.csv");
    fp = fopen(fname, "w");
    for (int j=0; j<NY; j++){
        Y[j] = j*DY;
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

    //cs_cl_print(Dxh, 1);

    // --- Define total linear operator ---
    printf("Defining linear operator Ae.\n");
    defineAe(k0, ky, Wy);

    // --- Define masking matrix ---
    printf("Defining masking matrix Q.\n");
    defineQ(SFx);

    // --- Use masking matrix to find source vector ---
    cs_cl *M = cs_cl_spalloc(Nc, Nc, 5*Nc, 1, 0);
    M = cs_cl_add(cs_cl_multiply(Q, Ae), cs_cl_multiply(Ae, Q), 1.0, -1.0);
    cs_cl_dropzeros(M);
    cs_cl_print(M, 1);
    cs_cl_gaxpy(M, fsrc, b);
    
    /* Simple bandaid
    cs_complex_t tmp;
    for (int j=0; j<NY; j++){
        tmp = b[(j+1)*NX-1];
        b[(j+1)*NX-1] = -1.0*creal(b[j*NX]) + I*cimag(b[j*NX]);
        b[j*NX] = -1.0*creal(tmp) + I*cimag(tmp);

        tmp = b[SFx-1 + j*NX];
        b[SFx-1 + j*NX] = -1.0*creal(b[SFx + j*NX]) + I*cimag(b[SFx + j*NX]);
        b[SFx + j*NX] = -1.0*creal(tmp) + I*cimag(tmp);
    }
    */

    printf("Source vector found. Printing.. \n");
    cvecprint("data/b.csv", b, Nc);
    printf("Finished printing. \n");

    #ifdef INCHECK
        printf("Input check only. Solver not invoked!\n");
        exit(EXIT_SUCCESS);
    #endif

    // --- Set tolerance and solve system with LU factorization ---
    double tol = 1.0;
    int cserr;
    printf("Solving system.\n");
    cserr = cs_cl_lusol(3, Ae, b, tol);
    if (cserr != 1) {fprintf(stderr, "Error encountered in solver.\n"); exit(EXIT_FAILURE);}
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

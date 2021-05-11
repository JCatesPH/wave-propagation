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

const int Nc = NX*NY;
cs_complex_t b[NX*NY], fsrc[NX*NY]; // Source vector
float X[NX]; // x grid points
float Y[NY]; // y grid points
float _Complex sx[LX];
cs_complex_t epszz[NX*NY], muxx[NX*NY], muyy[NX*NY]; 
cs_cl *Dxe, *Dye, *Dxh, *Dyh, *Ae, *Q; 


// --- Define source function ---
void fsource(float kx, float ky){
    printf("Defining source.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            fsrc[i*NY+j] = cexpf(I*(kx*X[i] + ky*Y[j]));
        }
    }
    printf("Finished defining source.\n");
}

// --- Define permittivity function ---
void epszzfunc(){
    float sigmax = -(sx_p + 1) * log(sx_R) / (2 * ETA0 * LX);
    float frac, s0x, sigx;

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
        }
    }
    printf("Finished defining permeability tensors.\n");
}

// --- Define finite-difference matrices ---
void finitediffs(float k0, float ky, float Wy) {

    // Need to set the diagonals for finite differences.
    printf("Defining finite-difference matrices.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            // Set diagonal elements
            if (i==j){
                cs_cl_entry(Dxe, i, j, -1/(k0*dx));
                cs_cl_entry(Dye, i, j, -1/(k0*dy));
            }
            // Set first superdiagonal
            else if (j==i+1){
                cs_cl_entry(Dxe, i, j, 1/(k0*dx));
            }
            // Set second band for Dye
            else if (j==i+NX){
                cs_cl_entry(Dye, i, j, 1/(k0*dx));
            }
            // Impose Floquet boundary condition
            else if (j==i-NX*(NY-1)) {
                cs_cl_entry(Dye, i, j, cexpf(I*ky*Wy));
            }
        }
    }
    Dxe = cs_cl_compress(Dxe);
    Dye = cs_cl_compress(Dye);
    printf("Finished defining finite-difference matrices.\n");
}

// --- Define full linear operator ---
void defineAe() {
    // Create local matrices for pieces of product
    cs_cl *muxxmat = cs_cl_spalloc(Nc, Nc, Nc, 0, 0);
    cs_cl *muyymat = cs_cl_spalloc(Nc, Nc, Nc, 0, 0);
    cs_cl *epszzmat = cs_cl_spalloc(Nc, Nc, Nc, 0, 0);
    // Define diagonal entries with mu vectors
    for (int i=0; i<Nc; i++){
        cs_cl_entry(muxxmat, i, i, 1/muxx[i]);
        cs_cl_entry(muyymat, i, i, 1/muyy[i]);
        cs_cl_entry(epszzmat, i, i, epszz[i]);
    }
    // Compress format for multiplication
    muxxmat = cs_cl_compress(muxxmat);
    muyymat = cs_cl_compress(muyymat);

    // Set vector to be all ones
    Ae = cs_cl_multiply(Dxh, muyymat);
    Ae = cs_cl_multiply(Ae, Dxe);
    Ae = cs_cl_add(Ae, cs_cl_multiply(Dyh, cs_cl_multiply(muxxmat, Dye)), 1.0, 1.0);
    Ae = cs_cl_add(Ae, epszzmat, 1.0, 1.0);

    // Clean up everything
    cs_cl_spfree(Dxh);
    cs_cl_spfree(Dxe);
    cs_cl_spfree(Dyh);
    cs_cl_spfree(Dye);
    cs_cl_spfree(epszzmat);
    cs_cl_spfree(muxxmat);
    cs_cl_spfree(muyymat);
}


// --- Define masking matrix for SF/TF ---
void defineQ(int SFx){
    for (int i=0; i<Nc; i++) {
        if ((i % NY) < SFx) {
            cs_cl_entry(Q, i, i, 1.0+I*0.0);
        }
        else {
            cs_cl_entry(Q, i, i, 0.0+I*0.0);
        }
    }
}

// ------- MAIN -------
int main ()
{
    // --- Allocate memory for variables ---
    Dxe = cs_cl_spalloc(Nc, Nc, 3*Nc, 0, 0);
    Dye = cs_cl_spalloc(Nc, Nc, 3*Nc, 0, 0);
    Ae = cs_cl_spalloc(Nc, Nc, 3*Nc, 0, -1);
    Q = cs_cl_spalloc(Nc, Nc, Nc, 0, 0);

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
    for (int i=0; i<Nc; i++) {
        fprintf(fp, "%12.9e, %+12.9e\n", creal(epszz[i]), cimag(epszz[i]));
    }


    // --- Generate and check the source field ---
    fsource(kx, ky);
    sprintf(fname, "data/src.csv");
    fp = fopen(fname, "w");
    for (int i=0; i<Nc; i++) {
        fprintf(fp, "%12.9e, %+12.9e\n", creal(fsrc[i]), cimag(fsrc[i]));
    }

    // --- Generate the finite-difference operators ---
    finitediffs(k0, ky, Wy);
    Dxh = cs_cl_transpose(Dxe, 0);
    Dyh = cs_cl_transpose(Dye, 0);


    // --- Define total linear operator ---
    defineAe();

    // --- Define masking matrix ---
    defineQ(SFx);


    // --- Use masking matrix to find source vector ---
    cs_cl_gaxpy(cs_cl_add(cs_cl_multiply(Q, Ae), cs_cl_multiply(Q, Ae), 1.0, -1.0), fsrc, b);
    printf("Source vector found.\n");
    sprintf(fname, "data/b.csv");
    fp = fopen(fname, "w");
    for (int i=0; i<Nc; i++) {
        fprintf(fp, "%12.9e, %+12.9e\n", creal(b[i]), cimag(b[i]));
    }

    // --- Set tolerance and solve system with LU factorization ---
    double tol = 1e-10;
    cs_cl_lusol(0, Ae, b, tol);
    printf("System solved.\n");

    // --- Print output for further plotting/analysis ---
    sprintf(fname, "data/Ez.csv");
    fp = fopen(fname, "w");
    for (int i=0; i<Nc; i++) {
        fprintf(fp, "%12.9e, %+12.9e\n", creal(b[i]), cimag(b[i]));
    }

    cs_cl_spfree(Ae);
    //return 0;
    printf("Exiting.\n\n");
}
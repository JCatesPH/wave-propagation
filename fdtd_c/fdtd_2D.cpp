/* From the Sullivan FDTD book */
#include "include/parameters.h"
#include "include/helpers.h"

// Define global variables
float gi2[NX], gi3[NX];                     // PML parameters (see Sullivan, Ch. 3)
float gj2[NY], gj3[NY];
float fi1[NX], fi2[NX], fi3[NX];
float fj1[NY], fj2[NY], fj3[NY];
float ga[NX][NY];                           // Dz to Ez  
float Ihx[NX][NY], Ihy[NX][NY];
float Dz[NX][NY], Hx[NX][NY], Hy[NX][NY];   // Field arrays
float Ez[NX][NY];                           // Actual field in z-direction   

int main()
{
    /* ---- Define variables ---- */
    int n,i,j; // indices
    FILE *fp;
    char fname[40];

    /* - Integration parameters - */
    int nsteps[NOUT]; // Number of time steps
    float dt = MIN(DX,DY)/6e8; // Time step (dx/2*c) 
    float T = 0.0; // Initial time

    for (n=0; n<NOUT; n++){
        nsteps[n] = n*DOUT;
    }

    /* ---- Print information ---- */
    printf("FDTD parameters: \n");
    printf("  Number of time steps : ");
    for (n=0; n<NOUT; n++){
        printf("%4d ", nsteps[n]);
    }
    printf("\n");
    printf("  Cell size : %.3f \n", DX);
    printf("  Time step : %.3g \n\n", dt);

    /* - Save x-y grid - */
    savegrid();
    /* - Set PML variables - */
    pmldef();

    /* ---- Initialize arrays ---- */
    for (int j=0; j<NY; j++){
        for (int i=0; i<NX; i++){
            Dz[i][j] = 0.0;
            Hx[i][j] = 0.0;
            Hy[i][j] = 0.0;
            Ihx[i][j] = 0.0;
            Ihy[i][j] = 0.0;
            ga[i][j] = 1.0 / epsz;
        }
    }

    /* ---- Main FDTD Loop ---- */
    printf("=============================\n");
    printf("Starting FDTD simulation. \n");
    printf("=============================\n");
    printf("  n  |  nsteps(n) \n");
    for (n=0; n<NOUT; n++){

        printf("%3d  |  %5d\n", n, nsteps[n]);
        T = nfdtdsteps(nsteps[n], T, dt);
        
        sprintf(fname, "data/Ez%04d.csv", nsteps[n]);
        fp = fopen(fname, "w");
        for (j=0; j<NY; j++){
            for (i=0; i<NX-1; i++){
                fprintf(fp, "%12.9f, ", Ez[i][j]);
            }
            fprintf(fp, "%12.9f\n", Ez[NX-1][j]);
        }
        fclose(fp);
        /* --- End of main loop --- */
    }
    printf("=============================\n");
    printf("Simulation finished successfully. \n");
    return 0;
}

/* From the Sullivan FDTD book */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define IE 200       // Number of cells in x
#define JE 200      // Number of cells in y
#define NOUT 20      // Number of time steps where output is stored
#define PI 3.14159
#define EPS0 8.85e-12

/* ---- Define some variables as global ---- */
float ga[IE][JE], dz[IE][JE], ez[IE][JE], hx[IE][JE], hy[IE][JE];
float ez_inc[JE], hx_inc[JE];
float ez_inc_low_m1, ez_inc_low_m2, ez_inc_high_m1, ez_inc_high_m2;

float epszz[IE][JE];
/* - Source parameters - */
const float t0 = 20.0; // Time of pulse peak
const float taup = 8.0; // Pulse duration
const int ic = IE/2 - 5; // First index of source location
const int jc = 12; // Second index of source location

/* Source function */
int src(float t) {
    ez_inc[jc] = exp(-4*log(2)*pow((t-t0)/taup, 2.0));
    return 0;
}

/* Permittivity function */
int perm() {
    for (int j=0; j<JE; j++){
        for (int i=0; i<IE; i++){
            epszz[i][j] = 1.0; 
        }
    }
    return 0;
}

int main()
{
    /* ---- Define variables ---- */
    int n,i,j,k; // indices
    float xn,xxn,xnum,curl_e;
    float gi2[IE], gi3[IE];
    float gj2[JE], gj3[JE];
    float fi1[IE], fi2[IE], fi3[IE];
    float fj1[JE], fj2[JE], fj3[JE];
    float ihx[IE][JE], ihy[IE][JE];
    float emin, emax; // Just for plotting purposes later
    emin = 0.0; emax = 0.0;
    FILE *fp;
    char fname[40];

    /* -- Integration parameters -- */
    int nsteps[NOUT]; // Number of time steps
    int npml = 10; // Number of cells in PML
    float ddx = 0.01; // Cell size
    float dt = ddx/6e8; // Time step (dx/2*c) 
    float T = 0.0; // Initial time

    /* - Total/scattered field - */ 
    int ia = 15;
    int ib = IE - ia - 1;
    int ja = 15;
    int jb = JE - ja - 1;

    float epsz = 1.0; // Permittivity

    for (n=0; n<NOUT; n++){
        nsteps[n] = n*2;
    }

    /* ---- Print information ---- */
    printf("FDTD parameters: \n");
    printf("  Number of time steps : ");
    for (n=0; n<NOUT; n++){
        printf("%4d ", nsteps[n]);
    }
    printf("\n");
    printf("  Cell size : %.3f \n", ddx);
    printf("  Time step : %.3g \n\n", dt);

    /* - Save x-y grid - */
    fp = fopen("data/X.csv", "w");
    for (i=0; i<IE-1; i++){
        fprintf(fp, "%12.9f, ", ddx*i);
    }
    fprintf(fp, "%12.9f", ddx*(IE-1));
    fclose(fp);
    fp = fopen("data/Y.csv", "w");
    for (j=0; j<JE-1; j++){
        fprintf(fp, "%12.9f, ", ddx*j);
    }
    fprintf(fp, "%12.9f", ddx*(JE-1));
    fclose(fp); 

    /* ---- Initialize arrays ---- */
    for (j=0; j<JE; j++){
        for (i=0; i<IE; i++){
            dz[i][j] = 0.0;
            hx[i][j] = 0.0;
            hy[i][j] = 0.0;
            ihx[i][j] = 0.0;
            ihy[i][j] = 0.0;
            ga[i][j] = 1.0 / epsz;
        }
        ez_inc[j] = 0.0;
        hx_inc[j] = 0.0;
    }
    ez_inc_low_m1 = 0; ez_inc_low_m2 = 0; 
    ez_inc_high_m1 = 0; ez_inc_high_m2 = 0;

    /* ---- Calculate the PML parameters ---- */
    for (i=0; i<IE; i++){
        gi2[i] = 1.0;
        gi3[i] = 1.0;
        fi1[i] = 0.0;
        fi2[i] = 1.0;
        fi3[i] = 1.0;
    }
    for (j=0; j<JE; j++){
        gj2[j] = 1.0;
        gj3[j] = 1.0;
        fj1[j] = 0.0;
        fj2[j] = 1.0;
        fj3[j] = 1.0;
    }

    for (i=0; i<=npml; i++){
        xnum = npml - i;
        xxn = xnum / npml;
        xn = 0.33*pow(xxn,3.0);
        gi2[i] = 1.0 / (1.0 + xn);
        gi2[IE-1-i] = 1.0 / (1.0 + xn);
        gi3[i] = (1.0 - xn) / (1.0 + xn);
        gi3[IE-1-i] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / npml;
        xn = 0.33*pow(xxn, 3.0);
        fi1[i] = xn;
        fi1[IE-2-i] = xn;
        fi2[i] = 1.0 / (1.0 + xn);
        fi2[IE-2-i] = 1.0 / (1.0 + xn);
        fi3[i] = (1.0 - xn) / (1.0 + xn);
        fi3[IE-2-i] = (1.0 - xn) / (1.0 + xn);
    }
    for (j=0; j<=npml; j++){
        xnum = npml - j;
        xxn = xnum / npml;
        xn = 0.33*pow(xxn,3.0);
        gj2[j] = 1.0 / (1.0 + xn);
        gj2[JE-1-j] = 1.0 / (1.0 + xn);
        gj3[j] = (1.0 - xn) / (1.0 + xn);
        gj3[JE-1-j] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / npml;
        xn = 0.33*pow(xxn, 3.0);
        fj1[j] = xn;
        fj1[JE-2-j] = xn;
        fj2[j] = 1.0 / (1.0 + xn);
        fj2[JE-2-j] = 1.0 / (1.0 + xn);
        fj3[j] = (1.0 - xn) / (1.0 + xn);
        fj3[JE-2-j] = (1.0 - xn) / (1.0 + xn);
    }

    /* ---- Main FDTD Loop ---- */
    printf("=============================\n");
    printf("Starting FDTD simulation. \n");
    printf("=============================\n");
    printf("  n  |  nsteps(n) \n");
    for (n=0; n<NOUT; n++){

        printf("%3d  |  %5d\n", n, nsteps[n]);

        for (k=1; k<=nsteps[n]; k++){
            T = T + 1;

            for (j=1; j<JE; j++){
                ez_inc[j] = ez_inc[j] + 0.5*(hx_inc[j-1]-hx_inc[j]);
            }

            /* ABC for the incident buffer */
            ez_inc[0] = ez_inc_low_m2;
            ez_inc_low_m2 = ez_inc_low_m1;
            ez_inc_low_m1 = ez_inc[1];

            ez_inc[JE-1] = ez_inc_high_m2;
            ez_inc_high_m2 = ez_inc_high_m1;
            ez_inc_high_m1 = ez_inc[JE-2];


            /* - Update flux density */
            for (j=1; j<JE; j++){
                for (i=1; i<IE; i++){
                    dz[i][j] = gi3[i] * gj3[j] * dz[i][j] 
                        + 0.5 * gi2[i] * gj2[j] * (hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);
                }
            }

            /* - Source - */
            src(T);
            /* - Incident Dz values - */ 
            //      Corrects TF/SF boundary derivatives
            for (i=ia; i<=ib; i++){
                dz[i][ja] = dz[i][ja] + 0.5*hx_inc[ja-1];
                dz[i][jb] = dz[i][jb] - 0.5*hx_inc[jb];
            }

            /* - Calculate Ez - */
            for (j=0; j<JE; j++){
                for (i=0; i<IE; i++){
                    ez[i][j] = ga[i][j] * dz[i][j];
                    if (ez[i][j] > emax){
                        emax = ez[i][j];
                    }
                    if (ez[i][j] < emin){
                        emin = ez[i][j];
                    }
                }
            }

            /* - Set Ez edges to 0, as part of the PML - */
            for (j=0; j<JE-1; j++){
                ez[0][j] = 0.0;
                ez[IE-1][j] = 0.0;
            }
            for (i=0; i<IE-1; i++){
                ez[i][0] = 0.0;
                ez[i][JE-1] = 0.0;
            }

            /* - Calculate incident Hx - */
            for (j=0; j<JE; j++){
                hx_inc[j] = hx_inc[j] + 0.5 * (ez_inc[j] - ez_inc[j+1]);
            }

            /* - Calculate Hx - */
            for (j=0; j<JE-1; j++){
                for (i=0; i<IE; i++){
                    curl_e = ez[i][j] - ez[i][j+1];
                    ihx[i][j] = ihx[i][j] + fi1[i] * curl_e;
                    hx[i][j] = fj3[j] * hx[i][j] + 0.5 * fj2[j] * (curl_e + ihx[i][j]);
                }
            }

            /* - Incident Hx values - */ 
            //      Corrects TF/SF boundary derivatives
            for (i=ia; i<=ib; i++){
                hx[i][ja-1] = hx[i][ja-1] + 0.5*ez_inc[ja];
                hx[i][jb] = hx[i][jb] - 0.5*ez_inc[jb];
            }

            /* - Calculate Hy - */
            for (j=0; j<JE-1; j++){
                for (i=0; i<IE-1; i++){
                    curl_e = ez[i+1][j] - ez[i][j];
                    ihy[i][j] = ihy[i][j] + fj1[j] * curl_e;
                    hy[i][j] = fi3[i] * hy[i][j] + 0.5 * fi2[i] * (curl_e + ihy[i][j]);
                }
            }
            /* - Incident Hy values - */ 
            //      Corrects TF/SF boundary derivatives
            for (j=ja; j<=jb; j++){
                hy[ia-1][j] = hy[ia-1][j] - 0.5*ez_inc[j];
                hy[ib][j] = hy[ib][j] + 0.5*ez_inc[j];
            }

        }
        
        sprintf(fname, "data/Ez%04d.csv", nsteps[n]);
        fp = fopen(fname, "w");
        for (j=0; j<JE; j++){
            for (i=0; i<IE-1; i++){
                fprintf(fp, "%12.9f, ", ez[i][j]);
            }
            fprintf(fp, "%12.9f\n", ez[IE-1][j]);
        }
        fclose(fp);
        /* --- End of main loop --- */
    }
    printf("=============================\n");
    printf("Simulation finished successfully. \n");
    printf("  min(Ez) = %6f \n  max(Ez) = %6f \n\n", emin, emax);
    return 0;
}
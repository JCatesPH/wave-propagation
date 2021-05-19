#include "include/helpers.h"


/* Source function */
float src(float t) {
    return sqrt(2*I0/(C0*EPS0)) * (sqrt(1-R)*exp(-2*log(2)*pow((t-t0)/taup, 2.0)) * cos(omeg0*(t-t0)) 
        + sqrt(R) * exp(-8*log(2)*pow((t-t0)/taup, 2.0)) * cos(2*omeg0*(t-t0) + phi) );
}


void pmldef() {
    float xn,xxn,xnum;

    /* ---- Calculate the PML parameters ---- */
    #pragma omp parallel for
    for (int i=0; i<NX; i++){
        gi2[i] = 1.0;
        gi3[i] = 1.0;
        fi1[i] = 0.0;
        fi2[i] = 1.0;
        fi3[i] = 1.0;
    }

    #pragma omp parallel for
    for (int j=0; j<NY; j++){
        gj2[j] = 1.0;
        gj3[j] = 1.0;
        fj1[j] = 0.0;
        fj2[j] = 1.0;
        fj3[j] = 1.0;
    }

    #pragma omp parallel for
    for (int i=0; i<=LX; i++){
        xnum = LX - i;
        xxn = xnum / LX;
        xn = 0.33*pow(xxn,3.0);
        gi2[i] = 1.0 / (1.0 + xn);
        gi2[NX-1-i] = 1.0 / (1.0 + xn);
        gi3[i] = (1.0 - xn) / (1.0 + xn);
        gi3[NX-1-i] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / LX;
        xn = 0.33*pow(xxn, 3.0);
        fi1[i] = xn;
        fi1[NX-2-i] = xn;
        fi2[i] = 1.0 / (1.0 + xn);
        fi2[NX-2-i] = 1.0 / (1.0 + xn);
        fi3[i] = (1.0 - xn) / (1.0 + xn);
        fi3[NX-2-i] = (1.0 - xn) / (1.0 + xn);
    }

    #pragma omp parallel for
    for (int j=0; j<=LY; j++){
        xnum = LY - j;
        xxn = xnum / LY;
        xn = 0.33*pow(xxn,3.0);
        gj2[j] = 1.0 / (1.0 + xn);
        gj2[NY-1-j] = 1.0 / (1.0 + xn);
        gj3[j] = (1.0 - xn) / (1.0 + xn);
        gj3[NY-1-j] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / LY;
        xn = 0.33*pow(xxn, 3.0);
        fj1[j] = xn;
        fj1[NY-2-j] = xn;
        fj2[j] = 1.0 / (1.0 + xn);
        fj2[NY-2-j] = 1.0 / (1.0 + xn);
        fj3[j] = (1.0 - xn) / (1.0 + xn);
        fj3[NY-2-j] = (1.0 - xn) / (1.0 + xn);
    }
}


void savegrid() {
    FILE *fp;

    fp = fopen("data/X.csv", "w");
    for (int i=0; i<NX-1; i++){
        fprintf(fp, "%12.9f, ", DX*i);
    }
    fprintf(fp, "%12.9f", DX*(NX-1));
    fclose(fp);

    fp = fopen("data/Y.csv", "w");
    for (int j=0; j<NY-1; j++){
        fprintf(fp, "%12.9f, ", DY*j);
    }
    fprintf(fp, "%12.9f", DY*(NY-1));
    fclose(fp); 

}

float nfdtdsteps(int nsteps, float T, float dt) {
    float curl_e;

    for (int k=1; k<=nsteps; k++){
        
        T = T + dt;
        /* - Update flux density */
        #pragma omp parallel for
        for (int j=1; j<NY; j++){
            for (int i=1; i<NX; i++){
                Dz[i][j] = gi3[i] * gj3[j] * Dz[i][j] 
                    + 0.5 * gi2[i] * gj2[j] * (Hy[i][j] - Hy[i-1][j] - Hx[i][j] + Hx[i][j-1]);
            }
        }

        /* - Source - */
        Dz[ic][jc] = src(T);

        /* - Calculate Ez - */
        #pragma omp parallel for
        for (int j=0; j<NY; j++){
            for (int i=0; i<NX; i++){
                Ez[i][j] = ga[i][j] * Dz[i][j];
            }
        }

        /* - Set Ez edges to 0, as part of the PML - */
        #pragma omp parallel for
        for (int j=0; j<NY-1; j++){
            Ez[0][j] = 0.0;
            Ez[NX-1][j] = 0.0;
        }
        #pragma omp parallel for
        for (int i=0; i<NX-1; i++){
            Ez[i][0] = 0.0;
            Ez[i][NY-1] = 0.0;
        }

        /* - Calculate Hx - */
        #pragma omp parallel for
        for (int j=0; j<NY-1; j++){
            for (int i=0; i<NX; i++){
                curl_e = Ez[i][j] - Ez[i][j+1];
                Ihx[i][j] = Ihx[i][j] + fi1[i] * curl_e;
                Hx[i][j] = fj3[j] * Hx[i][j] + 0.5 * fj2[j] * (curl_e + Ihx[i][j]);
            }
        }

        /* - Calculate Hy - */
        #pragma omp parallel for
        for (int j=0; j<NY-1; j++){
            for (int i=0; i<NX-1; i++){
                curl_e = Ez[i+1][j] - Ez[i][j];
                Ihy[i][j] = Ihy[i][j] + fj1[j] * curl_e;
                Hy[i][j] = fi3[i] * Hy[i][j] + 0.5 * fi2[i] * (curl_e + Ihy[i][j]);
            }
        }

    }

    return T;
}

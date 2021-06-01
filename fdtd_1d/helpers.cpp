#include "include/helpers.h"

using namespace std;

/* Source function */
double src(double t) {
    return sqrt(2*I0/(C0*EPS0)) * (sqrt(1-R)*exp(-2*log(2)*pow((t-t0)/taup, 2.0)) * cos(omeg0*(t-t0)) 
        + sqrt(R) * exp(-8*log(2)*pow((t-t0)/taup, 2.0)) * cos(2*omeg0*(t-t0) + phi) );
}


void pmldef() {
    double xn,xxn,xnum;
    /* ---- Calculate the PML parameters ---- */
    #pragma omp parallel for
    for (int i=0; i < NX; i++){
        gi2[i] = 1.0;
        gi3[i] = 1.0;
        fi1[i] = 0.0;
        fi2[i] = 1.0;
        fi3[i] = 1.0;
    }

    for (int i=0; i <= LX; i++){
        xnum = LX - i;
        xxn = xnum / LX;
        xn = 0.33*pow(xxn,SX_P);
        gi2[i] = 1.0 / (1.0 + xn);
        gi2[NX-1-i] = 1.0 / (1.0 + xn);
        gi3[i] = (1.0 - xn) / (1.0 + xn);
        gi3[NX-1-i] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / LX;
        xn = 0.33*pow(xxn, SX_P);
        fi1[i] = xn;
        fi1[NX-2-i] = xn;
        fi2[i] = 1.0 / (1.0 + xn);
        fi2[NX-2-i] = 1.0 / (1.0 + xn);
        fi3[i] = (1.0 - xn) / (1.0 + xn);
        fi3[NX-2-i] = (1.0 - xn) / (1.0 + xn);
    }

}


void savegrid() {
    FILE *fp;
    fp = fopen("data/X.csv", "w");
    for (int i=0; i<NX-1; i++){
        fprintf(fp, "%12.9f\n", DX*i);
    }
    fprintf(fp, "%12.9f", DX*(NX-1));
    fclose(fp);

}

void initArrays(double dt) {
    ofstream freqfile;
    char valstr[20];

    #pragma omp parallel for
    for (int i=0; i<NX; i++){
        Dx[i] = 0.0;
        Hy[i] = 0.0;
        Ix[i] = 0.0;

        if (i > idie) {
            gb[i] = sigma * dt / epsz;
            ga[i] = 1.0 / (epsz + gb[i]);
        }
        else {
            gb[i] = 0.0;
            ga[i] = 1.0;
        }
        
        for (int k=0; k < NF; k++) {
            Ew_re[k][i] = 0.0;
            Ew_im[k][i] = 0.0;
        }
    }

    //freqfile.open("data/omeg.csv");
    for (int k=0; k < NF; k++) {
        Fsrc_re[k] = 0.0;
        Fsrc_im[k] = 0.0;
        omeg[k] = 2*PI/((NOUT*DOUT-1)*dt) * k;
        //sprintf(valstr, "%12.9f", omeg[k]);
        //freqfile << valstr << endl;
    }
    //freqfile.close();
}

double nfdtdsteps(int N, double T, double dt) {
    double cr2 = alpha * chi3 * ETA0 * ETA0;
    double r2;
    for (int n=1; n<=N; n++){
        T = T + dt;
        /* - Update flux density */
        #pragma omp parallel for
        for (int i=1; i<NX; i++){
            Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]);
            //Dx[i] = Dx[i] + 0.5 * (Hy[i-1] - Hy[i]);
        }

        /* - Source - */
        // Dx[isrc] = Dx[isrc] + src(T); // soft source
        Dx[isrc] += src(T); // hard source

        /* - Calculate Ez - */
        #pragma omp parallel for
        for (int i=0; i < NX-1; i++){
            // Ex[i] = ga[i] * (Dx[i] - Ix[i]); // original
            r2 = cr2 * pow(Ex[i], 2);
            Ex[i] = (Dx[i] + 2 * r2) / (1.0/ga[i] + 3 * r2 * Ex[i]); // Kerr effect included
            Ix[i] = Ix[i] + gb[i] * fi1[i] * Ex[i];
            //Ix[i] = Ix[i] + gb[i] * Ex[i];
        }
        /* - Calculate Fourier transform - */
        #pragma omp parallel for
        for (int k=0; k < NF; k++) {
            for (int i=1; i<NX; i++){
                Ew_re[k][i] += cos(omeg[k]*T) * Ex[i];
                Ew_im[k][i] += sin(omeg[k]*T) * Ex[i];
            }
            if (T < 2*t0) {
                Fsrc_re[k] += cos(omeg[k]*T) * Ex[isrc];
                Fsrc_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
        }
        
        /* - Set Ez edges to 0, as part of the PML - */
        Ex[0] = 0.0;
        Ex[NX-1] = 0.0;

        /* - Calculate Hx - */
        #pragma omp parallel for
        for (int i=0; i<NX; i++){
            Hy[i] = fi3[i] * Hy[i] + 0.5 * fi2[i] * (Ex[i] - Ex[i+1]);
            //Hy[i] = Hy[i] + 0.5 * (Ex[i] - Ex[i+1]);
        }

    }

    return T;
}

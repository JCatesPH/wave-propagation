/* From the Sullivan FDTD book */
#include "include/parameters.h"
#include "include/helpers.h"

using namespace std;

// Define global variables
double gi2[NX], gi3[NX];                     // PML parameters (see Sullivan, Ch. 3)
double fi1[NX], fi2[NX], fi3[NX];            // PML parameters (see Sullivan, Ch. 3)
double ga[NX], gb[NX];                               // Dz to Ez  
double Ix[NX];
double Dx[NX], Hy[NX];           // Field arrays
double Ex[NX];                           // Actual field in z-direction   
double chi3;
double Ew_re[NF][NX], Ew_im[NF][NX]; // Real and imaginary frequency domain output.
double Fsrc_re[NF], Fsrc_im[NF];     // Real and imaginary frequency domain source.
double omeg[NF];                     // Frequencies for DFT.
double amp_in[NF];

int main()
{
    auto t1 = Clock::now();
    /* ---- Define variables ---- */
    int n,i,j; // indices
    ofstream outfile;
    string fname;

    /* - Integration parameters - */
    int nsteps[NOUT]; // Number of time steps
    double dt = DX / 6e8; // Time step (dx/2*c) 
    double T = 0.0; // Initial time

    for (n=0; n<NOUT; n++){
        nsteps[n] = n*DOUT;
    }

    chi3 = 4 * EPS0 * C0 * n2 / 3;

    /* ---- Print information ---- */
    cout << "FDTD parameters: " << endl;
    cout << "  Number of time steps : ";
    for (n=0; n<NOUT; n++){
        printf("%4d ", nsteps[n]);
    }
    cout << endl << "-----------------------------" << endl;
    cout << " Grid parameters:" << endl;
    cout << "   Num. grid     : " << NX << endl;
    cout << "   Cells in PML  : " << LX << endl;
    cout << "   Buffer region : " << BX << endl;
    cout << "   Grid spacing  : " << DX << endl;
    cout << "   Time Step     : " << dt << endl;
    cout << "   Total time    : " << NOUT*DOUT*dt << endl;
    cout << "-----------------------------" << endl;
    cout << " PML parameters:" << endl;
    cout << "   Polynomial    : " << SX_P << endl;
    cout << "   Reflectivity  : " << SX_R << endl;
    cout << "   Smax          : " << SX_M << endl;
    cout << "-----------------------------" << endl;
    cout << " Beam parameters: (two-color Gaussian pulse)" << endl;
    cout << "   Source index  : " << isrc << endl;
    cout << "   Intensity     : " << I0 << endl;
    cout << "   Pulse peak    : " << t0 << endl;
    cout << "   FWHM dur.     : " << taup << endl;
    cout << "   Frequency     : " << omeg0 << endl;
    cout << "   2nd Ratio     : " << R << endl;
    cout << "   Phase offset  : " << phi << endl;
    cout << "-----------------------------" << endl;
    cout << " Material parameters:" << endl;
    cout << "   Rel. Perm.    : " << epsz << endl;
    cout << "   Conductivity  : " << sigma << endl;
    cout << "   Nonlinear ind.: " << n2 << endl;
    cout << "   Chi3          : " << chi3 << endl;



    /* - Save x-y grid - */
    savegrid();
    /* - Set PML variables - */
    pmldef();

    /* ---- Initialize arrays ---- */
    initArrays(dt);

    char valstr[50];
    outfile.open("data/omeg.csv");
    for (int k=0; k < NF; k++) {
        sprintf(valstr, "%10.7e", omeg[k]);
        outfile << valstr << endl;
    }
    outfile.close();

    /* ---- Main FDTD Loop ---- */
    printf("=============================\n");
    printf("Starting FDTD simulation. \n");
    printf("=============================\n");
    printf("  n  |  nsteps(n)  |  T  \n");
    printf("-----------------------------\n");

    outfile.open("data/Ex.csv");
    for (n=0; n<NOUT; n++){

        printf("%3d  |  %5d  | %8.5e\n", n, nsteps[n], T);
        T = nfdtdsteps(DOUT, T, dt);
        
        for (i=0; i<NX-1; i++){
            sprintf(valstr, "%12.9f, ", Ex[i]);
            outfile << valstr;
        }
        sprintf(valstr, "%12.9f", Ex[NX-1]);
        outfile << valstr << endl;

    }
    outfile.close();
    printf("=============================\n");
    /* --- End of main loop --- */

    for (int k=0; k < NF; k++) {
        amp_in[k] = sqrt(pow(Fsrc_re[k],2) + pow(Fsrc_im[k],2));
    }

    ofstream outfile1, outfile2;
    outfile1.open("data/Ew_re.csv");
    outfile2.open("data/Ew_im.csv");
    for (int i=0; i<NX; i++){
        for (int k=0; k < NF-1; k++) {
            sprintf(valstr, "%12.9e, ", Ew_re[k][i] / amp_in[k]);
            outfile1 << valstr;
            sprintf(valstr, "%12.9e, ", Ew_im[k][i] / amp_in[k]);
            outfile2 << valstr;
        }
        sprintf(valstr, "%12.9e", Ew_re[NF-1][i] / amp_in[NF-1]);
        outfile1 << valstr << endl;
        sprintf(valstr, "%12.9e", Ew_im[NF-1][i] / amp_in[NF-1]);
        outfile2 << valstr << endl;
    }
    outfile1.close();
    outfile2.close();


    auto t2 = Clock::now();
    printf("Simulation finished successfully. \n");
    cout << "  Time to calculate: " 
        <<  (double) std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6
        << " [seconds]" << std::endl;
    return 0;
}

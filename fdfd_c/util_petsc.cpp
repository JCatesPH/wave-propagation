#include "include/util_petsc.h"

using namespace std;

Mat Dxe, Dye, Dxh, Dyh;
complex<double> sx[LX];

// --- Define source function ---
int fsource(double A0, double kx, double ky){
    PetscErrorCode ierr;
    FILE *fp;
    char fname[40];
    double X[NX]; // x grid points
    double Y[NY]; // y grid points

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

    // --- Define source on grid ---
    PetscPrintf(MPI_COMM_WORLD,"Defining source.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++){
            //fsrc[i+j*NX] = cexp(I*(kx*X[i] + ky*Y[j]));
            ierr = VecSetValue(fsrc, i+j*NX, A0*exp(-1i*(kx*X[i] + ky*Y[j])), INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    PetscPrintf(MPI_COMM_WORLD,"Finished defining source.\n");
    VecAssemblyBegin(fsrc);
    VecAssemblyEnd(fsrc);
    return ierr;
}

// --- Define permittivity function ---
int epszzfunc(double epsr){
    PetscErrorCode ierr;
    double sigmax = -(SX_P + 1) * log(SX_R) / (2 * ETA0 * LX);
    #ifdef VERBOSE
        PetscPrintf(MPI_COMM_WORLD,"sigmax = %f\n", sigmax);
    #endif
    double frac, s0x, sigx;

    // Define PML region and change eps accordingly.
    PetscPrintf(MPI_COMM_WORLD,"Defining permittivity tensor.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            //printf("%d, %d\n", i, j);
            if (i < LX) {
                frac = (double) (LX-i)/LX;
                s0x = 1.0 + SX_M * pow(frac, SX_P);
                sigx = sigmax * pow(sin(PI*frac/2.0), 2.0);
                sx[i] = s0x * (1.0 - 1i * ETA0 * sigx);
                // epszz[i+j*NX] = sx[i];
                ierr = VecSetValue(epszz, i+j*NX, sx[i], INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (i > NX-LX) {
                ierr = VecSetValue(epszz, i+j*NX, sx[NX-1-i], INSERT_VALUES); CHKERRQ(ierr);
            }
            // --- Sets the rel perm in slab ---
            else if (i > idx1 && i < idx2) {
                ierr = VecSetValue(epszz, i+j*NX, epsr, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }

    // --- Index matching ---
    double r1 = 1;
    for (int i = idx1 - 10; i < idx1 + 10; i++){
        for (int j=0; j<NY; j++) {
            // Sigmoid function: S[x_] := (eps2 - eps1)/(1 + Exp[-5 x]) + eps1
            ierr = VecSetValue(epszz, i+j*NX, 
                (epsr - 1) / ( 1 + exp(-r1 * (i - idx1)) ) + 1, 
                INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    for (int i = idx2 - 10; i < idx2 + 10; i++){
        for (int j=0; j<NY; j++) {
            // Sigmoid function: S[x_] := (eps2 - eps1)/(1 + Exp[-5 x]) + eps1
            ierr = VecSetValue(epszz, i+j*NX, 
                (epsr - 1) / ( 1 + exp(r1 * (i - idx2)) ) + 1, 
                INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    VecAssemblyBegin(epszz);
    VecAssemblyEnd(epszz);
    PetscPrintf(MPI_COMM_WORLD,"Finished defining permittivity tensor.\n");
    return ierr;
}


// --- Define permeability tensors ---
int mufunc() {
    PetscErrorCode ierr;
    // Define PML region and change mu accordingly.
    PetscPrintf(MPI_COMM_WORLD,"Defining permeability tensors.\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            /* 
                ! NOTE : Defined as 1/muxx and 1/muyy for simplicity later.
            */ 
            if (i < LX) {
                ierr = VecSetValue(muxx, i+j*NX,   sx[i], INSERT_VALUES); CHKERRQ(ierr);
                ierr = VecSetValue(muyy, i+j*NX, 1/sx[i], INSERT_VALUES); CHKERRQ(ierr);
            }
            else if (i > NX-LX) {
                ierr = VecSetValue(muxx, i+j*NX,   sx[NX-1-i], INSERT_VALUES); CHKERRQ(ierr);
                ierr = VecSetValue(muyy, i+j*NX, 1/sx[NX-1-i], INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }
    VecAssemblyBegin(muxx);
    VecAssemblyEnd(muxx);
    VecAssemblyBegin(muyy);
    VecAssemblyEnd(muyy);
    PetscPrintf(MPI_COMM_WORLD,"Finished defining permeability tensors.\n");
    return ierr;
}


// --- Define finite-difference matrices ---
int finitediffs(double k0, double ky, double Wy) {
    PetscErrorCode ierr;

    // Need to set the diagonals for finite differences.
    PetscPrintf(MPI_COMM_WORLD,"Defining finite-difference matrices.\n");
    PetscScalar tmp0x[1], tmp0y[1], tmp1[1], tmpNx[1], tmpNxNy1[1];
    tmp0x[0] = -1/(k0*DX); tmp0y[0] = -1/(k0*DY); tmp1[0] = 1/(k0*DX); tmpNx[0] = 1/(k0*DY); tmpNxNy1[0] = exp(1i*ky*Wy)/(k0*DY);
    PetscInt i0[1], i1[1], iNx[1], iNxNy1[1];

    for (int i=0; i<Nc; i++){
        // Set diagonal elements
        i0[0] = i; i1[0] = i+1; iNx[0] = i+NX; iNxNy1[0] = i-NX*(NY-1);
        ierr = MatSetValues(Dxe, 1, i0, 1, i0, tmp0x, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(Dye, 1, i0, 1, i0, tmp0y, INSERT_VALUES); CHKERRQ(ierr);

        // Set first superdiagonal
        if (i < Nc-1){
            ierr = MatSetValues(Dxe, 1, i0, 1, i1, tmp1, INSERT_VALUES); CHKERRQ(ierr);
        }
        // Set second band for Dye
        if (i < Nc-NX){
            ierr = MatSetValues(Dye, 1, i0, 1, iNx, tmpNx, INSERT_VALUES); CHKERRQ(ierr);
        }
        // Impose Floquet boundary condition
        if (i >= NX*(NY-1)) {
            ierr = MatSetValues(Dye, 1, i0, 1, iNxNy1, tmpNxNy1, INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    MatAssemblyBegin(Dxe, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Dxe, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Dye, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Dye, MAT_FINAL_ASSEMBLY);

    PetscPrintf(MPI_COMM_WORLD,"Finished defining finite-difference matrices.\n");
    return ierr;
}


// --- Define masking matrix for SF/TF ---
int defineQ(int SFx){
    PetscErrorCode ierr;
    Vec q;
    ierr = VecCreate(PETSC_COMM_WORLD,&q); CHKERRQ(ierr);
    ierr = VecSetSizes(q,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(q);CHKERRQ(ierr);

    ierr = VecSet(q, 0.0+1i*0.0);CHKERRQ(ierr);
    for (int i=0; i<NX; i++){
        for (int j=0; j<NY; j++) {
            if (i < SFx) {
                ierr = VecSetValue(q, i+j*NX, 1.0+1i*0.0, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }
    VecAssemblyBegin(q);
    VecAssemblyEnd(q);

    PetscViewer lab;
    PetscViewerCreate(PETSC_COMM_WORLD, &lab);
    PetscViewerSetType(lab, PETSCVIEWERASCII);
    PetscViewerFileSetMode(lab,FILE_MODE_WRITE);
    PetscViewerFileSetName(lab,"data/q.csv");
    VecView(q, lab);
    PetscViewerDestroy(&lab);

    ierr = MatDiagonalSet(Q, q, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    return ierr;
}


// --- Define full linear operator ---
int defineAe(Mat *Ap, double k0, double ky, double Wy) {
    //Mat A = *Ap;
    PetscErrorCode ierr;
    // --- Vector of ones for 'MatDiagonalScale' calls ---
    Vec Ones;
    ierr = VecCreate(PETSC_COMM_WORLD,&Ones); CHKERRQ(ierr);
    ierr = VecSetSizes(Ones,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(Ones);CHKERRQ(ierr);

    ierr = VecSet(Ones, 1.0+1i*0.0);CHKERRQ(ierr);

    // --- Declare finite-difference operators ---
    
    ierr = MatCreate(PETSC_COMM_WORLD, &Dxe); CHKERRQ(ierr);
    ierr = MatSetSizes(Dxe,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Dxe); CHKERRQ(ierr);
    ierr = MatSetUp(Dxe); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Dxh); CHKERRQ(ierr);
    ierr = MatSetSizes(Dxh,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Dxh); CHKERRQ(ierr);
    ierr = MatSetUp(Dxh); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Dye); CHKERRQ(ierr);
    ierr = MatSetSizes(Dye,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Dye); CHKERRQ(ierr);
    ierr = MatSetUp(Dye); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Dyh); CHKERRQ(ierr);
    ierr = MatSetSizes(Dyh,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Dyh); CHKERRQ(ierr);
    ierr = MatSetUp(Dyh); CHKERRQ(ierr);

    // 
    Mat K;
    ierr = MatCreate(PETSC_COMM_WORLD, &K); CHKERRQ(ierr);
    ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(K); CHKERRQ(ierr);
    ierr = MatSetUp(K); CHKERRQ(ierr);

    // --- Generate the finite-difference operators ---
    finitediffs(k0, ky, Wy);
    ierr = MatHermitianTranspose(Dxe, MAT_INITIAL_MATRIX, &Dxh); CHKERRQ(ierr);
    ierr = MatHermitianTranspose(Dye, MAT_INITIAL_MATRIX, &Dyh); CHKERRQ(ierr);
    ierr = MatScale(Dxh, -1.0+1i*0.0); CHKERRQ(ierr);
    ierr = MatScale(Dyh, -1.0+1i*0.0); CHKERRQ(ierr);

    // --- Carry out multiplications and additions ---
    ierr = MatDiagonalScale(Dxe, muyy, Ones); CHKERRQ(ierr);
    ierr = MatMatMult(Dxh, Dxe, MAT_INITIAL_MATRIX, PETSC_DEFAULT, Ap); CHKERRQ(ierr);
    ierr = MatDiagonalScale(Dye, muxx, Ones); CHKERRQ(ierr);
    ierr = MatMatMult(Dyh, Dye, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &K); CHKERRQ(ierr);
    ierr = MatAXPY(*Ap, 1.0, K, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatDiagonalSet(*Ap, epszz, ADD_VALUES); CHKERRQ(ierr);

    // Clean up everything
    ierr = VecDestroy(&Ones); CHKERRQ(ierr);
    ierr = MatDestroy(&K); CHKERRQ(ierr);

    ierr = MatDestroy(&Dxe); CHKERRQ(ierr);
    ierr = MatDestroy(&Dye); CHKERRQ(ierr);
    ierr = MatDestroy(&Dxh); CHKERRQ(ierr);
    ierr = MatDestroy(&Dyh); CHKERRQ(ierr);

    PetscPrintf(MPI_COMM_WORLD,"Finished defining Ae and freed memory.\n");
    return ierr;
}


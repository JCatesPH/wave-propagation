#include "include/util_petsc.h"

using namespace std;

//#define VERBOSE // Uncomment to print additional debugging messages.
//#define INCHECK // Uncomment to exit before solver is invoked. For checking inputs.

const int Nc = NX*NY;
Vec b_SFTF, fsrc; // Source vector
Vec b_NL1, b_NL2; // Nonlinear source terms
Vec ez_11, ez_21; // Nonlinear source terms from last iteration
Vec ez_1, ez_1cc, ez_2; // Electric field for fundamental and second harmonic
Vec res1, res2; // Residuals for first and second harmonic 
Vec epszz, muxx, muyy; 
Mat Ae1, Ae2; // Linear operators for first and second harmonic
Mat Q; // Masking matrix for SF/TF

static char help[] = "[Empty help message]\n\n";


// ------- MAIN -------
int main (int argc, char **args)
{
    auto t1 = Clock::now();
    KSP            ksp1, ksp2;          /* linear solver context */
    PetscErrorCode ierr;
    PetscMPIInt size;

    // --- Initialize PETSc ---
    ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRMPI(ierr);
    //if (size != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");

    // --- Create 'lab' for IO ---
    PetscViewer lab;
    PetscViewerCreate(PETSC_COMM_WORLD, &lab);
    PetscViewerSetType(lab, PETSCVIEWERASCII);
    PetscViewerFileSetMode(lab,FILE_MODE_WRITE);

    // --- Create vectors ---
    ierr = VecCreate(PETSC_COMM_WORLD,&b_SFTF); CHKERRQ(ierr);
    ierr = VecSetSizes(b_SFTF,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b_SFTF);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&b_NL1); CHKERRQ(ierr);
    ierr = VecSetSizes(b_NL1,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b_NL1);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&b_NL2); CHKERRQ(ierr);
    ierr = VecSetSizes(b_NL2,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b_NL2);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&fsrc); CHKERRQ(ierr);
    ierr = VecSetSizes(fsrc,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(fsrc);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&ez_1); CHKERRQ(ierr);
    ierr = VecSetSizes(ez_1,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ez_1);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&ez_11); CHKERRQ(ierr);
    ierr = VecSetSizes(ez_11,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ez_11);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&ez_1cc); CHKERRQ(ierr);
    ierr = VecSetSizes(ez_1cc,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ez_1cc);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&res1); CHKERRQ(ierr);
    ierr = VecSetSizes(res1,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(res1);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&ez_2); CHKERRQ(ierr);
    ierr = VecSetSizes(ez_2,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ez_2);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&res2); CHKERRQ(ierr);
    ierr = VecSetSizes(res2,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(res2);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&ez_21); CHKERRQ(ierr);
    ierr = VecSetSizes(ez_21,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ez_21);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&epszz); CHKERRQ(ierr);
    ierr = VecSetSizes(epszz,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(epszz);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&muxx); CHKERRQ(ierr);
    ierr = VecSetSizes(muxx,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(muxx);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&muyy); CHKERRQ(ierr);
    ierr = VecSetSizes(muyy,PETSC_DECIDE,Nc); CHKERRQ(ierr);
    ierr = VecSetFromOptions(muyy);CHKERRQ(ierr);

    ierr = VecSet(epszz, 1.0);CHKERRQ(ierr);
    ierr = VecSet(muxx, 1.0);CHKERRQ(ierr);
    ierr = VecSet(muyy, 1.0);CHKERRQ(ierr);

    // --- Allocate memory for variables ---
    ierr = MatCreate(PETSC_COMM_WORLD, &Ae1); CHKERRQ(ierr);
    ierr = MatSetSizes(Ae1,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Ae1); CHKERRQ(ierr);
    ierr = MatSetUp(Ae1); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Ae2); CHKERRQ(ierr);
    ierr = MatSetSizes(Ae2,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Ae2); CHKERRQ(ierr);
    ierr = MatSetUp(Ae2); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Q); CHKERRQ(ierr);
    ierr = MatSetSizes(Q,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Q); CHKERRQ(ierr);
    ierr = MatSetUp(Q); CHKERRQ(ierr);


    //Dxe = cs_cl_spalloc(Nc, Nc, 3*Nc, 1, 1);
    //Dye = cs_cl_spalloc(Nc, Nc, 3*Nc, 1, 1);
    //Ae = cs_cl_spalloc(Nc, Nc, 5*Nc, 1, 0);
    //Q = cs_cl_spalloc(Nc, Nc, Nc, 1, 1);

    // --- Find wavenumbers ---
    double k0 = omeg / C0; // source wave number
    double kx = k0 * cos(thetai); 
    double ky = k0 * sin(thetai);

    // --- Define derived domain variables ---
    double Wy = NY * DY; // width of y domain (for Floquet)
    int SFx = LX + BX;

    // --- Define and save the permittivity and permeability tensors ---
    epszzfunc(epsr1);
    mufunc();

    PetscViewerFileSetName(lab,"data/epszz.csv");
    VecView(epszz, lab);

    // --- Generate and check the source field ---
    fsource(A0, kx, ky);

    PetscViewerFileSetName(lab,"data/src.csv");
    VecView(fsrc, lab);

    // --- Define total linear operator ---
    PetscPrintf(MPI_COMM_WORLD,"Defining linear operator Ae1.\n");
    defineAe(&Ae1, k0, ky, Wy);
    
    PetscPrintf(MPI_COMM_WORLD,"Defining linear operator Ae2.\n");
    epszzfunc(epsr2);
    defineAe(&Ae2, 2*k0, ky, Wy);

    // --- Assemble operators ---
    ierr = MatAssemblyBegin(Ae1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Ae1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    ierr = MatAssemblyBegin(Ae2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Ae2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    // --- Destroy perm tensors ---
    ierr = VecDestroy(&epszz); CHKERRQ(ierr);
    ierr = VecDestroy(&muxx); CHKERRQ(ierr);
    ierr = VecDestroy(&muyy); CHKERRQ(ierr);

    // --- Define masking matrix ---
    PetscPrintf(MPI_COMM_WORLD,"Defining masking matrix Q.\n");
    defineQ(SFx);

    // --- Use masking matrix to find source vector ---
    Mat M1, M2;
    ierr = MatCreate(PETSC_COMM_WORLD, &M1); CHKERRQ(ierr);
    ierr = MatSetSizes(M1,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(M1); CHKERRQ(ierr);
    ierr = MatSetUp(M1); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &M2); CHKERRQ(ierr);
    ierr = MatSetSizes(M2,PETSC_DECIDE,PETSC_DECIDE,Nc,Nc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(M2); CHKERRQ(ierr);
    ierr = MatSetUp(M2); CHKERRQ(ierr);

    ierr = MatMatMult(Q, Ae1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &M1); CHKERRQ(ierr);
    ierr = MatMatMult(Ae1, Q, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &M2); CHKERRQ(ierr);
    ierr = MatAXPY(M1, -1.0, M2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = MatMult(M1, fsrc, b_SFTF); CHKERRQ(ierr);

    ierr = MatDestroy(&M1); CHKERRQ(ierr);
    ierr = MatDestroy(&M2); CHKERRQ(ierr);


    PetscPrintf(MPI_COMM_WORLD,"Source vector found. Printing.. \n");
    PetscViewerFileSetName(lab,"data/b.csv");
    VecView(b_SFTF, lab);
    PetscPrintf(MPI_COMM_WORLD,"Finished printing. \n");

    #ifdef INCHECK
        PetscPrintf(MPI_COMM_WORLD,"Input check only. Solver not invoked!\n");
        exit(EXIT_SUCCESS);
    #endif

    // --- Set tolerance and initialize solvers ---
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp1);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp1, Ae1, Ae1);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp1,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp1);CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp2); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp2, Ae2, Ae2); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp2,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp2); CHKERRQ(ierr);

    PetscPrintf(MPI_COMM_WORLD,"Solving linear system.. \n\n");

    ierr = KSPSolve(ksp1, b_SFTF, ez_1); CHKERRQ(ierr);
    ierr = KSPView(ksp1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    PetscPrintf(MPI_COMM_WORLD,"\n\nSystem solved. Finding generated harmonic field.\n");

    // --- Computes b_{NL,2\omega} = -d_{33} e_{z,\omega} \times e_{z,\omega}
    VecPointwiseMult(b_NL2, ez_1, ez_1);

    PetscViewerFileSetName(lab,"data/Ez_om.csv");
    VecView(ez_1, lab);

    VecScale(b_NL2, -d33);

    VecAssemblyBegin(b_NL2);
    VecAssemblyEnd(b_NL2);

    PetscPrintf(MPI_COMM_WORLD,"Saving nonlinear source vector, b_NL2.\n");
    
    PetscViewerFileSetName(lab,"data/b_NL2.csv");
    VecView(b_NL2, lab);

    // --- Solves for second-harmonic field ---
    PetscPrintf(MPI_COMM_WORLD,"\nSolving second-harmonic system.\n");
    ierr = KSPSolve(ksp2, b_NL2, ez_2); CHKERRQ(ierr);
    ierr = KSPView(ksp2, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    PetscPrintf(MPI_COMM_WORLD,"\n\nSystem solved. Finding b_NL1.\n");
    VecCopy(ez_1, ez_1cc);
    VecConjugate(ez_1cc);
    VecPointwiseMult(b_NL1, ez_2, ez_1cc);
    VecScale(b_NL1, -2*d33);

    VecAssemblyBegin(b_NL1);
    VecAssemblyEnd(b_NL1);

    PetscPrintf(MPI_COMM_WORLD,"Saving nonlinear source vector, b_NL1.\n");
    PetscViewerFileSetName(lab,"data/b_NL1.csv");
    VecView(b_NL1, lab);

    // --- Refine solution and check norm ---
    VecCopy(ez_1, ez_11);
    VecAXPY(b_NL1, 1.0, b_SFTF);
    VecAssemblyBegin(b_NL1);
    VecAssemblyEnd(b_NL1);
    ierr = KSPSolve(ksp1, b_NL1, ez_1); CHKERRQ(ierr);
    
    VecAXPY(ez_11, -1.0, ez_1);

    double resnorm;
    ierr = VecNorm(ez_11, NORM_INFINITY, &resnorm);
    PetscPrintf(PETSC_COMM_WORLD,"L_inf Norm of residual vector: %g\n", resnorm);
    ierr = VecNorm(ez_11, NORM_2, &resnorm);
    PetscPrintf(PETSC_COMM_WORLD,"L_2 Norm of residual vector: %g\n", resnorm);


    // --- Print output for further plotting/analysis ---
    PetscViewerFileSetName(lab,"data/Ez_om_2.csv");
    VecView(ez_1, lab);

    PetscViewerFileSetName(lab,"data/Ez_2om.csv");
    VecView(ez_2, lab);

    // --- Free remaining memory ---
    ierr = MatDestroy(&Ae1); CHKERRQ(ierr);
    ierr = MatDestroy(&Ae2); CHKERRQ(ierr);
    ierr = MatDestroy(&Q); CHKERRQ(ierr);
    ierr = VecDestroy(&b_SFTF); CHKERRQ(ierr);
    ierr = VecDestroy(&b_NL1); CHKERRQ(ierr);
    ierr = VecDestroy(&b_NL2); CHKERRQ(ierr);
    ierr = VecDestroy(&ez_1); CHKERRQ(ierr);
    ierr = VecDestroy(&ez_2); CHKERRQ(ierr);

    ierr = KSPDestroy(&ksp1); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp2); CHKERRQ(ierr);

    PetscPrintf(MPI_COMM_WORLD,"Simulation finished successfully. \n");

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
        auto t2 = Clock::now();
        std::cout << "  Time to calculate: " 
        <<  (double) std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6
        << " [seconds]" << std::endl
        <<  (double) std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6 / 60
        << " [minutes]" << std::endl;;
    }

    ierr = PetscFinalize();
    return ierr;
}

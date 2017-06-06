/* Solves the linear sistem of equations*/

#include "fermi.h"

int fersolv_ST(void)
{
    int    i;
    double error;

    EPSSetOperators(eps,A,B);
    EPSSolve(eps);
    EPSGetIterationNumber(eps,&its);
    EPSGetDimensions(eps,&nev,NULL,NULL);

    for(i=0;i<nev;i++)
    {
      EPSGetEigenpair(eps,i,&ikeff,NULL,phi_n,NULL);
      keff=1/ikeff;
      EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);
      PetscPrintf(FERMI_Comm,"Keff: %lf Its:%d error:%e\n",keff,its,error); 
    }

    return 0;
}

int fersolv_TR(void)
{
//    PetscViewerASCIIOpen(FERMI_Comm,"b.dat",&viewer);
//    VecView(b,viewer);
//    PetscViewerASCIIOpen(FERMI_Comm,"M.dat",&viewer);
//    MatView(M,viewer);
    KSPSetOperators(ksp,M,M);
    KSPSolve(ksp,b,phi_n);
    KSPGetIterationNumber(ksp,&its);
    return 0;
}

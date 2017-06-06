/* Rutines for assembly K matrix */

#include "fermi.h"

int ferass_TR(int step)
{
  int e;
  node_list_t *pr,*px,*pe;

  /* Matrix M is calculating only once */
  if(!step)
  {
    MatZeroEntries(M);
    for(e=0;e<mesh.nelemv;e++)
    {
      if(ferelem_M(e))
      {
        return 1;
      }
      MatSetValues(M,nke,idxm,nke,idxm,Me,ADD_VALUES);
    }/*elem loop*/
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  }

  /* Matrix A used for calculating RHS b = A*phi */
  MatZeroEntries(A);
  for(e=0;e<mesh.nelemv;e++)
  {
    if(ferelem_ABM(e))
    {
      return 1;
    }
    MatSetValues(A,nke,idxm,nke,idxm,Ae,ADD_VALUES);
    MatSetValues(A,nke,idxm,nke,idxm,Be,ADD_VALUES);
    MatSetValues(A,nke,idxm,nke,idxm,Me,ADD_VALUES);
  }/*elem loop*/
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  /*control rods pertubations*/
  pr=list_ctrlr.head;
  while(pr)
  {
    pe=((ctrlrod_t*)pr->data)->elemv.head;
    px=((ctrlrod_t*)pr->data)->xsa.head;
    while(pe)
    {
      if(ferelem_R(*(int*)pe->data,*(double*)px->data,-1.0))
      {
        return 1;
      }
      pe=pe->next;
      px=px->next;
      MatSetValues(A,nke,idxm,nke,idxm,Be,ADD_VALUES);
    }
    pr=pr->next;
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  
  /* b = A*phi_n */
  MatMult(A,phi_n,b);

  /* Dirichlet BC */
  VecSetValues(b,ndir,dirIndex,dirZeros,INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatSetOption(M,MAT_KEEP_NONZERO_PATTERN , PETSC_TRUE );
  MatZeroRows(M,ndir,dirIndex,1.0,NULL,NULL);
  return 0;
}

int ferass_TR_1(int step)
{
  /* Assemblies object for the TR simulation 
   * if step = 0 assemblies the mass matrix M 
   the final A = A + B/keff + M, 
   the perturbation matrix B 
   calculates b=(A+B)*phi_n
   assemblies precursors terms on b += be
   puts boundary conditions on M and b 

   if step > 0 assemblies the perturbation matrix B           
   calculates b=(A+B)*phi_n
   assemblies precursors terms on b += be
   puts boundary conditions on M and b 
   */
  int e;
  node_list_t *pr,*px,*pe;
  if(!step)
  {
    MatZeroEntries(M);
    for(e=0;e<mesh.nelemv;e++)
    {
      if(ferelem_M(e))
      {
        PetscPrintf(FERMI_Comm,"Problem calculating elemental objects TR steps.\n"); 
        return 1;
      }
      MatSetValues(M,nke,idxm,nke,idxm,Me,ADD_VALUES);
    }/*elem loop*/
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
    MatScale(A,-1.0);
    MatAXPY(A,ikeff,B,SAME_NONZERO_PATTERN);
    MatAXPY(A,1.0,M,SAME_NONZERO_PATTERN);
  }
  /*control rods pertubations*/
  MatZeroEntries(B);
  pr=list_ctrlr.head;
  while(pr)
  {
    pe=((ctrlrod_t*)pr->data)->elemv.head;
    px=((ctrlrod_t*)pr->data)->xsa.head;
    while(pe)
    {
      if(ferelem_R(*(int*)pe->data,*(double*)px->data,-1.0))
      {
        PetscPrintf(FERMI_Comm,"ferassm.c:problem calculating rod element perturbation.\n"); 
        return 1;
      }
      pe=pe->next;
      px=px->next;
      MatSetValues(B,nke,idxm,nke,idxm,Be,ADD_VALUES);
    }
    pr=pr->next;
  }
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  MatZeroEntries(K);
  MatAXPY(K,1.0,A,SAME_NONZERO_PATTERN);
//  MatAXPY(K,1.0,B,DIFFERENT_NONZERO_PATTERN);
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  MatZeroRows(K,ndir,dirIndex,1.0,NULL,NULL);
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
  /* b = A*phi_n */
  MatMult(K,phi_n,b);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  /* Neumann BC */

  /* Dirichlet BC */
  VecSetValues(b,ndir,dirIndex,dirZeros,INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatSetOption(M,MAT_KEEP_NONZERO_PATTERN , PETSC_TRUE );
  MatZeroRows(M,ndir,dirIndex,1.0,NULL,NULL);
  return 0;
}

int ferass_ST(void)
{

  /* Assemblies object for the ST simulation */
  int e;
  node_list_t *pr,*px,*pe;

  MatZeroEntries(A);
  MatZeroEntries(B);

  vol=0.0;

  for(e=0;e<mesh.nelemv;e++)
  {
    if(ferelem_AB(e)){
      PetscPrintf(FERMI_Comm,"Problem calculating elemental objects ST steps.\n"); 
      return 1;
    }
    MatSetValues(A,nke,idxm,nke,idxm,Ae,ADD_VALUES);
    MatSetValues(B,nke,idxm,nke,idxm,Be,ADD_VALUES);
  }/*elem loop*/

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

  /*control rods pertubations*/
  pr=list_ctrlr.head;
  while(pr)
  {
    pe=((ctrlrod_t*)pr->data)->elemv.head;
    px=((ctrlrod_t*)pr->data)->xsa.head;
    while(pe)
    {
      if(ferelem_R(*(int*)pe->data,*(double*)px->data,+1.0))
      {
        return 1;
      }
      pe=pe->next;
      px=px->next;
      MatSetValues(A,nke,idxm,nke,idxm,Be,ADD_VALUES);
    }
    pr=pr->next;
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  /* Dirichlet BC */
  MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
  MatZeroRows(A,ndir,dirIndex,1.0,NULL,NULL);

  MatSetOption(B,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
  MatZeroRows(B,ndir,dirIndex,0.0,NULL,NULL);

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  //    PetscViewerASCIIOpen(FERMI_Comm,"A.dat",&viewer);
  //    MatView(A,viewer);
  //    PetscViewerASCIIOpen(FERMI_Comm,"B.dat",&viewer);
  //    MatView(B,viewer);
  return 0;
}

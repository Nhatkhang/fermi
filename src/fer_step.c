/* 

   FERMI basic evolutions kinds

   ST: Static

   TR: Transient

*/

#include "fermi.h"

int ferstep_ST(void)
{

  if(fersrods(calcu.t)){
    PetscPrintf(FERMI_Comm,"Problem calculating control rods.\n"); 
    return 1;
  }
  if(ferass_ST()){
    PetscPrintf(FERMI_Comm,"Problem assembling the A & B."); 
    return 1;
  }
  if(fersolv_ST()){
    PetscPrintf(FERMI_Comm,"Problem solving the steady state.\n"); 
    return 1;
  }

  /* Normalize the flux according to the power */
  if(fer_norm()){
    PetscPrintf(FERMI_Comm,"Problem normalizing the power.\n"); 
    return 1;
  }

  return 0;
}

int ferstep_TR(int step)
{

  if(fersrods(calcu.t)){
    PetscPrintf(FERMI_Comm,"Problem calculating control rods.\n"); 
    return 1;
  }
  if(ferass_TR(step)){
    PetscPrintf(FERMI_Comm,"Problem assembling the A & b.\n"); 
    return 1;
  }
  if(fersolv_TR()){
    PetscPrintf(FERMI_Comm,"Problem solving the transient steps.\n"); 
    return 1;
  }

  return 0;
}

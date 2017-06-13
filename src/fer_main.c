/* 

FERMI main program

*/


#include "fermi.h"

int main(int argc,char **argv)
{

  int           step;
  char          nam[32];
  node_list_t * pNod;

  if(ferinit(argc,argv))
    goto error;

  calcu.t = calcu.t0;
  step=0;
  
  if(calcu.timedep == QS){

    pNod=calcu.time.head;
    while(pNod){

      dtn=((tcontrol_t*)pNod->data)->dt;

      // recv information
      fer_comm_step(COUPLE_RECV);

      // Assembly the system and solve with SLEPc
      if(ferstep_ST())
        goto error;

      sprintf(nam,"steady_r%d_t%d",rank,step);
      print_vtk(nam);
      print_out(&phi_n, step);

      // send information
      fer_comm_step(COUPLE_SEND);

      calcu.t=calcu.t + dtn;
      step ++;
      if(calcu.t>((tcontrol_t*)pNod->data)->tf-1.0e-8)
        pNod=pNod->next;

    }

  }

  else if(calcu.timedep == TR){

    //================================================== 
    // Transient simulation
    //
    // 1) Calculate steady state solving Ax = (1/k) Bx
    // 2) Recv information for coupling (if it required)
    // 3) Calculates a "dt" increase in the flux solving
    //    Ax = b
    // 4) Send information for coupling (if it required)
    // 5) Repeat from "2)" up to achieving final time "tf" 
    //
    PetscPrintf(FERMI_Comm,"calculating stationary state.\n");
    
    // Assembly the system and solve with SLEPc
    if(ferstep_ST())
      goto error;

    sprintf(nam,"steady_r%d",rank);
    print_vtk(nam);

    PetscPrintf(FERMI_Comm,"initial power:%e\n",power);
    PetscPrintf(FERMI_Comm,"time    power   its\n");

    pNod=calcu.time.head;
    while(pNod){

      dtn=((tcontrol_t*)pNod->data)->dt;

      // recv information
      fer_comm_step(COUPLE_RECV);

      // Assembly the system and solve with PETSc
      if(ferstep_TR(step))
        goto error;

      // Calculates power
      if(fer_pow(&power))
        goto error;

      PetscPrintf(FERMI_Comm,"%lf %e %d \n",calcu.t,power,its);

      sprintf(nam,"tran_rank%d_t%d",rank,step);
      print_vtk(nam);
      print_out(&phi_n, step);

      // send information
      fer_comm_step(COUPLE_SEND);

      calcu.t=calcu.t + dtn;
      step ++;
      if(calcu.t>((tcontrol_t*)pNod->data)->tf-1.0e-8)
        pNod=pNod->next;

    }
  }

error:

  ferfini();  

  return 0;
}


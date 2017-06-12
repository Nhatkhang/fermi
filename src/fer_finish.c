/* frees all the memory */

#include "fermi.h"

int ferfini(void)
{

  //================================================================================ 
  // CLOSING OUTPUT FILES (in those who need it)
  //================================================================================ 
  //    
  node_list_t  * pn;
  output_t    * po;
  pn = list_outpu.head;
  while(pn){
    po = (output_t *)pn->data;
    switch (po->kind){

      case 1:
	break;

      case 2:
	// power on physical entities as a function of time 
	fclose(po->kind_2.fp);
	break;

      default:
	return 1;

    }
    pn = pn->next;
  }
  //
  //================================================================================ 

  // eliminamos las estructuras de comunicadores
  // e intercomunicadores
  //    commdom_delete(); 

  MPI_Barrier(WORLD_Comm);
  SlepcFinalize();
  MPI_Finalize();

  return 0;
}

/* 

Communications routines for FERMI using PLEPP library

Authors: 
Federico Caccia  (fedecaccia32@gmail.com)
Guido Giuntoli   (guido.giuntoli@bsc.es)

-> FERMI can be connected to an "indefinite" number of process.

*/

#include <mpi.h>
#include "fermi.h"


int fer_comm_step(int order)
{

  node_list_t *pn;
  comm_t      *comm;
  int         count;
  int         tag;
  MPI_Status  status;

  // We travel all the communications in list_comm
  // and perform the communication with all of them
  while(pn)
  {

    // First and input order is received,
    // these order should be syncronize with
    // the execution of external codes

    comm = (comm_t*)pn->data;

    switch(comm->kind){

      case 1:

	if(order == COUPLE_RECV){

	  count = egn * nxs_mat;

	// we receive cross sections
	MPI_recv(
	    comm.comm_1.xs, 
	    count, 
	    MPI_DOUBLE,
	    comm.comm_1.rem_leader,
	    tag,
	    comm.comm_1.intercomm,
	    &status;
	    );
	}
	else if(order == COUPLE_SEND){

	}

	break;

      default:
	return 1;
    }
  
    pn = pn->next;

  }

         
  /**************************************************/

  return 0;
}


/**************************************************/

int fer_comm_init(MPI_Comm *world_comm, 
    MPI_Comm *FERMI_Comm, 
    MPI_Comm *INTER_Comm)
{

  /*
     We initialize the local communicator FERMI_Comm and 
     the inter-communicator array INTER_Comm[]

  */

  #ifdef COMMDOM 

      int   i;
      int   ierr;
      char  my_name[] = "fermi"; // name for PLEPP coupling scheme

      commdom_create();
      commdom_set_names(coupling.world, my_name);
      commdom_create_commij(world_comm, FERMI_Comm);

      // We travel all the friend of FERMI 
      // in order to get the inter-communicators
      // created by commdom_create_commij
      for(i=0; i < coupling.num_friends; i++){
	commdom_get_commij(coupling.friends[i],&INTER_Comm[i]);
      }

  #endif

  return 0;

}

/**************************************************/

int fer_corecv(MPI_Comm * couple_comm)
{
    int        ierr;
    MPI_Status status;

    // Receive cross sections data from all materials specified on Input
    ierr = MPI_Recv(&order, 1, MPI_INTEGER, 0, tag, couple_comm, &status)

    return 0;
} 

/**************************************************/

int fer_cosend(MPI_Comm * couple_comm, int * control_fg)
{
    // Sends data calculates and the waits for the server order for control flow

    int ierr, tag = 0;

    ierr = MPI_Send (&output_var, N_output_var,MPI_DOUBLE, 0, tag, couple_comm);

    // Receive control_fg instruction
    ierr = MPI_Recv(control_fg, 1, MPI_INTEGER, 0, tag, couple_comm, &status)

    return 0;
} 

/****************************************************************************************************/

int fer_coends(MPI_Comm * couple_comm)
{
    // Finish the communication with the server
    MPI_Comm_disconnect(couple_comm);

    return 0;
} 

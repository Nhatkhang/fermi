/* 

Communications routines for FERMI using PLEPP library

Authors: 
Federico Caccia  (fedecaccia32@gmail.com)
Guido Giuntoli   (guido.giuntoli@bsc.es)

-> FERMI can be connected to an "indefinite" number of process.

*/

#include "mpi.h"
#include "fermi.h"


int fer_comm_step(int order)
{

  node_list_t * pn, * pm;
  comm_t      * comm;
  int           count;
  int           tag=0;
  int           ierr;
  int           i, en;
  MPI_Status    status;
  pvl_t       * mat;

  // We travel all the communications in list_comm
  // and perform the communication with all of them
  pn = list_comms.head;
  while(pn)
  {

    // First and input order is received,
    // these order should be syncronize with
    // the execution of external codes

    comm = (comm_t*)pn->data;

    switch(comm->kind){

      case 1:

	if(order == COUPLE_RECV){

	  count = egn * nxs_mat * comm->comm_1.nphy;

	  // we receive cross sections
	  ierr = MPI_Recv( comm->comm_1.xs,  count,  MPI_DOUBLE, comm->comm_1.rem_leader,  tag,  *(comm->comm_1.intercomm), &status );
	  if(ierr){
	    return 1;
	  }
           
	  // put cross sections on the list_mater structure
	  // TODO: the cross section structure should be change
	  // on the future because is not efficient to localize cross
	  // sections in a list in this kind of code (arrays are better)
	  for( i=0 ; i < comm->comm_1.nphy ; i++ ){
            

	    // localizamos el material en la lista y copiamos
	    pm = list_mater.head;
	    while(pm){
	      mat = (pvl_t*)pm->data;

	      if(strcmp(mat->name, comm->comm_1.phys[i]) == 0){
		
		// copiamos D
		for(en=0 ; en < egn; en++){
		  mat->D[en]      = comm->comm_1.xs[i * nxs_mat + 0*egn + en];
		}

		// copiamos xsa
		for(en=0 ; en < egn; en++){
		  mat->xs_a[en]   = comm->comm_1.xs[i * nxs_mat + 1*egn + en];
		}

		// copiamos nxs_f
		for(en=0 ; en < egn; en++){
		  mat->nxs_f[en]  = comm->comm_1.xs[i * nxs_mat + 2*egn + en];
		}

		// copiamos exs_f
		for(en=0 ; en < egn; en++){
		  mat->exs_f[en]  = comm->comm_1.xs[i * nxs_mat + 3*egn + en];
		}

		// copiamos chi
		for(en=0 ; en < egn; en++){
		  mat->chi[en]    = comm->comm_1.xs[i * nxs_mat + 4*egn + en];
		}

		// copiamos xs_s
		for(en=0 ; en < egn*(egn-1); en++){
		  mat->xs_s[en]   = comm->comm_1.xs[i * nxs_mat + 5*egn + en];
		}

	      }

	      pm = pm->next;
	    }
	  }


	}
	else if(order == COUPLE_SEND){

	  count = comm->comm_1.nphy;

          // calculate powers on each physical entity
	  fer_pow_phys( comm->comm_1.nphy, comm->comm_1.ids, comm->comm_1.pow );

	  // we receive cross sections
	  ierr = MPI_Send( comm->comm_1.pow,  count,  MPI_DOUBLE, comm->comm_1.rem_leader,  tag,  *(comm->comm_1.intercomm));
	  if(ierr){
	    return 1;
	  }

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

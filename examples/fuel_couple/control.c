/*

   This program is to test the coupling of FERMI with other
   codes. In this case this is the master code which sends the 
   orders to FERMI for calculating the steady flux under different 
   XS distributions and values.

   Authors: 

   Miguel Zavala
   Guido Giuntoli

*/

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "commdom_wrapper.h"  

int main(int argc,char **argv)
{

  int    i;
  int    rank       = -1;
  int    size       = -1;
  int    local_rank = -1;
  int    local_size = -1;
  int    size_commij;

  char   file_c[64];

  
  MPI_Comm world_comm;   // Global communicator
  MPI_Comm CONTROL_Comm; // Local communicator
  MPI_Comm INTER_Comm;   // Inter-communicator

  world_comm = MPI_COMM_WORLD;
  CONTROL_Comm = MPI_COMM_NULL;
  INTER_Comm   = MPI_COMM_NULL;

  MPI_Init(&argc, &argv);

  char world[]   = "acople";
  char my_name[] = "control";
  char friend[]  = "fermi";

  commdom_create();

  commdom_set_names( world, sizeof(world), my_name, sizeof(my_name))
  commdom_create_commij(&world_comm, &CONTROL_Comm);
  MPI_Comm_rank(CONTROL_Comm, local_rank);
  MPI_Comm_size(CONTROL_Comm, local_size);
  commdom_get_commij_size(size_commij);
  
  MPI_Finalize();

  return 0;
}

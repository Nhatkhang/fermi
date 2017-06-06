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
  int    globa_rank = -1;
  int    globa_size = -1;
  int    local_rank = -1;
  int    local_size = -1;
  int    sw, sm;
  int    size_commij;

  char   file_c[64];

  
  MPI_Comm WORLD_Comm;   // Global communicator
  MPI_Comm CONTROL_Comm; // Local communicator
  MPI_Comm INTER_Comm;   // Inter-communicator

  MPI_Init(&argc, &argv);

  WORLD_Comm   = MPI_COMM_WORLD;
  CONTROL_Comm = MPI_COMM_NULL;
  INTER_Comm   = MPI_COMM_NULL;

  MPI_Comm_rank(WORLD_Comm, &globa_rank);
  MPI_Comm_size(WORLD_Comm, &globa_size);
  printf("control.c : globa_rank = %d globa_size = %d\n", globa_rank, globa_size);

  char world[]   = "acople";
  char my_name[] = "control";
  char Friend[]  = "fermi";

  commdom_create();

  sw = sizeof(world);
  sm = sizeof(my_name);
  commdom_set_names( world, &sw, my_name, &sm);
  commdom_create_commij((int*)&WORLD_Comm, (int*)&CONTROL_Comm);
  MPI_Comm_rank(CONTROL_Comm, &local_rank);
  MPI_Comm_size(CONTROL_Comm, &local_size);
  printf("control.c : local_rank = %d local_size = %d\n", local_rank, local_size);
  commdom_get_commij_size(&size_commij);
  
  MPI_Finalize();

  return 0;
}

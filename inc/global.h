/*

   FERMI global variable declaration

*/

#include <petscksp.h>
#include <slepceps.h>
#include "stdbool.h"
#include "list.h"
#include "types.h"
#include "mesh.h"

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#define DIM 3
#define NPE 8

MPI_Comm   WORLD_Comm; // global communicator
MPI_Comm   FERMI_Comm; // local  communicator
MPI_Comm * INTER_Comm; // array of inter-communicators to communicate
                       // with other codes
int      * remote_ranks; // remotes ranks in INTER_Comm to communicate
                         // with the others

int       globa_rank; // rank in WORLD_Comm
int       globa_size; // size in WORLD_Comm
int       local_rank; // rank in FERMI_Comm
int       local_size; // size in FERMI_Comm

enum {QS, TR};        /*Quasi Static, transient*/
enum {K1};            /*Calculation of elemental matrix by this K modes*/

bool       couple_fl;
coupling_t coupling;

PetscViewer kspview;
PetscViewer viewer;

int        rank;
int        nproc;   
int        nummat;                 
             
int      * loc2gold, *loc2gnew;
int        ndir;
int      * npp;
int        ntot;
int      * dirIndex;
double   * dirValue;
double   * dirZeros;
                    
int memory;

list_t list_nodes;
list_t list_ghost;
list_t list_elemv;
list_t list_elems;
list_t list_physe;
list_t list_mater;
list_t list_bound;
list_t list_fun1d; /* list of functions */
list_t list_ctrlr; /* list of control rods */
list_t list_outpu;
list_t list_comms;

mesh_t mesh;

int nke,nbe,nev,its;
double *Ae,*Be,*Me,*be;
Vec  phi_n,phi_o,b,b_a;
Vec  xlocal;
double ikeff,keff;
Mat A,B,M,K;
KSP ksp;
PC  pc;          /* preconditioner context */
EPS eps;
int Istart,Iend;

char inputfile[32];
char meshfile[32];
char epartfile[32];
char npartfile[32];
                    
calcu_t  calcu;

int      *idxm;
double   **der,***ode,**sh,**jac,**ijac,**coor,*wp;

int      egn, pgn;
int      nxs_mat;      // number of xs values per material
double   **phi;
double   power;
double   vol;

/* Precursors constants */
double *veloc,*lambda,*beta,*chi,beta_tot;      
double dtn;
#endif

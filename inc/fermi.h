/* 

   Structure types and constants

 */

#ifndef _FUNCT_H_
#define _FUNCT_H_

#include <slepceps.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "global.h"
#include "mesh.h"
#include "list.h"
#include "types.h"
#include "gmsh.h"
#include "fem.h"
#include "fun.h"
#include "utils.h"

#ifdef COMMDOM
  #include "commdom_wrapper.h"
#endif

/*

   Constants definition

*/

// Coupling orders used for external communication

#define COUPLE_INIT  0    
#define COUPLE_RECV  1
#define COUPLE_SEND  2
#define COUPLE_ENDS  3


//parser.c
int parse_input(void);
int parse_mesh(void);
int parse_mats(void);
int parse_mode(void);
int parse_func(void);
int parse_boun(void);
int parse_crod(void);
int parse_outp(void);
int parse_communication(void);
int parse_material(char *buff, pvl_t *mat);
int cmp_mat(void *a, void *b);
int parse_boundary(char *bufcpy, bound_t *bou);
int cmp_bou(void *a, void *b);
int cmp_time(void *a, void *b);
int get_int(char *buf, const char *name,int *a);
int get_char(char *buf, const char *name,char *a);
int parse_coupling(const char file_c[]);

int ferass_TR(int step);
int ferass_ST(void);
int ferinit(int argc,char **argv);
int ferfini(void);
int ferstep_ST(void);
int ferstep_TR(int step);
int fersolv_ST(void);
int fersolv_TR(void);
int ferelem_AB(int e);
int ferelem_ABM(int e);
int ferelem_M(int e);
int ferelem_R(int e,double xsa,double factor);
//lst2msh.c
int cpynode (node_list_t *node_nl, node_t *node);
int cpyelemv (node_list_t *elem_nl, elem_t *elemv);
int cpyelems (node_list_t *elem_nl, elem_t *elems);
//ferboun.c
int ferbouset(void);
int cmp_nod(void *a, void *b);

//fer_power.c
int fer_norm(void);
int fer_pow(double *fpower);
int fer_pow_phys(int n, int * ids, double *fpower);
int fer_pow_elem(int e, double *fpower);

//ferrods.c
int ferirods(void);
int fersrods(double t);
//output.c
int print_struct(int step);
int print_vtk(const char *name);
int printMatrixR(char *name,double *A, int m, int n);
int printMatrixRC(char *name,double **A, int m, int n);
int printVector(char *name,double *vec, int m);
int printMat(char *name,Mat *A, int m);
int printVec(char *name,Vec *vec, int m);
int vtkcode(int dim,int npe);
int print_out(Vec *phi, int step);

//fer_couple.c
int fer_couple(int order, MPI_Comm * couple_comm, char * server_n);
int fer_coinit(MPI_Comm * couple_comm, char * server_n);
int fer_corecv(MPI_Comm * couple_comm);
int fer_cosend(MPI_Comm * couple_comm, int * control_fg);
int fer_coends(MPI_Comm * couple_comm);

int fer_comm_init(void);
int fer_comm_step(int order);

#endif

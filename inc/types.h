
#ifndef _TYPES_H_
#define _TYPES_H_

#include "list.h"
#include "mesh.h"
#include "fun.h"

enum {SEQUENCIAL, PARALLEL};

typedef struct bit{

  unsigned int val : 1;

}bit_t;

typedef struct _propES_t{

  int gmshid;  /* ID that correspond to evewry element in the gmshfile */
  int elemv;   /* The elemv local numeration at which belongs */

}ps_t;

typedef struct _bound_t{

  char name[16];
  int kind;
  int order;
  int dimS;    
  list_t nodeL;    
  list_t elemsL;    
  list_t elemvL;    

}bound_t;


typedef struct _pv_t{

  char     * name;
  int        gmshid;

  double   * D;               // diffusion coeficient
  double   * xs_a;            // absortion XS  
  double   * nxs_f;           // nu x fission XS
  double   * exs_f;           // energy x fission XS
  double   * xs_s;            // scattering XS 
  double   * xs_r;            // remotion XS
  double   * chi;             // fission spectrum

  int        hasprec;
  double   * conc;            // Fission precursors concentration  ( I groups ) 

}pv_t;

typedef struct _pvl_t{

  char       name[16];

  double   * D;               // diffusion coeficient
  double   * xs_a;            // absortion XS  
  double   * nxs_f;           // nu x fission XS
  double   * exs_f;           // energy x fission XS
  double   * xs_s;            // scattering XS 
  double   * xs_r;            // remotion XS
  double   * chi;             // fission spectrum

  int        hasprec;

}pvl_t;

/*************************************************************/
typedef struct _kind_1_t{

  char     phys[16];

}kind_1_t;

/* this is used to print localized powers on physical entities on a file */
typedef struct _kind_2_t{

  FILE     * fp;
  char       file[16];
  char    ** phys;     // array of Physical Entities names
  int        nphy;
  int      * ids ;
  double   * pow;

}kind_2_t;

typedef struct _output_t{

  /* esto va a volar pronto */
  char       file[16];
  char       phys[16];
  double     norm[3];

  int        kind;

  kind_1_t   kind_1;
  kind_2_t   kind_2;

}output_t;

/*************************************************************/

/* Structures to communicate information with other codes using PLEPP */

typedef struct _comm_1_t{

  char         friend_name[64];
  int          nphy;
  char      ** phys;          // array of Physical Entities names
  int        * ids ;
  int          remote_rank;;

  /* recv */
  double     * xs;

  /* send */
  double     * pow;

  MPI_Comm   * intercomm;

}comm_1_t;

typedef struct _comm_t{

  int        kind;

  comm_1_t   comm_1;

}comm_t;

/*************************************************************/

typedef struct _tcontrol_t{

  double tf;
  double dt;

}tcontrol_t;

typedef struct _ctrlrod_t{

  char   name_ele[16]; /* elem physical entity name */
  char   name_nod[16]; /* node physical entity name */
  int    nfun;         /* function id number */
  f1d_t  *funins;      /* insertion value function of time */
  double n[3];         /* direction of control rod insertion */
  double p[3];         /* reference point from where insertion starts */
  double xsaval;       /* xsa value to perturb  */
  list_t elemv;        /* list of vol elem that the ctrlrod intersects*/
  list_t xsa;          /* list of absortion xs on element list*/

}ctrlrod_t;

typedef struct{

  double t0;
  double t;
  list_t time;

    int timedep;
    int kmode;  
    int mode;
    int exec;
    
}calcu_t;

/*************************************************************/

typedef struct _coupling_t{

  char     world[64];
  int      num_friends;
  char  ** friends;
    
}coupling_t;


#endif

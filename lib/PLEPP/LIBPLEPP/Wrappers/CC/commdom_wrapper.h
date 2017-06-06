//#include "commdom.hpp"

//=======================================================================||===//
//=======================================================================||===//
#ifndef COMMDOMM_WRAPPER_H
#define COMMDOMM_WRAPPER_H

//=======================================================================||===//
#ifdef __cplusplus
extern "C" {  // fortran visibility!!
#endif

  // loaders arrays
  void 
  commdom_create(); 

  void 
  commdom_delete(); 

  void
  commdom_get_argvs(char* ftype);


  void
  commdom_analyse_argvs(const char* ftype, int* ntype); 


  void 
  commdom_set_names(const char*, int*, const char*, int*);
  
  void 
  commdom_create_commij(int*, int*);

  void
  commdom_get_commij_size(int*);

  void
  commdom_get_commij(char*, int*, int*); 
/*
  void 
  commdom_sendrecv_int(int*, int*, int*, int*, int*, int*);

  void 
  commdom_sendrecv_real_(double*, int*, double*, int*, int*, int*);

  void
  commdom_who_iam(const char* fname, INTEGER* nname, INTEGER* ok); 

  void
  commdom_who_areyou(const char* fname, int* nname, int* ok); 
*/
#ifdef __cplusplus
}
#endif
//=======================================================================||===//


#endif //LOADER_ALYA_WRAPPER_H
//=======================================================================||===//
//=======================================================================||===//

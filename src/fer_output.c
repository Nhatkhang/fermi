/*
 * output.c - This files contains the output functions
 * 
 * print_struct - print internal structures of the code, useful tool
 *                for debugging
 *
 * save_vtk - Generates a vtk format file where are storage the mesh,
 *           the material properties distribution and the flux
 * 
 * Autor : Guido Giuntoli
 * Last Modification: 29/12/2015
 * 
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fermi.h"
#include "global.h"
#include "list.h"
#include "gmsh.h"

#define   VTK_POINT         1
#define   VTK_LINE          3
#define   VTK_TRIANGLE      5
#define   VTK_QUADRANGLE    9
#define   VTK_TETRAHEDRON   10
#define   VTK_HEXAHEDRON    12
#define   VTK_6N_PRISM      13

int print_struct(int step)
{

  char file_name[64];
  int  i,d,e,*pInt;
  gmshN_t * p_gmsh_node;
  gmshE_t * p_gmsh_elem;
  gmshP_t * p_gmsh_phys;
  node_list_t * onode,*pNod;
  ps_t    * ps;

  sprintf(file_name,"struct_r%d_s%d.dat",rank,step);

  FILE *fout = fopen(file_name,"w");

  if(!fout)
  {
    PetscPrintf(FERMI_Comm,"output.c: (output_print_structures) Error opening output file.");
    return 1;
  }

  fprintf(fout,"Structures - rank %d\n\n",rank);  

  fprintf(fout,"========================================\n");
  fprintf(fout,"list_nodes : \n");
  fprintf(fout,"sizelist : %d\n",list_nodes.sizelist);
  onode = list_nodes.head;
  for(i=0;i<list_nodes.sizelist;i++)
  {
    p_gmsh_node = (gmshN_t *)onode->data;
    fprintf(fout,"list_nodes - node : %d\n",i+1);
    fprintf(fout,"n : %d  coor :",p_gmsh_node->n);
    for(d=0;d<DIM;d++)
      fprintf(fout,"%4.3lf ",p_gmsh_node->coor[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"list_ghost : \n");
  fprintf(fout,"sizelist : %d\n",list_ghost.sizelist);
  onode = list_ghost.head;
  for(i=0;i<list_ghost.sizelist;i++)
  {
    p_gmsh_node = (gmshN_t *)onode->data;
    fprintf(fout,"list_ghost - node : %d\n",i+1);
    fprintf(fout,"n : %d  coor :",p_gmsh_node->n);
    for(d=0;d<DIM;d++)
      fprintf(fout,"%4.3lf ",p_gmsh_node->coor[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"list_elemv : \n");
  fprintf(fout,"sizelist : %d\n",list_elemv.sizelist);
  onode = list_elemv.head;
  for(i=0;i<list_elemv.sizelist;i++)
  {
    p_gmsh_elem = (gmshE_t *)onode->data;
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"npe : %d ",p_gmsh_elem->npe);
    fprintf(fout,"gmshid : %d ",p_gmsh_elem->gmshid);
    fprintf(fout,"node : ");
    for(d=0;d<p_gmsh_elem->npe;d++)
      fprintf(fout,"%3d ",p_gmsh_elem->node[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"list_elems : \n");
  fprintf(fout,"sizelist : %d\n",list_elems.sizelist);
  onode = list_elems.head;
  for(i=0;i<list_elems.sizelist;i++){
    p_gmsh_elem = (gmshE_t *)onode->data;
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"npe : %d ",p_gmsh_elem->npe);
    fprintf(fout,"gmshid : %d ",p_gmsh_elem->gmshid);
    fprintf(fout,"node : ");
    for(d=0;d<p_gmsh_elem->npe;d++)
      fprintf(fout,"%3d ",p_gmsh_elem->node[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"list_physe : \n");
  fprintf(fout,"sizelist : %d\n\n",list_physe.sizelist);
  onode = list_physe.head;
  for(i=0;i<list_physe.sizelist;i++)
  {
    p_gmsh_phys = (gmshP_t *)onode->data;
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"name : %s\n",p_gmsh_phys->name);
    fprintf(fout,"dim : %d\n",p_gmsh_phys->dim);
    fprintf(fout,"gmshid : %d\n\n",p_gmsh_phys->gmshid);
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"list_bound :\n");
  fprintf(fout,"sizelist : %d\n\n",list_bound.sizelist);
  onode = list_bound.head;
  for(i=0;i<list_bound.sizelist;i++)
  {
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"name : %s\n",((bound_t *)onode->data)->name); 
    fprintf(fout,"order: %d\n",((bound_t *)onode->data)->order); 
    fprintf(fout,"kind : %d\n",((bound_t *)onode->data)->kind); 
    fprintf(fout,"dimS : %d\n",((bound_t *)onode->data)->dimS); 
    fprintf(fout,"nnodes: % d ",((bound_t *)onode->data)->nodeL.sizelist);
    fprintf(fout,"nodes: ");
    pNod=((bound_t *)onode->data)->nodeL.head;
    while(pNod)
    {
      fprintf(fout,"%4d ",*(int*)pNod->data);
      pNod=pNod->next;
    }
    fprintf(fout,"\n\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");
  
  fprintf(fout,"========================================\n");
  fprintf(fout,"list_ctrlr :\n");
  fprintf(fout,"sizelist : %d\n\n",list_ctrlr.sizelist);
  onode = list_ctrlr.head;
  for(i=0;i<list_ctrlr.sizelist;i++)
  {
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"name_ele : %s\n",((ctrlrod_t *)onode->data)->name_ele); 
    fprintf(fout,"name_nod : %s\n",((ctrlrod_t *)onode->data)->name_nod); 
    fprintf(fout,"nfun : %d\n",((ctrlrod_t *)onode->data)->nfun); 
    fprintf(fout,"n : ");
    for(d=0;d<3;d++)
      fprintf(fout,"%lf ",((ctrlrod_t *)onode->data)->n[d]); 
    fprintf(fout,"\n");
    fprintf(fout,"p : "); 
    for(d=0;d<3;d++)
      fprintf(fout,"%lf ",((ctrlrod_t *)onode->data)->p[d]); 
    fprintf(fout,"\n");
    fprintf(fout,"xsaval : %lf\n",((ctrlrod_t *)onode->data)->xsaval); 
    fprintf(fout,"elemv : ");
    pNod=((ctrlrod_t *)onode->data)->elemv.head;
    while(pNod)
    {
      fprintf(fout,"%4d ",*(int*)pNod->data);
      pNod=pNod->next;
    }
    fprintf(fout,"\n");
    fprintf(fout,"xsa : ");
    pNod=((ctrlrod_t *)onode->data)->xsa.head;
    while(pNod)
    {
      fprintf(fout,"%lf ",*(double*)pNod->data);
      pNod=pNod->next;
    }
    fprintf(fout,"\n\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"Dirichlet indeces %3d : \n",ndir);
  for(i=0;i<ndir;i++)
    fprintf(fout,"%2d ",dirIndex[i]);
  fprintf(fout,"\n\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"npp : \n");
  for(i=0;i<nproc;i++)
    fprintf(fout,"rank %d : npp[%d] = %d \n",i,i,npp[i]);
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"mesh:\nnode:\n");
  fprintf(fout," local: \n");
  for(i=0;i<mesh.nnodes;i++){
    fprintf(fout,"%d ",i);
    for(d=0;d<3;d++)
      fprintf(fout,"%+4.3lf ",mesh.node[i].coor[d]);
    onode=mesh.node[i].elemvL.head;        
    fprintf(fout,"\nelemvL: "); 
    e=0;        
    while(onode)
    {
      pInt=(int*)onode->data;
      fprintf(fout,"%4d",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8)
    {
      fprintf(fout,"    "); 
      e++; 
    }
    onode=mesh.node[i].elemsL.head;        
    fprintf(fout,"\nelemsL: "); 
    e=0;        
    while(onode)
    {
      pInt=(int*)onode->data;
      fprintf(fout,"%4d ",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8)
    {
      fprintf(fout,"    "); 
      e++; 
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"ghosts:\n");
  for(i=0;i<mesh.nghost;i++){
    fprintf(fout,"%d ",i);
    for(d=0;d<3;d++)
      fprintf(fout,"%+4.3lf ",mesh.node[mesh.nnodes+i].coor[d]);
    onode=mesh.node[i].elemvL.head;        
    fprintf(fout,"\nelemvL: "); 
    e=0;        
    while(onode)
    {
      pInt=(int*)onode->data;
      fprintf(fout,"%4d",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8)
    {
      fprintf(fout,"    "); 
      e++; 
    }
    onode=mesh.node[mesh.nnodes+i].elemsL.head;        
    fprintf(fout,"\nelemsL: "); 
    e=0;        
    while(onode)
    {
      pInt=(int*)onode->data;
      fprintf(fout,"%4d ",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8)
    {
      fprintf(fout,"    "); 
      e++; 
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"elemv : \n");
  for(i=0;i<mesh.nelemv;i++){
    fprintf(fout,"e=%4d npe: %d\nnodel:",i,mesh.elemv[i].npe);
    for(d=0;d<mesh.elemv[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elemv[i].nodel[d]);
    fprintf(fout,"\nnodeg:");
    for(d=0;d<mesh.elemv[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elemv[i].nodeg[d]);
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"elems : \n");
  for(i=0;i<mesh.nelems;i++){
    fprintf(fout,"e=%4d npe: %d\nnodel:",i,mesh.elems[i].npe);
    for(d=0;d<mesh.elems[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elems[i].nodel[d]);
    fprintf(fout,"\nnodeg: ");
    for(d=0;d<mesh.elems[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elems[i].nodeg[d]);
    ps=(ps_t*)mesh.elems[i].prop;
    fprintf(fout,"gmshid: %4d \n",ps->gmshid);
  }
  fprintf(fout,"\n");

  fprintf(fout,"========================================\n");
  fprintf(fout,"loc2gold : \n");
  for(i=0;i<(mesh.nnodes+mesh.nghost);i++)
    fprintf(fout,"%5d %5d \n",i,loc2gold[i]);
  fprintf(fout,"\n");
  fprintf(fout,"loc2gnew : \n");
  for(i=0;i<(mesh.nnodes+mesh.nghost);i++)
    fprintf(fout,"%5d %5d \n",i,loc2gnew[i]);
  fprintf(fout,"\n");

  fclose(fout);     

  return 0;   
}

int print_out(Vec *phi, int step)
{

  /*
     It travels all the list.outpu and determines if it has to print something
     in this time step "step" the kind of the structure output_t means what it 
     has to print 

   */

  int           i;
  int           error;
  double        pow_tot = 0.0;
  Vec           x_local;
  node_list_t * pNod;
  output_t    * po;
        
  VecGhostUpdateBegin(*phi,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*phi,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*phi,&x_local);

  error = 0;

  pNod = list_outpu.head;

  while(pNod){

    po=(output_t *)pNod->data;

    switch(po->kind){

      case 1:

//        double        vloc[2], vglo[2];
//
//        bthdisf(po->phys, po->norm, vloc, u);
//        MPI_Reduce(vloc, vglo, 2, MPI_DOUBLE, MPI_SUM, 0, FERMI_Comm);
//
//        if(!rank){
//          fl=fopen(po->file,"a");
//        }
//
//        fprintf(fl,"%lf %lf \n",vglo[0],vglo[1]);
//
//        fclose(fl);
        break;

      case 2:

	/*
	   prints power on different physical entities on ASCII file

	   kind = 2  
	   file <file_name.dat>
	   nphy <# of physical entities>
	   <"phys_1"> <"phys_2"> ... <"phys_n">

	 */
        

	fer_pow_phys( po->kind_2.nphy, po->kind_2.ids, po->kind_2.pow );

	if(local_rank == 0){
	  for(i=0;i<po->kind_2.nphy;i++){
	    fprintf(po->kind_2.fp, "%e ",po->kind_2.pow[i]);
	    pow_tot += po->kind_2.pow[i];
	  }
	  // we print the sum of all the powers in the Physical Entities
	  fprintf(po->kind_2.fp, "%e\n",pow_tot);
	}

        break;

      default:
        return 1;

    }

    if(error){
      return 1;
    }

    pNod=pNod->next;
  }

  return 0;
}

/****************************************************************************************************/

int print_vtk(const char *name) {

  int n,e,d,count,index;
  double vald;
  char  filevtk[32],ending[16];
  FILE *vtkf;
  node_list_t *pn,*pe,*px;

  strcpy(filevtk,name);
  filevtk[strlen(name)]='\0';
  sprintf(ending,".vtk");
  strcat(filevtk,ending);
  vtkf = fopen(filevtk, "w");

  /************************************************************/
  /*Mesh geometry data                                        */
  /************************************************************/    
  fprintf(vtkf, "# vtk DataFile Version 2.0\n");
  fprintf(vtkf, "Fermi\n");
  fprintf(vtkf, "ASCII\n");
  fprintf(vtkf, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(vtkf, "POINTS %d double\n", (mesh.nnodes+mesh.nghost));

  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for(d=0;d<3;d++)
      fprintf(vtkf, "%lf ", mesh.node[n].coor[d]);
    fprintf(vtkf, "\n");
  }

  count=0;
  for(e=0;e<mesh.nelemv;e++)
    count += mesh.elemv[e].npe + 1;
  fprintf(vtkf, "CELLS %d %d\n", mesh.nelemv, count);
  for (e=0;e<mesh.nelemv;e++){
    fprintf(vtkf, "%d ", mesh.elemv[e].npe);
    for (n=0;n<mesh.elemv[e].npe;n++)
      fprintf(vtkf, "%d ", mesh.elemv[e].nodel[n]);
    fprintf(vtkf, "\n");
  }

  fprintf(vtkf, "CELL_TYPES %i\n", mesh.nelemv);
  for (e=0;e<mesh.nelemv;e++)
    fprintf(vtkf, "%d\n",vtkcode(DIM,mesh.elemv[e].npe));  


  /************************************************************/
  /* DATA IN NODES                                            */
  /************************************************************/
  fprintf(vtkf, "POINT_DATA %i\n",(mesh.nnodes+mesh.nghost));

  /* Ownership */
  fprintf(vtkf, "SCALARS ownership FLOAT\n");
  fprintf(vtkf, "LOOKUP_TABLE default\n");
  for(n=0;n<mesh.nnodes;n++)
    fprintf(vtkf,"%lf\n",1.0);
  for(n=0;n<mesh.nghost;n++)
    fprintf(vtkf,"%lf\n",0.0);

  /* Displacement */
  VecGhostUpdateBegin(phi_n,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(phi_n,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(phi_n,&xlocal);

  for(e=0;e<egn;e++)
  {
    fprintf(vtkf,"SCALARS flux_%d FLOAT\n",e);
    fprintf(vtkf, "LOOKUP_TABLE default\n");
    for(n=0;n<(mesh.nnodes+mesh.nghost);n++)
    {
      index=n*egn+e;
      VecGetValues(xlocal,1,&index,&vald);
      fprintf(vtkf,"%lf\n",vald);
    }        
  }

  /************************************************************/
  /* DATA IN CELLS                                            */
  /************************************************************/      
  fprintf(vtkf, "CELL_DATA %i\n",mesh.nelemv);

  pn=list_ctrlr.head;  
  while(pn)
  {
    fprintf(vtkf, "SCALARS rod_%s FLOAT\n",((ctrlrod_t*)pn->data)->name_ele);
    fprintf(vtkf, "LOOKUP_TABLE default\n");
    pe=((ctrlrod_t*)pn->data)->elemv.head;
    px=((ctrlrod_t*)pn->data)->xsa.head;
    for (e=0;e<mesh.nelemv;e++)
    {
      if(pe)
      {
        if(e==*(int*)pe->data)
        {
          fprintf(vtkf, "%lf\n",*(double*)px->data);
          pe=pe->next;
          px=px->next;
        }else{
          fprintf(vtkf, "%lf\n",0.0);
        }
      }else
      {
        fprintf(vtkf, "%lf\n",0.0);
      }
    }
    pn=pn->next;
  }

  for (e=0;e<egn;e++)
  {
    fprintf(vtkf, "SCALARS xs_a%d FLOAT\n",e);
    fprintf(vtkf, "LOOKUP_TABLE default\n");
    for (d=0;d<mesh.nelemv;d++)
      fprintf(vtkf, "%lf\n", ((pv_t*)mesh.elemv[d].prop)->xs_a[e]);

    fprintf(vtkf, "SCALARS nxs_f%d FLOAT\n",e);
    fprintf(vtkf, "LOOKUP_TABLE default\n");
    for (d=0;d<mesh.nelemv;d++)
      fprintf(vtkf, "%lf\n", ((pv_t*)mesh.elemv[d].prop)->nxs_f[e]);

    fprintf(vtkf, "SCALARS D%d FLOAT\n",e);
    fprintf(vtkf, "LOOKUP_TABLE default\n");
    for (d=0;d<mesh.nelemv;d++)
      fprintf(vtkf, "%lf\n", ((pv_t*)mesh.elemv[d].prop)->D[e]);
  } 
  fclose(vtkf);
  return 0;

}   

int printMatrixR(char *name,double *A, int m, int n){

  int i,j;

  FILE *f=fopen(name,"w");

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      fprintf(f,"%lf ",A[i*n+j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);   
  return 0;
}

int printMatrixRC(char *name,double **A, int m, int n){

  int i,j;

  FILE *f=fopen(name,"w");

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      fprintf(f,"%lf ",A[i][j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);   
  return 0;
}

int printMat(char *name,Mat *A,int n){

  int i,j;
  double val;

  FILE *f=fopen(name,"w");

  MatGetOwnershipRange(*A,&Istart,&Iend);
  for(i=Istart;i<Iend;i++){
    for(j=0;j<n;j++){
      MatGetValues(*A,1,&i,1,&j,&val);
      fprintf(f,"%lf ",val);
    }
    fprintf(f,"\n");
  }
  fclose(f);   
  return 0;
}

int printVector(char *name,double *vec, int m){

  int j;

  FILE *f=fopen(name,"w");

  for(j=0;j<m;j++){
    fprintf(f,"%lf %lf\n ",vec[j],vec[j]);
  }

  fclose(f);   
  return 0;
}

int printVec(char *name,Vec *vec, int m){

  /* Prints the values of a Vec, it can include ghost values*/

  int j;
  double val;

  FILE *f=fopen(name,"w");

  VecGhostGetLocalForm(*vec,&xlocal);

  for(j=0;j<m;j++){
    VecGetValues(xlocal,1,&j,&val);
    fprintf(f,"%lf %lf\n ",val,val);
  }

  fclose(f);   
  return 0;
}

int vtkcode(int dim,int npe){

  switch(dim){
    case 1:
      switch(npe){
        case 2 :
          return VTK_LINE;
        default:
          return -1;
      }
    case 2:
      switch(npe){
        case 3 :
          return VTK_TRIANGLE;
        case 4 :
          return VTK_QUADRANGLE;
        default:
          return -1;
      }
    case 3:
      switch(npe){
        case 4 :
          return VTK_TETRAHEDRON;
        case 6 :
          return VTK_6N_PRISM;  
        case 8 :
          return VTK_HEXAHEDRON;  
        default:
          return -1;
      }
    default:
      return -1;
  }
}

/* Rutines for calculating control rods movements*/

#include "fermi.h"

int ferirods(void)
{
  /* determines the elems & elemv list for each control rod in the list 
   *
   * elemv list: are those volumetric elements that are contain in the 
   * volume with physical entity name_ele
   */
  node_list_t *pr,*pp;
  int e,d,error,gmshid,flag=0;
  double val;

  pr=list_ctrlr.head;
  while(pr)
  {
    /* read elements */
    pp=list_physe.head;
    while(pp){
      if(!strcmp(((gmshP_t*)pp->data)->name,((ctrlrod_t*)pr->data)->name_ele))
        break;
      pp=pp->next;
    }
    if(!pp)
    {
      PetscPrintf(FERMI_Comm,"ferrods.c:the rod %s doesn't have physe.\n",((ctrlrod_t*)pr->data)->name_ele); 
      return 1;
    }
    gmshid=((gmshP_t*)pp->data)->gmshid;
    list_init(&((ctrlrod_t*)pr->data)->elemv,sizeof(int),NULL);
    list_init(&((ctrlrod_t*)pr->data)->xsa,sizeof(double),NULL);
    for(e=0;e<mesh.nelemv;e++)
    {
      if(((pv_t*)mesh.elemv[e].prop)->gmshid==gmshid)
      {
        val=0.0;
        list_insertlast(&((ctrlrod_t*)pr->data)->elemv,(void*)&e);
        list_insertlast(&((ctrlrod_t*)pr->data)->xsa,(void*)&val);
      }
    }
    /* read reference node */
    pp=list_physe.head;
    while(pp){
      if(!strcmp(((gmshP_t*)pp->data)->name,((ctrlrod_t*)pr->data)->name_nod))
        break;
      pp=pp->next;
    }
    if(!pp)
    {
      PetscPrintf(FERMI_Comm,"ferrods.c:the rod %s doesn't have physe.\n",((ctrlrod_t*)pr->data)->name_nod); 
      return 1;
    }
    gmshid=((gmshP_t*)pp->data)->gmshid;
    for(e=0;e<mesh.nelems;e++)
    {
      if(((ps_t*)mesh.elems[e].prop)->gmshid==gmshid)
      {
        for(d=0;d<3;d++)
          ((ctrlrod_t*)pr->data)->p[d]=mesh.node[mesh.elems[e].nodel[0]].coor[d];
        for(d=0;d<nproc;d++)
          error=MPI_Send((void*)((ctrlrod_t*)pr->data)->p,3,MPI_DOUBLE,d,0,FERMI_Comm);
        flag=1;
        break;
      }
    }

    if(!flag)
    {
      error=MPI_Recv((void*)((ctrlrod_t*)pr->data)->p,3,MPI_DOUBLE,1,0,FERMI_Comm,MPI_STATUS_IGNORE);
      if(error)
        return 1;
    }

    pr=pr->next;
    MPI_Barrier(FERMI_Comm);
  }
  return 0;
}

int fersrods(double t){
  /*fils the list xsa of the ctrlrod_t structure according to the value of time*/
  double h,bari[3],dprod;
  int n,d,npe;
  node_list_t *pr,*pe,*px;
  pr=list_ctrlr.head;
  while(pr)
  {
    f1d_eval(t,((ctrlrod_t*)pr->data)->funins,&h);
    pe=((ctrlrod_t*)pr->data)->elemv.head;
    px=((ctrlrod_t*)pr->data)->xsa.head;
    while(pe)
    {
      memset(bari,0.0,3*sizeof(double));
      npe=mesh.elemv[*(int*)pe->data].npe;
      for(n=0;n<npe;n++)
      {
        for(d=0;d<3;d++)
          bari[d]+=mesh.node[mesh.elemv[*(int*)pe->data].nodel[n]].coor[d];
      }
      for(d=0;d<3;d++)
        bari[d]/=npe;
      for(d=0;d<3;d++)
        bari[d]-=((ctrlrod_t*)pr->data)->p[d];
      dprod=0.0;
      for(d=0;d<3;d++)
      {
        dprod+=bari[d]*((ctrlrod_t*)pr->data)->n[d];
      }
      if(dprod<h)
      {
        *(double*)px->data=((ctrlrod_t*)pr->data)->xsaval;
      }else{
        *(double*)px->data=0.0;
      }
      pe=pe->next;
      px=px->next;
    }
    pr=pr->next;
  }
  return 0;
}

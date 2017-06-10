/* 
 * parser.c : parser rutine for FERMI inputs
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "fermi.h"
#include "list.h"
#include "types.h"
#include "matlaw.h"
#include "mesh.h"

#define NBUF 128

int parse_input(void){

  if(access(inputfile,F_OK) == -1){
    PetscPrintf(FERMI_Comm,"parser.c: file %s not found.\n",inputfile);
    return 1;
  }

  if(parse_mesh())return 1; 
  if(parse_mats())return 2;
  if(parse_mode())return 3;
  if(parse_func())return 4; 
  if(parse_boun())return 5; 
  if(parse_crod())return 6; 
  if(parse_outp())return 7; 
  if(parse_communication())return 8; 
  return 0;
}

int parse_mesh(void)
{

  FILE *file= fopen(inputfile,"r");
  char *data,buf[NBUF];
  int  fl=0,com=0,ln=0;

  while(fgets(buf,NBUF,file))
  {
    ln++;
    data=strtok(buf," \n");
    if(data)
    {
      if(!strcmp(data,"$Mesh"))
      {
        fl=1;
      }else if(!strcmp(data,"$EndMesh"))
      {
        if(com==7)
        { 
          return 0;
        }else if(nproc==1 && com>=1){
          return 0;
        }else if(nproc>1&&(!(com&4)||!(com&8)) )
        {
          PetscPrintf(FERMI_Comm,"parser.c:part file NF.\n");
          return 1;
        }else{
          PetscPrintf(FERMI_Comm,"parser.c:$Mesh sect BF.\n");
          return 1;    
        }
      }
      if(fl==2){
        if(strcmp(data,"mesh_file") == 0){    
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:meshfile line %d.\n",ln); 
            return 1;
          }
          strcpy(meshfile,data);   
          if(access(meshfile,F_OK) == -1){
            PetscPrintf(FERMI_Comm,"parser.c:%s NF.\n",meshfile); 
            return 1;
          }
          com=com|1;  
        }else if(strcmp(data,"parfile_e") == 0){    
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:partfile line %d.\n",ln);
            return 1;
          }
          strcpy(epartfile,data);
          if(access(epartfile,F_OK) == -1){
            PetscPrintf(FERMI_Comm,"parser.c:%s NF.\n",epartfile); 
            return 1;
          }
          com=com|2;
        }else if(strcmp(data,"parfile_n") == 0){    
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:partfile line %d.\n",ln);
            return 1;
          }
          strcpy(npartfile,data);
          if(access(npartfile,F_OK) == -1){
            PetscPrintf(FERMI_Comm,"parser.c:%s NF.\n",npartfile); 
            return 1;
          }
          com=com|4;     
        }else if(strcmp(data,"#")!=0){    
          PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
          return 1;
        }    
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 1;
}

int parse_mats(void)
{
  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF],bufcpy[NBUF];
  int  i,fl=0,ln=0,com=0;
  pvl_t mat;

  while(fgets(buf,NBUF,file))
  {
    ln++;
    strcpy(bufcpy,buf);
    data=strtok(buf," \n");
    if(data)
    {
      if(!strcmp(data,"$Xs"))
      {
        fl=1;
      }else if(!strcmp(data,"$EndXs"))
      {
        if(!list_mater.sizelist)
        {
          PetscPrintf(FERMI_Comm,"parser.c:no materials specified.\n");
          return 1;     
        }
        return 0;
      }
      if(fl==2)
      {
        if(!strcmp(data,"egn"))
        {
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
            return 1;
          }    
          egn=atoi(data);
          if(egn<1)
          {
            PetscPrintf(FERMI_Comm,"parser.c:egn should positive at line %d.\n",ln); 
            return 1;
          }
          veloc=(double*)calloc(egn,sizeof(double));
	  switch(egn){
	    case 1:
	      nxs_mat=5;
	      break;
	    case 2:
	      nxs_mat=12;
	      break;
	    default:
	      return 1;
	  }
          com=com|1;
        }else if(!strcmp(data,"pgn"))
        {
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
            return 1;
          }    
          pgn=atoi(data);
          if(pgn<1)
          {
            PetscPrintf(FERMI_Comm,"parser.c:pgn should positive at line %d.\n",ln); 
            return 1;
          }
          beta=(double*)calloc(pgn,sizeof(double));
          lambda=(double*)calloc(pgn,sizeof(double));
          chi=(double*)calloc(pgn,sizeof(double));
          com=com|2;
        }else if(!strcmp(data,"vel"))
        {
          for(i=0;i<egn;i++)
          {
            data=strtok(NULL," \n");    
            if(!data){
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
              return 1;
            }    
            veloc[i]=atof(data);
          }
          com=com|4;
        }else if(!strcmp(data,"kyn"))
        {
          for(i=0;i<pgn;i++)
          {
            data=strtok(NULL," \n");    
            if(!data){
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
              return 1;
            }    
            beta[i]=atof(data);
          }
          for(i=0;i<pgn;i++)
          {
            data=strtok(NULL," \n");    
            if(!data){
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
              return 1;
            }    
            lambda[i]=atof(data);
          }
          for(i=0;i<pgn;i++)
          {
            data=strtok(NULL," \n");    
            if(!data){
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
              return 1;
            }    
            chi[i]=atof(data);
          }
          com=com|8;
        }else if(data[0]!='#')
        {
          if((com&1)!=1)
          {
            PetscPrintf(FERMI_Comm,"parser.c:egn should be before xs values at line %d.\n",ln); 
            return 1;
          }    
          mat.D    =(double*)calloc(egn,sizeof(double));
          mat.xs_a =(double*)calloc(egn,sizeof(double));
          mat.xs_s =(double*)calloc(egn*(egn-1),sizeof(double));
          mat.nxs_f=(double*)calloc(egn,sizeof(double));
          mat.exs_f=(double*)calloc(egn,sizeof(double));
          mat.chi  =(double*)calloc(egn,sizeof(double));
          if(parse_material(bufcpy,&mat))
          {
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln); 
            return 1;
          }
          if(list_insert_se(&list_mater,(void*)&mat))
            return 1;
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 1;
}

int parse_mode(void)
{
  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF];
  int  fl=0,com=0,ln=0,auxi;
  double auxd;
  list_t list_ts, list_tf, list_dt;        

  list_init(&calcu.time,sizeof(tcontrol_t),cmp_time);
  list_init(&list_tf,sizeof(double),cmp_int);
  list_init(&list_dt,sizeof(double),cmp_int);
  list_init(&list_ts,sizeof(double),cmp_int);

  while(fgets(buf,NBUF,file))
  {
    ln++;
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Mode")){
        fl=1;
      }else if(!strcmp(data,"$EndMode")){
        fclose(file);
        if((com&64)!=64)
        {
          power=1.0e6;
        }else{
          if(power<=1.0e-10)
          {
            PetscPrintf(FERMI_Comm,"parser.c:p0 <= 0 line %d.\n",ln); 
            return 1;
          }
        }
        node_list_t *pna=list_tf.head,*pnb=list_dt.head;
        tcontrol_t time;
        while(pna && pnb){
          time.tf=*(double*)pna->data;
          time.dt=*(double*)pnb->data;
          if(list_insert_se(&calcu.time,(void*)&time))
            return 1;
          pna=pna->next;
          pnb=pnb->next;
        }
        if(pna!=NULL||pnb!=NULL){
          PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
          return 1;
        }else if( (com&24)==24 ){
          PetscPrintf(FERMI_Comm,"parser.c:dt and ts specified at the same time.\n");
        }
        return 0;    
      }
      if(fl==2){
        if(strcmp(data,"timedep") == 0){
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
          }
          if(!strcmp(data,"QSTATIC")){    
            calcu.timedep=QS;
          }else if(!strcmp(data,"DYNAMIC")){    
            calcu.timedep=TR;
          }else{    
            PetscPrintf(FERMI_Comm,"parser.c:Invalid option at line %d.\n",ln);
            return 1;
          }
          com=com|2;
        }else if(strcmp(data,"t0") == 0)
        {
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          calcu.t0=atof(data);
        }else if(!strcmp(data,"tf"))
        {
          data = strtok(NULL," \n");
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          while(data){
            auxd=atof(data);
            list_insertlast(&list_tf,(void*)&auxd);
            data=strtok(NULL," \n");
          }
          com=com|4;
        }else if(strcmp(data,"dt") == 0){
          data = strtok(NULL," \n");
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          while(data){
            auxd=atof(data);
            list_insertlast(&list_dt,(void*)&auxd);
            data=strtok(NULL," \n");
          }
          com=com|8;
        }else if(strcmp(data,"tsteps") == 0){
          data = strtok(NULL," \n");
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          while(data){
            auxi=atoi(data);
            list_insertlast(&list_ts,(void*)&auxi);
            data=strtok(NULL," \n");
          }
          com=com|16;
        }else if(!strcmp(data,"kmode"))
        {
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          if(!strcmp(data,"K1")){    
            calcu.kmode=K1;
          }else{    
            PetscPrintf(FERMI_Comm,"parser.c:Invalid option at line %d.\n",ln);
            return 1;
          }
          com=com|32;
        }else if(!strcmp(data,"p0"))
        {
          data = strtok(NULL," \n");    
          if(!data)
          {
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          power=atof(data);
          com=com|64;
        }else if(strcmp(data,"#") != 0){    
          PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
          return 1;
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  return 1;    
}

int parse_crod(void){

  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF],bufcpy[NBUF];
  int fl=0,ln=0,i;
  ctrlrod_t ctrlrod;
  node_list_t *nod;

  while(fgets(buf,NBUF,file)){
    ln++;
    strcpy(bufcpy,buf);
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Ctrlrod")){
        fl=1;
      }else if(!strcmp(data,"$EndCtrlrod")){
        list_insertlast(&list_ctrlr,(void*)&ctrlrod);
        fl=0;
      }
      if(fl==2){
        if(data[0]!='#'){
          if(!strcmp(data,"name_ele")){
            data=strtok(NULL," \n");
            if(!data)
            {
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
              return 1;
            }
            strcpy(ctrlrod.name_ele,data);
          }else if(!strcmp(data,"name_nod")){
            data=strtok(NULL," \n");
            if(!data)
            {
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
              return 1;
            }
            strcpy(ctrlrod.name_nod,data);
          }else if(!strcmp(data,"func")){
            data=strtok(NULL," \n");
            if(!data)
            {
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
              return 1;
            }
            ctrlrod.nfun=atoi(data);
            nod=list_fun1d.head;
            while(nod)
            {
              if(ctrlrod.nfun==((f1d_t*)nod->data)->fnum)
                break;
              nod=nod->next;
            }
            if(!nod)
            {
              PetscPrintf(FERMI_Comm,"parser.c:func %d NF line %d.\n",ctrlrod.nfun,ln);
              return 1;
            }
            ctrlrod.funins=nod->data;
          }else if(!strcmp(data,"norm")){
            for(i=0;i<3;i++)
            {
              data=strtok(NULL," \n");
              if(!data)
              {
                PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
                return 1;
              }
              ctrlrod.n[i]=atoi(data);
            }
          }else if(!strcmp(data,"xsa")){
            data=strtok(NULL," \n");
            if(!data)
            {
              PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
              return 1;
            }
            ctrlrod.xsaval=atof(data);
          }
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 0;
}

int parse_func(void)
{
  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF];
  int fl=0,com=0,ln=0,i,fs;
  double *xy;
  f1d_t * f1d;
  list_t list_xy;

  list_init(&list_xy, 2*sizeof(double), cmp_dou);
  while(fgets(buf,NBUF,file))
  {
    ln++;
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Function"))
      {
        fl=1;
        fs=0;
      }else if(!strcmp(data,"$EndFunction"))
      {
        if(com==15){ 
          com=0;
          fl=0;
        }else if(fl==0){
          PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
          return 1;
        }else{
          PetscPrintf(FERMI_Comm,"parser.c:incomplete $Function section.\n");
          return 1;
        }
      }
      if(fl==2){
        if(strcmp(data,"funcknd") == 0){
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          if(!strcmp(data,"1D")){
            f1d=(f1d_t*)calloc(1,sizeof(f1d_t));
            com=com|1;
          }
        }else if(strcmp(data,"funcnum") == 0){
          if(!f1d)
            return 1;
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          f1d->fnum=atoi(data);
          com=com|2;
        }else if(strcmp(data,"funcint") == 0){
          if(!f1d)
            return 1;
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          if(!strcmp(data,"INTER1")){
            f1d->inter=atoi(data);
          }    
          com=com|4;
        }else if(strcmp(data,"start") == 0){
          fs=1;
          if( (com & 1 ) != 1){
            PetscPrintf(FERMI_Comm,"parser.c:<funckind> should be < start>.\n");  
            return 1;
          }
        }else if(strcmp(data,"end") == 0){
          if(fs==0){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          if(!list_xy.sizelist){
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          f1d->n=list_xy.sizelist;
          f1d->x=(double*)calloc(list_xy.sizelist,sizeof(double));
          f1d->y=(double*)calloc(list_xy.sizelist,sizeof(double));
          for(i=0;i<f1d->n;i++){
            f1d->x[i]=*(((double*)(list_xy.head->data))+0);
            f1d->y[i]=*(((double*)(list_xy.head->data))+1);
            list_delfirst(&list_xy);
          }
          list_insert_se(&list_fun1d,(void*)f1d);
          com=com|8;
        }else if(fs){
          xy=(double*)calloc(2,sizeof(double));
          xy[0]=atof(data);
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          xy[1]=atof(data);
          list_insert_se(&list_xy,(void*)xy);
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 0;    
}       

int parse_boun(void){

  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF],bufcpy[NBUF];
  int fl=0,ln=0;
  bound_t bou;

  while(fgets(buf,NBUF,file)){
    ln++;
    strcpy(bufcpy,buf);
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Boundary")){
        fl=1;
      }else if(!strcmp(data,"$EndBoundary")){
        if(!list_bound.sizelist){ 
        PetscPrintf(FERMI_Comm,"parser.c: No boundaries in %s.\n",inputfile);
          return 1;
        }
        return 0;    
      }
      if(fl==2){
        if(data[0]!='#'){
          if(parse_boundary(bufcpy,&bou))
          {
            PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          list_insert_se(&list_bound,(void*)&bou);
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 1;
}

int parse_material(char *buff, pvl_t *mat){

  char *data;
  int i;

  if(!strlen(buff))
    return 1;
  data = strtok(buff," \n");
  if(!data)
    return 1;
  if(data[0]!='\"'||data[strlen(data)-1]!='\"')
    return 1;
  strcpy(mat->name,data);

  /* second column is 0|1 which means if the mat has or not has precursors */ 
  data = strtok(NULL," \n");
  if(!data)
    return 1;
  if(atoi(data)==1 || atoi(data)==0){
    mat->hasprec=atoi(data);
  }

  /* diffusion coeficients */ 
  for(i=0;i<egn;i++){
    data = strtok(NULL," \n");
    if(!data)
      return 1;
    mat->D[i]=atof(data);
  }

  /* absortion cross sections */ 
  for(i=0;i<egn;i++){
    data = strtok(NULL," \n");
    if(!data)
      return 1;
    mat->xs_a[i]=atof(data);
  } 


  /* scattering cross sections */ 
  for(i=0;i<egn*(egn-1);i++){
    data = strtok(NULL," \n");
    if(!data)
      return 1;
    mat->xs_s[i]=atof(data);
  } 

  /* fission cross sections */ 
  for(i=0;i<egn;i++){
    data = strtok(NULL," \n");
    if(!data)
      return 1;
    mat->nxs_f[i]=atof(data);
  } 

  /* energy cross sections */ 
  for(i=0;i<egn;i++){
    data = strtok(NULL," \n");
    if(!data)
      return 1;
    mat->exs_f[i]=atof(data);
  }

  /* fission spectrum */ 
  for(i=0;i<egn;i++){
    data = strtok(NULL," \n");
    if(!data)
      return 1;
    mat->chi[i]=atof(data);
  }
  return 0;
}

/*****************************************************************************************************/

int parse_outp(void)
{

  /*
     Parse for $Output word on parser
     it could be repeated, each one has a different kind
     which means what they want to output
   */

  FILE    * file= fopen(inputfile,"r");
  char    * data,buf[NBUF], buf_a[NBUF];
  int       fl, ln, i;
  output_t  output;

  ln = 0;
  fl = 0;
  while(fgets(buf,NBUF,file)){
    ln++;
    strcpy(buf_a, buf);
    data=strtok(buf_a," \n");

    if(data){ // if the line is not empty

      if(!strcmp(data,"$Output")){
	fl=1;
      }
      else if(!strcmp(data,"$EndOutput")){
	if(output.kind==1){ 

	  list_insertlast(&list_outpu,(void*)&output);

	}
	else if(output.kind==2){ 

	  list_insertlast(&list_outpu,(void*)&output);

	}
	else{
	  PetscPrintf(FERMI_Comm,"parser.c:BF line %d.\n",ln);
	  return 1;
	}
	fl = 0;
      }

      if(fl==2 && data[0]!='#'){

	// we are in the line after $Output
	if( get_int(buf,"kind",&output.kind)) 
	  return 1;  

	switch(output.kind){

	  case 1:
	    break;

	  case 2:

	    /*
	       kind = 2  

	       power on different physical entities on ASCII file
	       file <file.dat>
	       nphy <num_of_phys>
	       "phys_1" "phys_2" ... "phys_n"

	     */

	    if(!fgets(buf,NBUF,file)) 
	      return 1;
	    ln ++;

	    if( get_char(buf,"file",output.kind_2.file)) 
	      return 1;

	    if(!fgets(buf,NBUF,file)) 
	      return 1;
	    ln ++;

	    if( get_int(buf,"nphy",&output.kind_2.nphy)) 
	      return 1;  
	    output.kind_2.phys = malloc(output.kind_2.nphy*sizeof(char*));
	    output.kind_2.ids  = malloc(output.kind_2.nphy*sizeof(int*));
	    output.kind_2.pow  = malloc(output.kind_2.nphy*sizeof(double*));
	    for(i = 0; i < output.kind_2.nphy; i++){
	      output.kind_2.phys[i] = (char *)malloc(16*sizeof(char));
	    }

	    if(!fgets(buf,NBUF,file)) 
	      return 1;
	    ln ++;

	    // now we read the physical entities names
	    data = strtok(buf," \n");
	    i = 0;
	    while(i < output.kind_2.nphy && data != NULL){
	      strcpy( output.kind_2.phys[i],data);
	      data = strtok(NULL," \n");
	      i++;
	    }
	    if( i != output.kind_2.nphy ){
	      return 1;
	    }

	    break;

	  default:
	    break;
	}

      }
      if(fl==1)
	fl=2;
    }
  }
  fclose(file);
  return 0;
}

/*****************************************************************************************************/

int parse_communication(void)
{

  /*
     Parse for $Communication word on input file

     The user is able to specify what "kind" 
     of communication wants to perform 
     
     kind 1: means (recv) cross sections (send) powers

     Notes:

     -> $Communication could be repeated various time on input 
     file with different arguments (and maybe the same probably as 
     an error)

   */

  FILE      *file = fopen(inputfile,"r");
  char      *data,buf[NBUF], buf_a[NBUF];
  int       fl=0, ln=0, i;
  comm_t    comm;

  while(fgets(buf,NBUF,file))
  {

    ln++;
    strcpy(buf_a, buf);
    data=strtok(buf_a," \n");

    if(data){ // if the line is not empty

      if(!strcmp(data,"$Communication")){
	fl=1;
      }
      else if(!strcmp(data,"$EndCommunication")){
	list_insertlast(&list_comms,(void*)&comm);
	fl = 0;
      }
    }

    if(fl==2 && data[0]!='#'){

      // we are in the line after $Communication
      if( get_int(buf,"kind",&comm.kind)) 
	return 1;  

      switch(comm.kind){

	case 1:

	  /*
	     kind     1
	     friend   <friend_name>
	     nphy     <num_of_phys>
	     <phys_1> <phys_2>      ... <phys_n>

	     -> friend is the program to communicate information
	     -> xs    on different physical are going to be recv
	     -> power on different physical are going to be send

	   */

          // Search for its friend's name
	  if(!fgets(buf,NBUF,file)) 
	    return 1;
	  ln ++;
	  if( get_char(buf,"friend",comm.comm_1.friend_name)) 
	    return 1;

	  if(!fgets(buf,NBUF,file)) 
	    return 1;
	  ln ++;
	  if( get_int(buf,"nphy",&comm.comm_1.nphy)) 
	    return 1;  

	  comm.comm_1.phys = malloc(comm.comm_1.nphy*sizeof(char**));
	  comm.comm_1.ids  = malloc(comm.comm_1.nphy*sizeof(int*));
	  comm.comm_1.pow  = malloc(comm.comm_1.nphy*sizeof(double*));
	  comm.comm_1.xs   = malloc(comm.comm_1.nphy * \
	      nxs_mat *sizeof(double*));

	  if(!fgets(buf,NBUF,file)) 
	    return 1;
	  ln ++;

	  // now we read the physical entities names
	  data = strtok(buf," \n");
	  i = 0;
	  while(i < comm.comm_1.nphy && data != NULL){
	    comm.comm_1.phys[i] = strdup(data);
	    data = strtok(NULL," \n");
	    i++;
	  }
	  if( i != comm.comm_1.nphy ){
	    return 1;
	  }

	  break;

	default:
	  break;
      }

    }
    if(fl==1)
      fl=2;
  }
  fclose(file);
  return 0;
}

/*****************************************************************************************************/

int get_int(char *buf, const char *name,int *a)
{

  /*
     Looks in "buf" for :
     "name" <(int) a>
     returns 1 if this is not correct
             0 if this is 0K
   */

  char *data;

  data = strtok(buf," \n");    

  if(!data){
    PetscPrintf(FERMI_Comm,"%s expected.\n",name);
    return 1;
  }

  if(strcmp(data,name)){
    PetscPrintf(FERMI_Comm,"%s expected.\n",name);
    return 1;
  }

  data = strtok(NULL," \n");    
  if(!data){
    PetscPrintf(FERMI_Comm,"%s value expected.\n",name);
    return 1;
  }
  *a = atoi(data);

  return 0;
}

/*****************************************************************************************************/

int get_char(char *buf, const char *name,char *a)
{

  /*
     Looks in "buf" for :
     "name" <(char) a>
     returns 1 if this is not correct
             0 if this is 0K
   */

  char *data;

  data = strtok(buf," \n");    

  if(!data){
    PetscPrintf(FERMI_Comm,"%s expected.\n",name);
    return 1;
  }

  if(strcmp(data,name)){
    PetscPrintf(FERMI_Comm,"%s expected.\n",name);
    return 1;
  }

  data = strtok(NULL," \n");    
  if(!data){
    PetscPrintf(FERMI_Comm,"%s value expected.\n",name);
    return 1;
  }
  strcpy(a,data);

  return 0;
}

/*****************************************************************************************************/

int cmp_mat(void *a, void *b){

  return (strcmp(((pvl_t*)a)->name,((pvl_t*)b)->name)==0)?0:-1;
}

/*****************************************************************************************************/

int parse_coupling(const char file_c[])
{

  FILE * file = fopen(file_c,"r");
  char * data, buf[NBUF];
  int    ln=0, fl=0, i;

  while(fgets(buf,NBUF,file)){

    ln++;
    data=strtok(buf," \n");

    if(data){
      if(data[0]!='#'){

	switch(fl){

	  case 0:

	    // read world name
	    strcpy(coupling.world, data);
	    fl ++;
	    break;

	  case 1:

	    // read number of friend number
	    coupling.num_friends = atoi(data);
	    fl ++;
	    break;

	  case 2:

	    // read friend names
	    coupling.friends = malloc(coupling.num_friends * sizeof(char*));
	    for(i=0;i<coupling.num_friends;i++){

	      coupling.friends[i] = malloc(64 * sizeof(char));

	      if(data == NULL){
		PetscPrintf(FERMI_Comm,"you should give at least %d friend names at line %d in %s.\n",
		    coupling.num_friends, ln, file_c);
		return 1;
	      }
	      strcpy( coupling.friends[i], data);
	      data=strtok(NULL," \n");

	    }
	    fl ++;
	    break;

	  case 3:
	    return 0;

	  default:
	    return 1;

	}
      }
    }
  }

  fclose(file);
  return 0;
}

/*****************************************************************************************************/

int parse_boundary(char *buff, bound_t *bou)
{

  char *data;    

  if(!strlen(buff))
    return 1;

  data = strtok(buff," \n");
  if(!data)
    return 1;
  if(data[0]!='\"'||data[strlen(data)-1]!='\"')
    return 1;
  strcpy(bou->name,data);

  data = strtok(NULL," \n");
  if(!data)
    return 1;
  bou->order=atoi(data);

  data = strtok(NULL," \n");
  if(!data)
    return 1;
  bou->kind=atoi(data);

  list_init(&bou->nodeL,sizeof(int),cmp_int);
  list_init(&bou->elemsL,sizeof(int),cmp_int);

  return 0;
}

/*****************************************************************************************************/

int cmp_bou(void *a, void *b){

  if( ((bound_t*)a)->order>((bound_t*)b)->order){
    return 1;
  }else if( ((bound_t*)a)->order==((bound_t*)b)->order){
    return 0;
  }else{
    return -1;
  }
}


int cmp_time(void *a, void *b){

  if( ((tcontrol_t*)a)->tf>((tcontrol_t*)b)->tf){
    return 1;
  }else if( ((tcontrol_t*)a)->tf==((tcontrol_t*)b)->tf){
    return 0;
  }else{
    return -1;
  }
}

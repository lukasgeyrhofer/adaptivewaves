#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>


char filename_2d_in[128],filename_1d_out[128],filename_u[128];

int icount,dcount;
int *ival;
double *dval;

double dx = -1;
int space,space0;

double **c2;
double *c;
double *u;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}





void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfiles = 0;
  while((c = getopt(argn, argv,"i:u:o:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename_2d_in,optarg);
		haveinfiles += 1;
		break;
      case 'u':	strcpy(filename_u,optarg);
		haveinfiles += 1;
		break;
      case 'o':	strcpy(filename_1d_out,optarg);
		haveinfiles += 1;
		break;
    }
  }
  if(haveinfiles <= 2)print_error("need input files! use option -i and -u");
}


int read_2d_data() {
  FILE *fp;
  int i;
  fp=fopen(filename_2d_in,"rb");

  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(icount >= 3) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(&ival[0],icount,sizeof(int),fp);
    if(ival[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
    ival[2] = 0;
  }else{
    print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
  }
  if(dcount > 0) {
    dval = (double*)malloc(dcount*sizeof(double));
    fread(&dval[0],dcount,sizeof(double),fp);
  }
  
  c2=(double**)malloc(space*sizeof(double*));
  for(i=0;i<space;i++) {
    c2[i] = (double*)malloc(space*sizeof(double));
    fread(&c2[i][0],space,sizeof(double),fp);
  }
  
  fclose(fp);
}

int read_u_data() {
  FILE *fp;
  int icount_u,dcount_u;
  int *ival_u;
  double *dval_u;
  
  double dx_u;
  int space_u,space0_u;
  
  fp=fopen(filename_u,"rb");

  if(fp == NULL)print_error("input file not found!");
  fread(&icount_u,1,sizeof(int),fp);
  fread(&dcount_u,1,sizeof(int),fp);
  
  fread(&dx_u,1,sizeof(double),fp);
  fread(&space_u,1,sizeof(int),fp);
  fread(&space0_u,1,sizeof(int),fp);
  
  if(icount_u > 0) {
    ival_u = (int*)malloc(icount_u*sizeof(int));
    fread(&ival_u[0],icount_u,sizeof(int),fp);
    free(ival_u);
  }
  if(dcount_u > 0) {
    dval_u = (double*)malloc(dcount_u*sizeof(double));
    fread(&dval_u[0],dcount_u,sizeof(double),fp);
    free(dval_u);
  }
  
  if((space_u!=space)||(space0 != space0_u))print_error("lattice does not match");
  
  u=(double*)malloc(space*sizeof(double));
  fread(&u[0],space,sizeof(double),fp);

  fclose(fp);
}



void integrate_2d_data() {
  int i,j;
  c = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++) {
    c[i] = 0.;
    for(j=0;j<space;j++) c[i] += c2[i][j]*u[j];
    c[i] *= dx;
  }
}

int write_1d_data() {
  FILE *fp;
  
  fp = fopen(filename_1d_out,"wb");
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    fwrite(&ival[0],icount,sizeof(int),fp);
  }
  if(dcount > 0) {
    fwrite(&dval[0],dcount,sizeof(double),fp);
  }
  
  fwrite(&c[0],space,sizeof(double),fp);
  
  fclose(fp);
}



void get_1d_profiledata() {
  int i;
  double xg=0.,xc=0.,popsize=0.,xxc=0.;
  double x;
  double xmean,varx;
  
  for(i=0;i<space;i++) {
    x = (i-space0)*dx;
    xg += x*c[i]*u[i]*dx;
    xc += x*c[i]*dx;
    xxc+= x*x*c[i]*dx;
    popsize += c[i]*dx;
  }
  xmean = xc/popsize;
  varx = xxc/popsize - xmean*xmean;
  printf("%17.10e\t%17.10e\t%17.10e\t%17.10e\n",popsize,xg,xmean,varx);
}
  

void cleanup() {
  int i;
  for(i=0;i<space;i++)free(c2[i]);
  free(c2);
  free(c);
  free(u);
  free(ival);
  free(dval);
}



int main(int argn, char* argv[]) {
  int i;
  
  parsecommandline(argn,argv);
  read_2d_data();
  read_u_data();
  integrate_2d_data();
  write_1d_data();
  get_1d_profiledata();
  cleanup();
  return 0;
}

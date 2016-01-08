#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

double r = -1.;
double dx = -1;
int p_int,p_double;
int *p_int_v;
double *p_double_v;
char filename_in[128],filename_out[128],filename_u[128];
int space,space0;

int offset_output;
int outputinterval = 1;

double starttime;
double **c2;
double *u;

int precision = 6;
int latticeprecision;

int minindex = -1, maxindex = -1;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  int havefiles = 0;
  while((c = getopt(argn, argv,"i:o:u:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename_in,optarg);
		havefiles ++;
		break;
      case 'o': strcpy(filename_out,optarg);
		havefiles ++;
		break;
      case 'u':	strcpy(filename_u,optarg);
		havefiles++;
		break;
    }
  }
  if(havefiles < 3)print_error("need three files, options -i -o -u");
}


int read_c2_data() {
  FILE *fp;
  int i;
  fp=fopen(filename_in,"rb");

  fread(&p_int,1,sizeof(int),fp);
  fread(&p_double,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(p_int >= 3) {
    p_int_v = (int*)malloc(p_int*sizeof(int));
    fread(&p_int_v[0],p_int,sizeof(int),fp);
    if(p_int_v[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
  }else{
    print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
  }
  if(p_double > 0) {
    p_double_v = (double*)malloc(p_double*sizeof(double));
    fread(&p_double_v[0],p_double,sizeof(double),fp);
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
  int icount,dcount;
  int *ival;
  double *dval;
  
  double dx_u;
  int space_u,space0_u;
  
  fp=fopen(filename_u,"rb");

  if(fp == NULL)print_error("input file not found!");
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx_u,1,sizeof(double),fp);
  fread(&space_u,1,sizeof(int),fp);
  fread(&space0_u,1,sizeof(int),fp);
  
  if(icount > 0) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(&ival[0],p_int,sizeof(int),fp);
    free(ival);
  }
  if(dcount > 0) {
    dval = (double*)malloc(dcount*sizeof(double));
    fread(&dval[0],dcount,sizeof(double),fp);
    free(dval);
  }
  
  if((space_u!=space)||(space0 != space0_u))print_error("lattice does not match");
  
  u=(double*)malloc(space*sizeof(double));
  fread(&u[0],space,sizeof(double),fp);

  fclose(fp);
}



void adjust_constraint() {
  int i,j;
  double sum = 0.;
  double inv;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      sum += u[i]*u[j]*c2[i][j];
    }
  }
  sum *= dx*dx;
  fprintf(stderr,"%lf\n",sum);
  inv  = 1./sum;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      c2[i][j] *= inv;
    }
  }
}
  

void write_c2_data() {
  int i;
  FILE *fp;
  fp = fopen(filename_out,"wb");
  fwrite(&p_int,1,sizeof(int),fp);
  fwrite(&p_double,1,sizeof(int),fp);
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  fwrite(&p_int_v[0],p_int,sizeof(int),fp);
  fwrite(&p_double_v[0],p_double,sizeof(double),fp);
  for(i=0;i<space;i++)fwrite(&c2[i][0],space,sizeof(double),fp);
  fclose(fp);
}

void cleanup() {
  int i;
  for(i=0;i<space;i++)free(c2[i]);
  free(c2);
  free(u);
}

int main(int argn, char *argv[]) {
  parsecommandline(argn,argv);
  read_c2_data();
  read_u_data();
  adjust_constraint();
  write_c2_data();
  cleanup();
  return 0;
}

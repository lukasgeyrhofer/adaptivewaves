// adaptivewaves_readfile.c
// v0.01,	120907,		first (numbered) version
//				works with arbitrary format = (icount,dcount)
//				output:
//				  stderr: options
//				  stdout: density

// v0.02	121023		option -P: precision (=number of digits after comma) of density output 

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
char filename[128];
int space,space0;

int offset_output;
int outputinterval = 1;

double starttime;
double **c2;

int precision = 6;
int latticeprecision;


int filetype = 0;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"i:O:P:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename,optarg);
		if(strcmp(filename,"-") == 0) {
		  filetype = 2;
		}else{
		  filetype = 1;
		}
		break;
      case 'O':	outputinterval = atoi(optarg);
		break;
      case 'P':	precision=atoi(optarg);
		break;
    }
  }
  if(filetype == 0)print_error("need input file! use option -i");
}


int read_data() {
  FILE *fp;
  int i;
  if(filetype == 2) {
    fp = stdin;
  }else if(filetype == 1){
    fp=fopen(filename,"rb");
    if(fp == NULL)print_error("input file not found!");
  }
  fread(&p_int,1,sizeof(int),fp);
  fread(&p_double,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
//   printf("p_int = %d\n",p_int);
  if(p_int >= 3) {
    p_int_v = (int*)malloc(p_int*sizeof(int));
    fread(&p_int_v[0],p_int,sizeof(int),fp);
//     for(i=0;i<p_int;i++)printf("I[%d] = %d\n",i,p_int_v[i]);
// 	if(p_int_v[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
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

  if(filetype == 1)fclose(fp);
}

void print_options() {
  int i;
  fprintf(stderr,"# reading binary 2d data file for stochastic simulation of adaptive waves...\n");
  fprintf(stderr,"#   format        = (%d,%d)\n",p_int,p_double);
  fprintf(stderr,"#   dx            = %lf\n",dx);
  fprintf(stderr,"#   space         = %d\n",space);
  fprintf(stderr,"#   space0        = %d\n",space0);
  for(i=0;i<p_int;i++) {
    fprintf(stderr,"#   IVAL[%02d]      = %d\n",i,p_int_v[i]);
  }
  for(i=0;i<p_double;i++) {
    fprintf(stderr,"#   DVAL[%02d]      = %g\n",i,p_double_v[i]);
  }
  
  if(precision != 6) {
    fprintf(stderr,"#   outputprecision = %d\n",precision);
  }
  
}


void print_data() {
  int i,j;
  latticeprecision = (int)(ceil(log(1/dx)/log(10.)) + 1);
//   printf("latticeprecision = %d\n",latticeprecision);
//   exit(1);
  offset_output = space0 % outputinterval;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      if(((i+offset_output)%outputinterval == 0)&&((j+offset_output)%outputinterval == 0))printf("%*.*lf %*.*lf %*.*e\n",latticeprecision+5,latticeprecision,(i-space0)*dx,latticeprecision+5,latticeprecision,(j-space0)*dx,precision+7,precision,c2[i][j]);
    }
    if((i+offset_output)%outputinterval == 0)printf("\n");
  }
}


void cleanup() {
  int i;
  for(i=0;i<space;i++)free(c2[i]);
  free(c2);
}

int main(int argn, char *argv[]) {
  parsecommandline(argn,argv);
  read_data();
  print_options();
  print_data();
  cleanup();
  return 0;
}

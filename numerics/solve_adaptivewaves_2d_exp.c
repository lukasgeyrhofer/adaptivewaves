#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>


char u_infile[128];
char c_outfile[128];
int have_c_outfile = 0,have_u_infile = 0;


int space,space0;
double dx;

double *u;
double **c2,*c1;
double **c2_new;
double *x;


double *mutationkernel;
double mutationrate = 1e-5;
double mutationsigma = 1e-2;

double wavespeed;

double alpha = 1.;
int maxSteps = 1000;


int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"u:o:a:S:")) != -1){
    switch(c) {
      case 'u':	strcpy(u_infile,optarg);
		have_u_infile = 1;
		break;
      case 'o':	strcpy(c_outfile,optarg);
		have_c_outfile = 1;
		break;
      case 'a':	alpha = atof(optarg);
		break;
      case 'S':	maxSteps = atoi(optarg);
		break;
    }
  }
  if(have_u_infile == 0)print_error("need u* profile, option -u FILENAME");
}



int read_u_file(int importparameters) {
  int i;
  int icount,dcount;
  int *ival;
  double *dval;
  
  int u_space,u_space0;
  double u_dx;
  
  FILE *fp;
  
  fp = fopen(u_infile,"rb");
  if(fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    fread(&u_dx,sizeof(double),1,fp);
    fread(&u_space,sizeof(int),1,fp);
    fread(&u_space0,sizeof(int),1,fp);
    
    if(importparameters == 1) {
      space = u_space;
      space0 = u_space0;
      dx = u_dx;
    }else{
      if(( space != u_space) || (space0 != u_space0))print_error("lattice does not match");
    }
    
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(double),dcount,fp);
      if(dcount >= 4) {
	mutationrate = dval[2];
	wavespeed = dval[3];
      }
      free(dval);
    }
    u = (double*)malloc(space*sizeof(double));
    fread(u,sizeof(double),space,fp);
  }else{
    print_error("could not open u* file");
  }
}
    
  


double get_gaussian(double x,double sigma2) {
  return gsl_sf_exp(-x*x/(2.*sigma2))*0.398942/sqrt(sigma2);
}


void check_constraint() {
  int i,j;
  double sum = 0.,invsum;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      sum += u[i]*u[j]*c2[i][j];
    }
  }
  invsum = 1./sum;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      c2[i][j] *= invsum;
    }
  }
}

void initialize() {
  int i,j;
  double xx,yy;
  
  read_u_file(1);
  
  c2     = (double**)malloc(space*sizeof(double**);
  c2_new = (double**)malloc(space*sizeof(double**); 
  x      = (double*)malloc(space*sizeof(double*);
  mutationkernel = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++) {
    c2[i] = (double*)malloc(space*sizeof(double));
    c2_new[i] = (double*)malloc(space*sizeof(double));
    x[i] = (i-space0)*dx;
    mutationkernel[i] = gsl_sf_exp(-i*dx/mutationsigma)*mutationrate/mutationsigma*dx
  }
  
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      c2[i][j] = get_gaussian(x[i],wavespeed)*get_gaussian(x[j],wavespeed);
    }
  }
  
  check_constraint();
}





void iterate_c2(int timestep) {
  int i,j;
  
  
  
  
}



void cleanup() {
  int i;
  for(i=0;i<space;i++) {
    free(c2[i]);
    free(c2_new[i]);
  }
  free(c2);
  free(c2_new);
  free(u);
}



int main(int argn, char *argv[]) {
  int i;
  parsecommandline(argn,argv);
  initialize();
  
  
  
  
  
  cleanup();
}
  
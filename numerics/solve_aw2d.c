#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

char infilename[128],outfilename[128],ufilename[128];

double **c2,**c2new;
double *c1;
double *u;

double *coeff_m,*coeff_0,*coeff_p;

int space,space0;
double dx;

int maxsteps = 10000;
int adjuststeps = 1000;
double alpha = 1.;

double wavespeed;
double mutation_rate;
int cmdlinepar_wavespeed = 0;
int cmdlinepar_mutation_rate = 0;



double popsize;

// ===========================================================================================
// print_error
// ===========================================================================================

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


// ===========================================================================================
// parsecommandline
// ===========================================================================================

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0, haveufile = 0, haveoutfile = 0;
  while((c = getopt(argn, argv,"i:o:u:v:D:S:a:A:")) != -1){
    switch(c) {
      case 'i':	strcpy(infilename,optarg);
		haveinfile = 1;
		break;
      case 'o':	strcpy(outfilename,optarg);
		haveoutfile = 1;
		break;
      case 'u':	strcpy(ufilename,optarg);
		haveufile = 1;
		break;
      case 'v':	wavespeed = atof(optarg);
		cmdlinepar_wavespeed = 1;
		break;
      case 'D':	mutation_rate = atof(optarg);
		cmdlinepar_mutation_rate = 1;
		break;
      case 'S':	maxsteps = atoi(optarg);
		break;
      case 'a':	alpha = atof(optarg);
		break;
      case 'A':	adjuststeps = atoi(optarg);
		break;
    }
  }
  if(haveinfile + haveoutfile + haveufile < 3)print_error("need infile (-i FILENAME), outfile (-o FILENAME) and ufile (-u FILENAME)");
}


void read_c2_file() {
  int i;
  int icount,dcount;
  int *ival;
  double *dval;
  FILE *fpc2;
  
  fpc2 = fopen(infilename,"rb");
  if(fpc2 == NULL)print_error("could not open infile");
  fread(&icount,1,sizeof(int),fpc2);
  fread(&dcount,1,sizeof(int),fpc2);
  fread(&dx,1,sizeof(double),fpc2);
  fread(&space,1,sizeof(int),fpc2);
  fread(&space0,1,sizeof(int),fpc2);
  if(icount > 0) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(&ival[0],icount,sizeof(int),fpc2);
    if(icount < 3)print_error("infile contains no 2dim density - 1");
//     if(ival[2] != 1)print_error("infile contains no 2dim density - 2");
    free(ival);
  }
  if(dcount > 0) {
    dval = (double*)malloc(dcount*sizeof(double));
    fread(&dval[0],dcount,sizeof(double),fpc2);
    if(dcount >= 4) {
      if(cmdlinepar_wavespeed == 0)wavespeed = dval[3];
      if(cmdlinepar_mutation_rate == 0)mutation_rate = dval[2];
    }else{
      if((cmdlinepar_mutation_rate == 0)||(cmdlinepar_wavespeed==0))print_error("mutation rate and wavespeed not specified");
    }
    free(dval);
  }
  c2 = (double**)malloc(space*sizeof(double*));
  for(i=0;i<space;i++) {
    c2[i] = (double*)malloc(space*sizeof(double));
    fread(&c2[i][0],space,sizeof(double),fpc2);
  }
  fclose(fpc2);
}


void read_u_file() {
  int icount,dcount;
  int *ival;
  double *dval;
  FILE *fpu;
  double u_dx;
  int u_space,u_space0;
  
  fpu = fopen(ufilename,"rb");
  if(fpu == NULL)print_error("could not open ufile");
  fread(&icount,1,sizeof(int),fpu);
  fread(&dcount,1,sizeof(int),fpu);
  fread(&u_dx,1,sizeof(double),fpu);
  fread(&u_space,1,sizeof(int),fpu);
  fread(&u_space0,1,sizeof(int),fpu);
  if(icount > 0) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(&ival[0],icount,sizeof(int),fpu);
    free(ival);
  }
  if(dcount > 0){
    dval = (double*)malloc(dcount*sizeof(double));
    fread(&dval[0],dcount,sizeof(double),fpu);
    free(dval);
  }
  if((u_space != space) || (u_space0 != space0) || (u_dx/dx < .9) || (u_dx/dx > 1.1))print_error("lattice does not match");
  u = (double*)malloc(space*sizeof(double));
  fread(&u[0],space,sizeof(double),fpu);
  fclose(fpu);
}


void write_c2_file() {
  int i;
  int icount = 3,dcount=4;
  int ival[3];
  double dval[4];
  FILE *fpc2;
  
  fpc2 = fopen(outfilename,"wb");
  fwrite(&icount,1,sizeof(int),fpc2);
  fwrite(&dcount,1,sizeof(int),fpc2);
  fwrite(&dx,1,sizeof(double),fpc2);
  fwrite(&space,1,sizeof(int),fpc2);
  fwrite(&space0,1,sizeof(int),fpc2);
  ival[0] = 1;  ival[1] = 0;  ival[2] = 1;
  dval[0] = 0.; dval[1] = 0.; dval[2] = mutation_rate; dval[3] = wavespeed;
  fwrite(&ival[0],3,sizeof(int),fpc2);
  fwrite(&dval[0],4,sizeof(double),fpc2);
  for(i=0;i<space;i++) {
    fwrite(&c2[i][0],space,sizeof(double),fpc2);
  }
  fclose(fpc2);
}



void initialize() {
  int i;
  double x;
  
  coeff_m = (double*)malloc(space*sizeof(double));
  coeff_0 = (double*)malloc(space*sizeof(double));
  coeff_p = (double*)malloc(space*sizeof(double));
  
  c1      = (double*)malloc(space*sizeof(double));
  c2new   = (double**)malloc(space*sizeof(double*));
  
  for(i=0;i<space;i++) {
    c2new[i] = (double*)malloc(space*sizeof(double));
    
    x = (i-space0)*dx;
    coeff_m[i] = - 0.5*wavespeed/dx + mutation_rate/(dx*dx);
    coeff_0[i] = x - 4.*u[i] - 2.*mutation_rate/(dx*dx);
    coeff_p[i] = 0.5*wavespeed/dx + mutation_rate/(dx*dx);
    
    c2[i][0] = 0.;
    c2[i][space-1] = 0.;
    c2[0][i] = 0.;
    c2[space-1][i] = 0.;
  }
}

void adjust_constaint() {
  int i,j;
  double sum=0.,inv;
  for(i=0;i<space;i++)for(j=0;j<space;j++)sum+=u[i]*u[j]*c2[i][j];
  sum *= dx*dx;
  inv = 1./sum;
  for(i=0;i<space;i++)for(j=0;j<space;j++)c2[i][j]*=inv;
}


int iterate_c2(int step) {
  int i,j;
  double f,fc;
  
  
  for(i=1;i<space-1;i++) {
    
    c1[i] = 0.;
    for(j=0;j<space;j++) {
      c1[i] += c2[i][j]*u[j];
    }
    c1[i] *= dx;
    
    for(j=1;j<space-1;j++) {
      f  = coeff_m[i] * c2[i-1][j] + coeff_0[i] * c2[i][j] + coeff_p[i] * c2[i+1][j];
      f += coeff_m[j] * c2[i][j-1] + coeff_0[j] * c2[i][j] + coeff_p[j] * c2[i][j+1];
      fc = coeff_0[i] + coeff_0[j];
      if(i==j) {
	f  += 2.*c1[i]/dx;
	fc += 2.*u[i];
      }
      c2new[i][j] = c2[i][j] - alpha*f/fc;
    }
  }
  
  for(i=1;i<space-1;i++)memcpy(&c2[i][1],&c2new[i][1],(space-2)*sizeof(double));
}
  



void cleanup() {
  int i;
  for(i=0;i<space;i++) {
    free(c2[i]);
    free(c2new[i]);
  }
  free(c2);
  free(c2new);
  free(u);
  free(c1);
  free(coeff_m);
  free(coeff_0);
  free(coeff_p);
}


int main(int argn, char *argv[]) {
  int i;
  parsecommandline(argn,argv);
  read_c2_file();
  read_u_file();
  
  initialize();
  
  for(i=1;i<=maxsteps;i++) {
    iterate_c2(i);
    if(i%adjuststeps==0)adjust_constaint();
  }
  
  write_c2_file();
  cleanup();
  return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

double dx_1,dx_2;
int icount_1,dcount_1,icount_2,dcount_2;
int *ival_1,*ival_2;
double *dval_1,*dval_2;
char filename_1[128],filename_2[128];
int space_1,space0_1,space_2,space0_2;

int minindex = -1,maxindex = -1;

double **c2_1,**c2_2;

double *c1_1,*c1_2;

double compare = 0.;

int haveufile = 0;
char ufilename[128];
double *u;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfiles = 0;
  while((c = getopt(argn, argv,"i:I:m:M:u:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename_1,optarg);
		haveinfiles += 1;
		break;
      case 'I':	strcpy(filename_2,optarg);
		haveinfiles += 1;
		break;
      case 'm':	minindex = atoi(optarg);
		break;
      case 'M': maxindex = atoi(optarg);
		break;
      case 'u':	haveufile = 1;
		strcpy(ufilename,optarg);
		break;
		
    }
  }
  if(haveinfiles <= 1)print_error("need input files! use option -i and -I");
}


int read_data_1() {
  FILE *fp;
  int i;
  fp=fopen(filename_1,"rb");

  fread(&icount_1,1,sizeof(int),fp);
  fread(&dcount_1,1,sizeof(int),fp);
  
  fread(&dx_1,1,sizeof(double),fp);
  fread(&space_1,1,sizeof(int),fp);
  fread(&space0_1,1,sizeof(int),fp);
  
  if(icount_1 >= 3) {
    ival_1 = (int*)malloc(icount_1*sizeof(int));
    fread(&ival_1[0],icount_1,sizeof(int),fp);
    if(ival_1[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
    free(ival_1);
  }else{
    print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
  }
  if(dcount_1 > 0) {
    dval_1 = (double*)malloc(dcount_1*sizeof(double));
    fread(&dval_1[0],dcount_1,sizeof(double),fp);
    free(dval_1);
  }
  
  c2_1=(double**)malloc(space_1*sizeof(double*));
  for(i=0;i<space_1;i++) {
    c2_1[i] = (double*)malloc(space_1*sizeof(double));
    fread(&c2_1[i][0],space_1,sizeof(double),fp);
  }
  
  if(minindex < 0)minindex = 0;
  if(minindex > space_1-1)minindex = space_1-1;

  if(maxindex < 0)maxindex = space_1;
  if(maxindex > space_1-1)maxindex = space_1;
  
  fclose(fp);
}



int read_data_2() {
  FILE *fp;
  int i;
  fp=fopen(filename_2,"rb");

  fread(&icount_2,1,sizeof(int),fp);
  fread(&dcount_2,1,sizeof(int),fp);
  
  fread(&dx_2,1,sizeof(double),fp);
  fread(&space_2,1,sizeof(int),fp);
  fread(&space0_2,1,sizeof(int),fp);
  
  if(icount_2 >= 3) {
    ival_2 = (int*)malloc(icount_2*sizeof(int));
    fread(&ival_2[0],icount_2,sizeof(int),fp);
    if(ival_2[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
    free(ival_2);
  }else{
    print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
  }
  if(dcount_2 > 0) {
    dval_2 = (double*)malloc(dcount_2*sizeof(double));
    fread(&dval_2[0],dcount_2,sizeof(double),fp);
    free(dval_2);
  }
  
  c2_2=(double**)malloc(space_2*sizeof(double*));
  for(i=0;i<space_2;i++) {
    c2_2[i] = (double*)malloc(space_2*sizeof(double));
    fread(&c2_2[i][0],space_2,sizeof(double),fp);
  }

  fclose(fp);
}


void read_u_file() {
  FILE *fpu;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(ufilename,"rb");
  if(fpu != NULL) {
    fread(&icount,1,sizeof(int),fpu);
    fread(&dcount,1,sizeof(int),fpu);
    fread(&u_dx,1,sizeof(double),fpu);
    fread(&u_space,1,sizeof(int),fpu);
    fread(&u_space0,1,sizeof(int),fpu);
    if(icount > 0) {
      tmpi = (int*)malloc(icount*sizeof(int));
      fread(&tmpi[0],icount,sizeof(int),fpu);
      free(tmpi);
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpu);
      free(tmpd);
    }
    if (( space_1 != u_space) || (space0_1 != u_space0)) {
      print_error("lattice does not match");
    }
    u=(double*)malloc(space_1*sizeof(double));
    fread(&u[0],space_1,sizeof(double),fpu);
    fclose(fpu);
  }
}

double compare_c2() {
  int i,j;
  double retval = 0.;
  double curval1,curval2;
  double indexdiff2_inv = 1./(1.*(maxindex-minindex)*(maxindex-minindex));
  if((space_1 != space_2)||(space0_1 != space0_2))print_error("lattice does not match");
  for(i=minindex;i<maxindex;i++) {
    for(j=minindex;j<maxindex;j++) {
      if(c2_1[i][j] > 0.) {
	curval1 = (c2_1[i][j] - c2_2[i][j]);
	curval2 = (c2_1[i][j] - c2_2[i][j]);
	retval += curval1*curval2*indexdiff2_inv;
      }
    }
  }
  return retval;
}

void get_density_parameters(double compare) {
  int i,j;
  double xg1 = 0.,xc1 = 0.,popsize1 = 0.;
  double xg2 = 0.,xc2 = 0.,popsize2 = 0.;
  double x;
  
  c1_1 = (double*)malloc(space_1*sizeof(double));
  c1_2 = (double*)malloc(space_1*sizeof(double));
  
  for(i=0;i<space_1;i++) {
    c1_1[i] = 0.;
    c1_2[i] = 0.;
    for(j=0;j<space_1;j++) {
      c1_1[i] += c2_1[i][j]*u[j];
      c1_2[i] += c2_2[i][j]*u[j];
    }
    c1_1[i] *= dx_1;
    c1_2[i] *= dx_1;
    
    x = (i-space0_1)*dx_1;
    
    xg1 += x*c1_1[i]*u[i];
    xg2 += x*c1_2[i]*u[i];
    
    xc1 += x*c1_1[i];
    xc1 += x*c1_2[i];
    
    popsize1 += c1_1[i];
    popsize2 += c1_2[i];
  }
  xg1 *= dx_1;
  xc1 /= popsize1;
  xg2 *= dx_1;
  xc2 /= popsize2;
  popsize1 *= dx_1;
  popsize2 *= dx_1;
  
  printf("%17.10lf\t%17.10lf\t%17.10lf\t%17.10lf\t%17.10lf\t%17.10lf\t%17.10lf\n",compare,xg1,xc1,popsize1,xg2,xc2,popsize2);
  
  
}
  


void cleanup() {
  int i;
  for(i=0;i<space_1;i++) {
    free(c2_1[i]);
    free(c2_2[i]);
  }
  free(c2_1);
  free(c2_2);
  if(haveufile == 1) {
    free(c1_1);
    free(c1_2);
    free(u);
  }
}

int main(int argn, char *argv[]) {
  double comp;
  parsecommandline(argn,argv);
  read_data_1();
  read_data_2();
  
  comp = compare_c2();
  if(haveufile == 1) {
    read_u_file();
    get_density_parameters(comp);
  }else{
    printf("%14.10e\n",comp);
  }
  cleanup();
  return 0;
}

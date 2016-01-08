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
char filename_u[128],filename_c2[128];
int space,space0;

int offset_output;
int outputinterval = 1;

double starttime;
double **c2;
double *c1;

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
  int haveinfiles = 0;
  while((c = getopt(argn, argv,"i:u:m:M:O:P:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename_c2,optarg);
		haveinfiles ++;
		break;
      case 'u': strcpy(filename_u,optarg);
		haveinfiles ++;
		break;
      case 'm':	minindex = atoi(optarg);
		break;
      case 'M':	maxindex = atoi(optarg);
		break;
      case 'O':	outputinterval = atoi(optarg);
		break;
      case 'P':	precision=atoi(optarg);
		break;
    }
  }
  if(haveinfiles < 2)print_error("need both input files, options -i and -u");
}


int read_c2_data() {
  FILE *fp;
  int i;
  fp=fopen(filename_c2,"rb");

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
  fp=fopen(filename_u,"rb");
  
  int space_u,space0_u;
  double dx_u;

  if(fp == NULL)print_error("input file not found!");
  fread(&p_int,1,sizeof(int),fp);
  fread(&p_double,1,sizeof(int),fp);
  
  fread(&dx_u,1,sizeof(double),fp);
  fread(&space_u,1,sizeof(int),fp);
  fread(&space0_u,1,sizeof(int),fp);
  
  if((space!=space_u)||(space0 != space0_u))print_error("lattice does not match");
  
  if(p_int > 0) {
    p_int_v = (int*)malloc(p_int*sizeof(int));
    fread(&p_int_v[0],p_int,sizeof(int),fp);
  }
  if(p_double > 0) {
    p_double_v = (double*)malloc(p_double*sizeof(double));
    fread(&p_double_v[0],p_double,sizeof(double),fp);
  }
  
  u=(double*)malloc(space*sizeof(double));
  fread(&u[0],space,sizeof(double),fp);

  fclose(fp);
}


void get_c1() {
  int i,j;
  c1 = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++) {
    c1[i] = 0.;
//     printf("u[%05d] = %lf\n",i,u[i]);
    for(j=0;j<space;j++) {
      c1[i] += u[j]*c2[i][j];
    }
    c1[i] *= dx;
  }
}



void print_data() {
  int i,j;
  double covij,vari,varj;
  
  latticeprecision = (int)(ceil(log(1/dx)/log(10.)));
//   printf("latticeprecision = %d\n",latticeprecision);
//   exit(1);
  offset_output = space0 % outputinterval;
  for(i=minindex;i<=maxindex;i++) {
    for(j=minindex;j<=maxindex;j++) {
      if(((i+offset_output)%outputinterval == 0)&&((j+offset_output)%outputinterval == 0)) {
	covij = c2[i][j] - c1[i]*c1[j];
	vari = c2[i][i] - c1[i]*c1[i];
	varj = c2[j][j] - c1[j]*c1[j];
	printf("%*.*lf %*.*lf %*.*e %*.*e %*.*e %*.*e\n",	latticeprecision+5,latticeprecision,(i-space0)*dx,
								latticeprecision+5,latticeprecision,(j-space0)*dx,
								precision+7,precision,covij,
								precision+7,precision,covij/sqrt(vari*varj),
								precision+7,precision,c2[i][j],
								precision+7,precision,c1[i]*c1[j]);
      }
    }
    if((i+offset_output)%outputinterval == 0)printf("\n");
  }
}


double check_constraint() {
  int i,j;
  double sum = 0.;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      sum += c2[i][j] * u[i] * u[j];
    }
  }
  sum *= dx*dx;
  return sum;
}
      


void cleanup() {
  int i;
  for(i=0;i<space;i++)free(c2[i]);
  free(c2);
  free(c1);
}

int main(int argn, char *argv[]) {
  parsecommandline(argn,argv);
  
  read_c2_data();
  read_u_data();  
  get_c1();
  
  
  fprintf(stderr,"integrated c2 profile: %lf\n",check_constraint());

  if(minindex > space-1)minindex = space-1;
  if(maxindex > space-1)maxindex = space-1;
  if(minindex < 0)minindex=0;
  if(maxindex < 0)maxindex=space-1;

  
//   printf("minindex = %d, maxindex = %d\n",minindex,maxindex);
  
//   print_options();
  print_data();
  cleanup();
  return 0;
}

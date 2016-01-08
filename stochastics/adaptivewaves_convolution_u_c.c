#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

char ufilename[128];
char cfilename[128];

double *u;
double *c;

int space,space0;
double dx;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


int parsecommandline(int argn, char* argv[]) {
  char c;
  int havecinfile = 0, haveuinfile = 0;
  while((c = getopt(argn,argv,"c:u:")) != -1) {
    switch(c) {
      case 'c':	strcpy(cfilename,optarg);
		havecinfile = 1;
		break;
      case 'u':	strcpy(ufilename,optarg);
		haveuinfile = 1;
		break;
    }
  }
  if(havecinfile + haveuinfile <2)print_error("need both options -c and -u");
}

void read_c() {
  int i;
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fp;
  
  fp = fopen(cfilename,"rb");
  if(fp==NULL)print_error("could not open c-infile");
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    tmpi = (int*)malloc(icount*sizeof(int));
    fread(&tmpi[0],icount,sizeof(int),fp);
    free(tmpi);
  }
  
  if(dcount > 0) {
    tmpd = (double*)malloc(dcount*sizeof(double));
    fread(&tmpd[0],dcount,sizeof(double),fp);
    free(tmpd);
  }
  c = (double*)malloc(space*sizeof(double));
  fread(&c[0],space,sizeof(double),fp);
  fclose(fp);
}


void read_u() {
  int i;
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fp;
  int u_space,u_space0;
  double u_dx;
  
  fp = fopen(ufilename,"rb");
  if(fp==NULL)print_error("could not open c-infile");
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&u_dx,1,sizeof(double),fp);
  fread(&u_space,1,sizeof(int),fp);
  fread(&u_space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    tmpi = (int*)malloc(icount*sizeof(int));
    fread(&tmpi[0],icount,sizeof(int),fp);
    free(tmpi);
  }
  
  if(dcount > 0) {
    tmpd = (double*)malloc(dcount*sizeof(double));
    fread(&tmpd[0],dcount,sizeof(double),fp);
    free(tmpd);
  }
  

  if((u_space!=space)||(u_space0!=space0)||(u_dx/dx < .9)||(u_dx/dx > 1.1)) {
    print_error("lattice does not match");
  }
  
  u = (double*)malloc(space*sizeof(double));
  fread(&u[0],space,sizeof(double),fp);
  fclose(fp);
}


double convolution(int offset) {
  int i;
  double sum =0.;
  
  for(i=offset;i<space-offset;i++) {
    sum+=c[i]*u[i+offset];
  }
  
  return sum*dx;
}

void cleanup() {
  free(u);
  free(c);
}



int main(int argn, char* argv[]) {
  int i;
  double conv;
  
  parsecommandline(argn,argv);
  read_c();
  read_u();
  
  for(i=-50;i<50;i++) {
    conv = convolution(i);
    printf("%14.10lf %g\n",i*dx,conv);
  }    
  
  cleanup();
  return 0;
}


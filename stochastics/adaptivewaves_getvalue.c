#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

char infile[128];
int space,space0;
double dx;

double *c;

double x;

int exponential = 0;

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


int parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  int havexvalue = 0;
  while((c = getopt(argn, argv,"i:x:e")) != -1){
    switch(c) {
      case 'i':	strcpy(infile,optarg);
		haveinfile = 1;
		break;
      case 'x':	x = atof(optarg);
		havexvalue = 1;
		break;
      case 'e':	exponential = 1;
		break;
    }
  }
  if (haveinfile+havexvalue<2)print_error("need both infile (-i) and xvalue (-x)");
  return 0;
}



void read_c_file() {
  FILE *fpc;
  int icount,dcount;
  int *iparams;
  double *dparams;
  
  fpc = fopen(infile,"rb");
  if(fpc != NULL) {
    fread(&icount,1,sizeof(int),fpc);
    fread(&dcount,1,sizeof(int),fpc);
    fread(&dx,1,sizeof(double),fpc);
    fread(&space,1,sizeof(int),fpc);
    fread(&space0,1,sizeof(int),fpc);
    if(icount > 0) {
      iparams=(int*)malloc(icount*sizeof(int));
      fread(&iparams[0],icount,sizeof(int),fpc);
      free(iparams);
    }
    if(dcount > 0) {
      dparams = (double*)malloc(dcount*sizeof(double));
      fread(&dparams[0],dcount,sizeof(double),fpc);
      free(dparams);
    }
    c=(double*)malloc(space*sizeof(double));
    fread(&c[0],space,sizeof(double),fpc);
    fclose(fpc);
  }
}

double get_x_value(double xx) {
  int index;
  double x1,x2,y1,y2,yy;
  
  if((x<-space0*dx)||(x>(space-space0)*dx))print_error("xvalue not in lattice");
  index = (int)floor(xx/dx)+space0;
  
  x1 = (index-space0)*dx;
  x2 = (index+1-space0)*dx;
  y1 = c[index];
  y2 = c[index+1];
  
  yy = y1 + (y2-y1)*(xx-x1)/(x2-x1);
  return yy;
}

void cleanup() {
  free(c);
}

int main(int argn,char *argv[]) {
  double xv;
  parsecommandline(argn,argv);
  read_c_file();
  xv=get_x_value(x);
  if(exponential == 1) {
    printf("%14.6e\n",xv);
  }else{
    printf("%14.10lf\n",xv);
  }
  cleanup();
  return 0;
}

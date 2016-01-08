#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>


double dx = 0.01;
int space = 10001;
int space0 = 5000;
double N = 1000.;
double zeropos = 0.5;

char outfilename[128];

double *u;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}




void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveoutfile = 0;
  while((c = getopt(argn, argv,"o:N:z:d:s:")) != -1){
    switch(c) {
      case 'o':	strcpy(outfilename,optarg);
		haveoutfile= 1;
		break;
      case 'N':	N=atof(optarg);
		break;
      case 'z':	zeropos = atof(optarg);
		break;
      case 'd':	dx = atof(optarg);
		break;
      case 's':	space = atoi(optarg);
		break;
    }
  }
  if(haveoutfile==0)print_error("option -o necessary!");
}

int initialize() {
  int i;
  space0 = zeropos*space;
  u=(double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++)u[i] = 1./N;
}

void cleanup() {
  free(u);
}


int write_u() {
  FILE *fp;
  int i;
  int dcount = 0,icount = 0;
  
  fp  = fopen(outfilename,"wb");
  if(fp==NULL)print_error("could not open file");
  
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  fwrite(&u[0],space,sizeof(double),fp);
  
  fclose(fp);
}


int main(int argn, char* argv[] ){
  
  parsecommandline(argn,argv);

  initialize();
  write_u();
  cleanup();

  return 0;
}

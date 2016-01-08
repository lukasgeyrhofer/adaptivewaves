#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

int space=1000,space0=500;
double dx=1e-4;
double *c;

char filename[128];


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"s:z:d:o:")) != -1){
    switch(c) {
      case 'o':	strcpy(filename,optarg);
		break;
      case 's':	space = atoi(optarg);
		break;
      case 'z':	space0=atoi(optarg);
		break;
      case 'd': dx = atof(optarg);
                break;
    }
  }
}

int main(int argn, char* argv[]) {
    FILE *fp;
    int icount=0,dcount=0;
    parsecommandline(argn,argv);
    c = (double*)calloc(space,sizeof(double));
    c[space0]=1.;
    fp = fopen(filename,"wb");
    if(fp != NULL) {
        fwrite(&icount,sizeof(int),1,fp);
        fwrite(&dcount,sizeof(int),1,fp);
  
        fwrite(&dx,sizeof(double),1,fp);
        fwrite(&space,sizeof(int),1,fp);
        fwrite(&space0,sizeof(int),1,fp);
        fwrite(c,sizeof(double),space,fp);
        fclose(fp);
    }else{
        print_error("error opening file for writing");
    }

    return 0;
}


        

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>



char infile[128],outfile[128];
int havefiles = 0;



double *c,**c2;
int *iparams,*newiparams;
double *dparams;
int icount,dcount;

int space,space0;
double dx;

// ===========================================================================================
// print_error
// ===========================================================================================

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}



int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"i:o:")) != -1){
    switch(c) {
      case 'i':	strcpy(infile,optarg);
		havefiles += 1;
		break;
      case 'o':	strcpy(outfile,optarg);
		havefiles += 1;
		break;
    }
  }
  if (havefiles<2)print_error("need both infile and outfile");

  return 0;
}


void read_c_file() {
  FILE *fpc;
  int i;
  
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
    }
    if(dcount > 0) {
      dparams = (double*)malloc(dcount*sizeof(double));
      fread(&dparams[0],dcount,sizeof(double),fpc);
    }
    c=(double*)malloc(space*sizeof(double));
    fread(&c[0],space,sizeof(double),fpc);
    fclose(fpc);
  }
}


void convert() {
  int i,j;
  
  c2 =(double**)malloc(space*sizeof(double*));
  for(i=0;i<space;i++) {
    c2[i] = (double*)malloc(space*sizeof(double));
    for(j=0;j<space;j++) {
      c2[i][j] = c[i]*c[j];
    }
  }
  if(icount <= 2) {
    newiparams = (int*)malloc(3*sizeof(int));
    memcpy(&newiparams[0],&iparams[0],icount*sizeof(int));
    for(i=icount;i<=2;i++)iparams[i]=0;
    icount=3;
  }else{
    newiparams = iparams;
  }
  iparams[2] = 1;
}

void write_c2_file() {
  FILE *fpc;
  int i;
  fpc = fopen(outfile,"wb");
  fwrite(&icount,1,sizeof(int),fpc);
  fwrite(&dcount,1,sizeof(int),fpc);
  fwrite(&dx,1,sizeof(double),fpc);
  fwrite(&space,1,sizeof(int),fpc);
  fwrite(&space0,1,sizeof(int),fpc);
  fwrite(&newiparams[0],icount,sizeof(int),fpc);
  fwrite(&dparams[0],dcount,sizeof(double),fpc);
  for(i=0;i<space;i++)fwrite(&c2[i][0],space,sizeof(double),fpc);
  fclose(fpc);
}



int main(int argn, char *argv[]) {

  
  parsecommandline(argn,argv);
  
  read_c_file();
  convert();
  write_c2_file();
  
  

  return 0;
}


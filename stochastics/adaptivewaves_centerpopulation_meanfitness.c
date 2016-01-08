#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>


char infilename[128];
char outfilename[128];

double dx;
int space,space0;

double *dens;

int icount,dcount;
int *ivalues;
double *dvalues;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char* argv[]) {
  char c;
  int haveinfile = 0,haveoutfile = 0;
  while((c = getopt(argn,argv,"i:o:sbef:I:")) != -1) {
    switch(c) {
      case 'i':	strcpy(infilename,optarg);
		haveinfile = 1;
		break;
      case 'o': strcpy(outfilename,optarg);
		haveoutfile = 1;
		break;
    }
  }
  if(haveinfile<1)print_error("need input file, option -i");
  if(haveoutfile<1)print_error("need output file, option -o");
}


int read_file() {
  FILE *fp;
  
  fp = fopen(infilename,"rb");
  
  if(fp == NULL)print_error("could not read file");
  
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    ivalues = (int*)malloc(icount*sizeof(int));
    fread(&ivalues[0],icount,sizeof(int),fp);
  }
  
  if(dcount > 0) {
    dvalues = (double*)malloc(dcount*sizeof(double));
    fread(&dvalues[0],dcount,sizeof(double),fp);
  }
  
  dens = (double*)malloc(space*sizeof(double));
  fread(&dens[0],space,sizeof(double),fp);
  
  fclose(fp);
  return 0;
}


double get_meanfitness() {
  int i;
  double sum = 0.,xsum = 0.;
  for(i=0;i<space;i++) {
    sum += dens[i];
    xsum += dens[i]*(i-space0);
  }
  return xsum*dx/sum;
}


void adjust_space0(double meanfit) {
  space0 += (int)(meanfit/dx);
}

int write_file() {
  FILE *fp;
  
  fp = fopen(outfilename,"wb");
  
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  if(icount > 0) fwrite(&ivalues[0],icount,sizeof(int),fp);
  if(dcount > 0) fwrite(&dvalues[0],dcount,sizeof(double),fp);
  
  fwrite(&dens[0],space,sizeof(double),fp);
  
  fclose(fp);
  return 0;
}



void cleanup() {
  free(dens);
}


int main(int argn, char* argv[]) {
  double meanfit;
  
  parsecommandline(argn,argv);
  read_file();
  
  meanfit = get_meanfitness();
  
  printf("meanfit = %lf\n",meanfit);
  
  adjust_space0(meanfit);
  
  write_file();
  cleanup();
}
  
  

  


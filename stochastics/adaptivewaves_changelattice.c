#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>


char infilename[128];
char outfilename[128];

int smaller = 0;
int bigger = 0;
int extend = 0;
int newMin = 0;
int newMax = 0;

int factor = 1;

int interpolation = 1;

double *dens, *dens_new;

int space,space_new;
int space0,space0_new;
double dx,dx_new;


int dcount,icount;
double *dvalues;
int *ivalues;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char* argv[]) {
  char c;
  int haveinfile = 0,haveoutfile = 0;
  while((c = getopt(argn,argv,"i:o:sbemMf:I:")) != -1) {
    switch(c) {
      case 'i':	strcpy(infilename,optarg);
		haveinfile = 1;
		break;
      case 'o': strcpy(outfilename,optarg);
		haveoutfile = 1;
		break;
      case 's':	smaller = 1;
		bigger = 0;
		extend = 0;
		break;
      case 'b': bigger = 1;
		smaller = 0;
		extend = 0;
		break;
      case 'e': extend = 1;
		bigger = 0;
		smaller = 0;
		break;
      case 'm':	newMin = 1;
		break;
      case 'M':	newMax = 1;
		break;
      case 'f':	factor = atoi(optarg);
		break;
      case 'I':	if(strcmp("linear",optarg) == 0) {
		  interpolation = 1;
		}else{
		  interpolation = 0;
		}
		break;
    }
  }
  if(haveinfile<1)print_error("need input file, option -i");
  if(haveoutfile<1)print_error("need output file, option -o");
  if(smaller+bigger+extend+newMax+newMin<1)print_error("need to specify scale direction, option -b (bigger), option -s (smaller), option -e (extend), -option -m (new minimal index), option -M (new max index)");
  if(extend == 1)print_error("option -e not implemented yet");
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


int rescale() {
  int i,j;
  int offset = 0;
  double diff;
  if(smaller == 1) {
    offset = space0 % factor;
    space_new = space/factor;
    dx_new = factor * dx;
    dens_new = (double*)malloc(space_new*sizeof(double));
    j=0;
    for(i=0;i<space;i++) {
      if(i%factor == offset) {
	dens_new[j] = dens[i];
	if(i == space0)space0_new = j;
	j++;
      }
    }
  }
  
  if(bigger == 1) {
    space_new = space*factor+1;
    dens_new = (double*)malloc(space_new*sizeof(double));
    dx_new = dx/(1.*factor);
    if(interpolation == 1) {
      for(i=0;i<space-1;i++) {
	dens_new[i*factor] = dens[i];
	diff = (dens[i+1] - dens[i])/(1.*factor);
	for(j=1;j<factor;j++) {
	  dens_new[i*factor+j] = dens[i] + j*diff;
	}
      }
      dens_new[space_new-1] = dens[space-1];
    }
  }
  return 0;
}


void cropLattice() {
  if(newMax == 1) {
    if(factor < space) {
      space_new = factor;
      space0_new = space0;
      dx_new = dx;
      dens_new = (double*)malloc(space_new*sizeof(double));
      memcpy(&dens_new[0],&dens[0],space_new*sizeof(double));
    }
  }else if(newMin == 1) {
    if(factor < space) {
      space_new = space - factor;
      space0_new = space0 - factor;
      dx_new = dx;
      dens_new = (double*)malloc(space_new*sizeof(double));
      memcpy(&dens_new[0],&dens[factor],space_new*sizeof(double));
    }
  }
}



int write_file() {
  FILE *fp;
  
  fp = fopen(outfilename,"wb");
  
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx_new,1,sizeof(double),fp);
  fwrite(&space_new,1,sizeof(int),fp);
  fwrite(&space0_new,1,sizeof(int),fp);
  
  if(icount > 0) fwrite(&ivalues[0],icount,sizeof(int),fp);
  if(dcount > 0) fwrite(&dvalues[0],dcount,sizeof(double),fp);
  
  fwrite(&dens_new[0],space_new,sizeof(double),fp);
  
  fclose(fp);
  return 0;
}


void cleanup() {
  free(dens);
  free(dens_new);
  if(icount>0)free(ivalues);
  if(dcount>0)free(dvalues);
}


int main(int argn, char *argv[]) {
  

  parsecommandline(argn, argv);
  
  read_file();
  
  if(smaller + bigger + extend > 0) {
    rescale();
  }else if(newMax + newMin > 0) {
    cropLattice();
  }
  
  write_file();
  
  cleanup();


  return 0;
}

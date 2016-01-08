#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

char infilename[128];
char outfilename[128];

double dx;
int space,space0;


double dx_new;
int space_new,space0_new;

int icount,dcount;
int *ival;
double *dval;

int newspace_min = 0,newspace_max = -1;

double *dens, *dens_new;



int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile=0,haveoutfile=0;
  while((c = getopt(argn, argv,"i:o:m:M:")) != -1){
    switch(c) {
      case 'i':	strcpy(infilename,optarg);
		haveinfile = 1;
		break;
      case 'o':	strcpy(outfilename,optarg);
		haveoutfile = 1;
		break;
      case 'm': newspace_min = atoi(optarg);
		break;
      case 'M': newspace_max = atoi(optarg);
		break;
    }
  }
  if(haveinfile==0)print_error("need infile, option -i");
  if(haveoutfile==0)print_error("need outfile, option -o");
}


int read_data() {
  FILE *fp;

  fp=fopen(infilename,"rb");
  if(fp == NULL)print_error("input file not found!");
  
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(&ival[0],icount,sizeof(int),fp);
  }
  if(dcount > 0) {
    dval = (double*)malloc(dcount*sizeof(double));
    fread(&dval[0],dcount,sizeof(double),fp);
  }
  
  dens=(double*)malloc(space*sizeof(double));
  fread(&dens[0],space,sizeof(double),fp);

  fclose(fp);
}


int write_data() {
  FILE *fp;
  
  fp=fopen(outfilename,"wb");
  if(fp == NULL)print_error("could not open output file");
  
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space_new,1,sizeof(int),fp);
  fwrite(&space0_new,1,sizeof(int),fp);
  
  if(icount>0)fwrite(&ival[0],icount,sizeof(int),fp);
  if(dcount>0)fwrite(&dval[0],dcount,sizeof(double),fp);
  
  fwrite(&dens_new[0],space_new,sizeof(double),fp);
  fclose(fp);
  /*
  int i;
  for(i=0;i<space;i++) {
    printf("%d %g\n",i,dens[i]);
  }*/
}


int crop_data() {
  int i;
  int overlapp;
  if(newspace_max == -1)newspace_max = space;
  
  space_new = newspace_max - newspace_min;
  space0_new = space0 - newspace_min;
  
  
  dens_new = (double*)malloc(space_new*sizeof(double));
  for(i=0;i<space_new;i++)dens_new[i] = 0.;
  
  if(newspace_max < newspace_min)print_error("max < min");
  
  overlapp = 0;
  if(newspace_min < 0) {
    if(newspace_max > space)overlapp = space;
    if((newspace_max < space)&&(newspace_max > 0))overlapp = newspace_max;
    if(overlapp > 0)memcpy(&dens_new[-newspace_min],&dens[0],overlapp*sizeof(double));
  }else if(newspace_min<space) {
    if(newspace_max > space)overlapp = space - newspace_min;
    if(newspace_max < space)overlapp = space_new;
    if(overlapp > 0)memcpy(&dens_new[0],&dens[newspace_min],overlapp*sizeof(double));
  }
    
}


void cleanup() {
  free(dens);
  free(dens_new);
  if(icount>0)free(ival);
  if(dcount>0)free(dval);
}

int main(int argn, char *argv[]) {
  parsecommandline(argn,argv);

  
  read_data();
  crop_data();
  write_data();
  
  cleanup();
  return 0;
}

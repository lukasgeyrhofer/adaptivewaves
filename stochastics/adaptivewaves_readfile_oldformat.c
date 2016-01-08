#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

double r = -1.;
double dx = -1;
int version;
char filename[128];
int space,space0;

int offset_output;
int outputinterval = 1;

double starttime;

double *c;

int filetype = 0;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"i:O:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename,optarg);
		if(strcmp(filename,"-") == 0) {
		  filetype = 2;
		}else{
		  filetype = 1;
		}
		break;
      case 'O':	outputinterval = atoi(optarg);
		break;
    }
  }
  if(filetype == 0)print_error("need input file! use option -i");
}


int read_data() {
  FILE *fp;
  if(filetype == 2) {
    fp = stdin;
  }else if(filetype == 1){
    fp=fopen(filename,"rb");
    if(fp == NULL)print_error("input file not found!");
  }
  fread(&version,1,sizeof(int),fp);
  switch(version) {
    case 1:	fread(&r,1,sizeof(double),fp);
		fread(&dx,1,sizeof(double),fp);
		fread(&space,1,sizeof(int),fp);
		fread(&space0,1,sizeof(int),fp);
		fread(&starttime,1,sizeof(double),fp);
		break;
    default:	print_error("version of input file not correct");
		break;
  }
  c=(double*)malloc(space*sizeof(double));
  fread(&c[0],space,sizeof(double),fp);
  if(filetype == 1)fclose(fp);
}

void print_options() {
  fprintf(stderr,"reading binary data file for stochastic simulation of adaptive waves...\n");
  fprintf(stderr,"options:\n");
  fprintf(stderr,"  VERSION   = %d\n",version);
  fprintf(stderr,"  rho       = %6.4e\n",r);
  fprintf(stderr,"  dx        = %lf\n",dx);
  fprintf(stderr,"  space     = %d\n",space);
  fprintf(stderr,"  space0    = %d\n",space0);
  fprintf(stderr,"  starttime = %lf\n",starttime);
}


void print_data() {
  int i;
  offset_output = space0 % outputinterval;
  for(i=0;i<space;i++) {
    if((i+offset_output)%outputinterval == 0)printf("%lf %e\n",(i-space0)*dx,c[i]);
  }
}


void cleanup() {
  free(c);
}

int main(int argn, char *argv[]) {
  parsecommandline(argn,argv);
  read_data();
  print_options();
  print_data();
  cleanup();
  return 0;
}

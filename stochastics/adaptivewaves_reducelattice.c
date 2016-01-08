#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

char infile[128];
int haveinfile = 0;
char baseoutfile[128];
int haveoutfile = 0;

int space;
int space0;
double dx;

int rspace;
int rspace0;
double rdx;

int icount,dcount;
int *ival;
double *dval;

double *profile;
double *reducedprofile;

int iterations = 2;
int ratio = 2;
double currentratio = 1;



int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

int parsecommandline(int argn, char* argv[]) {
  char c;
  int pointpos;
  while((c = getopt(argn, argv,"i:I:R:o:")) != -1){
    switch(c) {
      case 'i':	strcpy(infile,optarg);
		haveinfile = 1;
		break;
      case 'I':	iterations = atoi(optarg);
		break;
      case 'R':	ratio = atoi(optarg);
		break;
      case 'o':	strcpy(baseoutfile,optarg);
		haveoutfile = 1;
		break;
    }
  }
  if(haveinfile == 0)print_error("need input file (-i FILENAME)");
  if(haveoutfile == 0) {
    pointpos = (int)(strrchr(infile,'.')-infile);
    strncpy(baseoutfile,infile,pointpos);
  } 
}


int read_file() {
  FILE *fp;
  fp=fopen(infile,"rb");
  if(fp == NULL)print_error("input file not found!");

  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(ival,sizeof(int),icount,fp);
  }
  if(dcount > 0) {
    dval = (double*)malloc(dcount*sizeof(double));
    fread(dval,sizeof(double),dcount,fp);
  }
  
  profile=(double*)malloc(space*sizeof(double));
  fread(profile,sizeof(double),space,fp);

  fclose(fp);
}


int reduce_profile(int cr) {
  int baseshift = cr/2;
  int i=0,j=0;
  double n=0.;
  while(i<space) {
    if(i%cr == cr-1) {
      reducedprofile[j] += profile[i];
      reducedprofile[j] /= n;
      n = 0.;
      j++;
    }else if (i%cr == 0) {
      reducedprofile[j] = profile[i];
    }else {
      reducedprofile[j] += profile[i];
    }
    n++;
    i++;
  }
}


int write_reduced_profile(int cr) {
  FILE *fp;
  char rfile[128];
  sprintf(rfile,"%s-%04d.conf",baseoutfile,cr);
  
  printf("  writing outfile '%s'\n",rfile);
  
  fp = fopen(rfile,"wb");
  if(fp != NULL) {
    fwrite(&icount,sizeof(int),1,fp);
    fwrite(&dcount,sizeof(int),1,fp);
    fwrite(&rdx,sizeof(double),1,fp);
    fwrite(&rspace,sizeof(int),1,fp);
    fwrite(&rspace0,sizeof(int),1,fp);
    if(icount>0)fwrite(ival,sizeof(int),icount,fp);
    if(dcount>0)fwrite(dval,sizeof(double),dcount,fp);
    fwrite(reducedprofile,sizeof(double),rspace,fp);
    fclose(fp);
  }else{
    print_error("could not open file for writing");
  }
}



int cleanup() {
  free(profile);
  if(icount > 0) free(ival);
  if(dcount > 0) free(dval);
}
 
  
int main(int argn, char *argv[]) {
  int i;
  parsecommandline(argn,argv);

  read_file();

  for(i=0;i<=iterations;i++) {
    printf("===================================\n");
    printf("  ratio = %.0lf\n",currentratio);
    rspace = (int)(space/currentratio);
    rspace0 = (int)(space0/currentratio);
    rdx = (dx*currentratio);
    
    printf("  new lattice\n    rspace  = %d\n    rspace0 = %d\n    rdx     = %e\n",rspace,rspace0,rdx);
    
    reducedprofile = (double*)malloc(rspace*sizeof(double));

    reduce_profile(currentratio);
    write_reduced_profile(currentratio);
    
    free(reducedprofile);

    currentratio *= ratio;
  }
  printf("===================================\n");

  cleanup();
  return 0;
}
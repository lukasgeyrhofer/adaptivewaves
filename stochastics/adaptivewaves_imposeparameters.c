#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

char ifilename[128],pfilename[128],ofilename[123];
int havefile[3] = {0,0,0};

int uselattice = 0;


int ispace,pspace;
int ispace0,pspace0;
double idx,pdx;

int icount,dcount;
int *ival;
double *dval;

double *profile;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  int havefiles = 0;
  while((c = getopt(argn, argv,"i:o:p:L")) != -1){
    switch(c) {
      case 'i':	strcpy(ifilename,optarg);
		havefile[0] = 1;
		break;
      case 'o':	strcpy(ofilename,optarg);
		havefile[1] = 1;
		break;
      case 'p':	strcpy(pfilename,optarg);
		havefile[2] = 1;
		break;
      case 'L':	uselattice = 1;
		break;
    }
  }
  if ( (havefile[0] == 0) || (havefile[1] == 0) || (havefile[2] == 0) ) print_error("need all three files -i -o -p");
}

int read_conf() {
  int ci,cd;
  int *itmp;
  double *dtmp;
  FILE *fp;
  fp = fopen(ifilename,"rb");
  if(fp != NULL) {
    fread(&ci,sizeof(int),1,fp);
    fread(&cd,sizeof(int),1,fp);

    fread(&idx,sizeof(double),1,fp);
    fread(&ispace,sizeof(int),1,fp);
    fread(&ispace0,sizeof(int),1,fp);
    
    if(ci>0) {
      itmp = (int*)malloc(ci*sizeof(int));
      fread(itmp,sizeof(int),ci,fp);
      free(itmp);
    }
    if(cd>0) {
      dtmp = (double*)malloc(cd*sizeof(double));
      fread(dtmp,sizeof(double),cd,fp);
      free(dtmp);
    }
    
    profile = (double*)malloc(ispace*sizeof(double));
    fread(profile,sizeof(double),ispace,fp);
    
    fclose(fp);
  }else{
    print_error("could not open file (-i)");
  }
}

int read_parameters() {
  FILE *fp;
  double *tmpprof;
  fp = fopen(pfilename,"rb");
  if(fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);

    fread(&pdx,sizeof(double),1,fp);
    fread(&pspace,sizeof(int),1,fp);
    fread(&pspace0,sizeof(int),1,fp);
    
    if(icount>0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
    }
    if(dcount>0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(double),dcount,fp);
    }
    
//     tmpprof = (double*)malloc(dspace*sizeof(double));
//     fread(tmpprof,sizeof(double),dspace,fp);
//     free(tmpprof);
    
    fclose(fp);
  }else{
    print_error("could not open file (-p)");
  }
  
  
}


int write_conf() {
  FILE *fp;
  
  fp = fopen(ofilename,"wb");
  if(fp != NULL) {
    fwrite(&icount,sizeof(int),1,fp);
    fwrite(&dcount,sizeof(int),1,fp);
    
    if(uselattice) {
      fwrite(&pdx,sizeof(double),1,fp);
      fwrite(&pspace,sizeof(int),1,fp);
      fwrite(&pspace0,sizeof(int),1,fp);
    }else{
      fwrite(&idx,sizeof(double),1,fp);
      fwrite(&ispace,sizeof(int),1,fp);
      fwrite(&ispace0,sizeof(int),1,fp);
    }
    
    if(icount>0)fwrite(ival,sizeof(int),icount,fp);
    if(dcount>0)fwrite(dval,sizeof(double),dcount,fp);
    
    if(uselattice) {
      fwrite(profile,sizeof(double),pspace,fp);
    }else{
      fwrite(profile,sizeof(double),ispace,fp);
    }
    
    fclose(fp);
  }else{
    print_error("could not open file (-o)");
  }
}


int main(int argn,char *argv[]) {
  parsecommandline(argn,argv);
  read_conf();
  read_parameters();
  write_conf();
}
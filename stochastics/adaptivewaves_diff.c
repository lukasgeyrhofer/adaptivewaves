#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

double dx;
int space,space0;

double *c1,*c2;
char filename1[128],filename2[128];

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"i:I:")) != -1){
    switch(c) {
      case 'i':	strcpy(filename1,optarg);
		break;
      case 'I':	strcpy(filename2,optarg);
		break;
    }
  }
}


int read_data1() {
  FILE *fp;
  int icount,dcount;
  int *ival;
  double *dval;
  
  fp=fopen(filename1,"rb");
  if(fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    
    fread(&dx,sizeof(double),1,fp);
    fread(&space,sizeof(int),1,fp);
    fread(&space0,sizeof(int),1,fp);
    
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(int),dcount,fp);
      free(dval);
    }
    
    c1=(double*)malloc(space*sizeof(double));
    fread(c1,sizeof(double),space,fp);

    fclose(fp);
  }else{
    print_error("input file not found!");
  }
}

int read_data2() {
  FILE *fp;
  int icount,dcount;
  int *ival;
  double *dval;
  
  fp=fopen(filename2,"rb");
  if(fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    
    fread(&dx,sizeof(double),1,fp);
    fread(&space,sizeof(int),1,fp);
    fread(&space0,sizeof(int),1,fp);
    
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(int),dcount,fp);
      free(dval);
    }
    
    c2=(double*)malloc(space*sizeof(double));
    fread(c2,sizeof(double),space,fp);

    fclose(fp);
  }else{
    print_error("input file not found!");
  }
}



int main(int argn, char *argv[] ){
  int i;
  double diff=0.;
  
  parsecommandline(argn,argv);
  read_data1();
  read_data2();

  
  for(i=0;i<space;i++) {
    diff += (c1[i] - c2[i])*(c1[i] - c2[i]);
  }
  
  printf("%.10e\n",diff);
  return 0;
}


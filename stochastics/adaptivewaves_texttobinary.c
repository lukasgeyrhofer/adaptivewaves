#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>



int iparametercount = 5,dparametercount = 5;

int icount = 0,dcount = 0;

double *dval;
int *ival;


int space,space0;
double dx;
double zeropos;

int haveinfile = 0;
char infilename[128];
char outfilename[128];

double *dens;

FILE *fpin;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveoutfile = 0;
  int latticeparameters = 0;
  double *dtmp;
  int *itmp;
  while((c = getopt(argn, argv,"i:o:I:D:s:d:z:")) != -1){
    switch(c) {
      case 'o': strcpy(outfilename,optarg);
		haveoutfile = 1;
		break;
      case 'i':	strcpy(infilename,optarg);
		haveinfile = 1;
		break;
      case 'I':	if(icount >= iparametercount) {
		  iparametercount*=2;
		  itmp=(int*)realloc(&ival[0],iparametercount*sizeof(int));
		  ival = itmp;
		}
		ival[icount] = atoi(optarg);
		icount++;
		break;
      case 'D': if(dcount >= dparametercount) {
		  dparametercount *= 2;
		  dtmp=(double*)realloc(&dval[0],dparametercount*sizeof(double));
		  dval = dtmp;
		}
		dval[dcount] = atof(optarg);
		dcount++;
		break;
      case 'd':	dx = atof(optarg);
		latticeparameters +=1;
		break;
      case 's':	space = atoi(optarg);
		latticeparameters +=1;
		break;
      case 'z':	space0 = atoi(optarg);	
		latticeparameters +=1;
		break;
    }
  }

//   printf("lp = %d\n",latticeparameters);
  if(latticeparameters<3)print_error("lattice not defined, options -d, -s, -z");
  if(haveoutfile<1)print_error("outfile not defined, option -o");
}

void initialize() {
  int i;
  
//   space0 = space*zeropos;
  dens = (double*)malloc(space*sizeof(double));
  
  if(haveinfile == 0) {
    fpin = stdin;
  }else{
    fpin = fopen(infilename,"r");
    if(fpin == NULL)print_error("could not open infile, option -i");
  }
    
  
  
}
  
  

void read_data() {
  int i;
  for(i=0;i<space;i++) {
    fscanf(fpin,"%lf",&dens[i]);
//     printf("dens[%05d] = %10.5lf\n",i,dens[i]);
  }
  if(haveinfile==1) {
    fclose(fpin);
  }
}



void write_data() {
  FILE *fp;
  
  fp = fopen(outfilename,"wb");
  if(fp == NULL)print_error("could not open outfile");
  
//   printf("icount = %d, dcount = %d\n",icount,dcount);
  
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  if(icount>0)fwrite(&ival[0],icount,sizeof(int),fp);
  if(dcount>0)fwrite(&dval[0],dcount,sizeof(double),fp);
  
  fwrite(&dens[0],space,sizeof(double),fp);
  fclose(fp);
}
  

  
  
void cleanup() {
  free(dens);
  free(ival);
  free(dval);
}
  

  
  
int main(int argn, char *argv[] ){
  
  int i;
  dval = (double*)malloc(dparametercount*sizeof(double));
  ival = (int*)malloc(iparametercount*sizeof(int));
  parsecommandline(argn,argv);
  
  initialize();
  
  read_data();
  
  write_data();
  
  cleanup();

  return 0;
}


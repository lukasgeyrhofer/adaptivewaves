// adaptivewaves_changebinaryfile.c
// v0.01,	120907,		first version
//				works with arbitrary format = (icount,dcount)
//				array index starts with 0!!
//				rescale possible with option -C (for closing at higher moments)
//				option -R removes all parameters (and overrides any other options)

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

double *dens;
double *dval;
int *ival;

void *tmp;

double dx;
int space,space0;

double densrescalefactor = 1.0;
int rescale = 0;

int dcount,icount;

int change_d = -1;
int change_i = -1;

int remove_parameters = 0;

double new_dval;
int new_ival;

char infile[128],outfile[128];

int havefiles = 0;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"i:o:D:I:V:C:R")) != -1){
    switch(c) {
      case 'i': strcpy(infile,optarg);
		havefiles++;
		break;
      case 'o': strcpy(outfile,optarg);
		havefiles++;
		break;
      case 'D': change_d = atoi(optarg);
		break;
      case 'I': change_i = atoi(optarg);
		break;
      case 'V': if(change_d>=0) {
		  new_dval = atof(optarg);
		}else if((change_i>=0)&&(change_d==-1)){
		  new_ival = atoi(optarg);
		}else{
		  print_error("could not set value (option -V)");
		}
		break;
      case 'R': remove_parameters = 1;
		break;
      case 'C': rescale = 1;
		densrescalefactor = atof(optarg);
		break;
    }
  }
  if(havefiles != 2) {
    print_error("infile and outfile not provided (options -i and -o)");
  }
  if((change_d>=0)&&(change_i>=0)) {
    print_error("can only change integer or double parameter!");
  }
}



int read_data() {
  FILE *fp;

  fp=fopen(infile,"rb");
  
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
  
  fp = fopen(outfile,"wb");
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    fwrite(&ival[0],icount,sizeof(int),fp);
  }
  if(dcount > 0) {
    fwrite(&dval[0],dcount,sizeof(double),fp);
  }
  
  fwrite(&dens[0],space,sizeof(double),fp);
  
  fclose(fp);
}

void adjust_i() {
  int i;
  if(change_i >= icount) {
    ival = (int*)realloc(&ival[0],(change_i+1)*sizeof(int));
    for(i=icount;i<=change_i;i++)ival[i]=0;
    icount = change_i+1;
  }
  ival[change_i] = new_ival;
}


void adjust_d() {
  int i;
  if(change_d >= dcount) {
    dval = (double*)realloc(&dval[0],(change_d+1)*sizeof(double));
    for(i=dcount;i<=change_d;i++)dval[i] = 0.;
    dcount = change_d+1;
  }
  dval[change_d] = new_dval;
}


void adjust_dens() {
  int i;
  for(i=0;i<space;i++) {
    dens[i] *= densrescalefactor;
  }
}


void print_data() {
  int i;
  fprintf(stderr,"# change of parameters in binary file\n");
  fprintf(stderr,"#   format   = (%d,%d)\n",icount,dcount);
  fprintf(stderr,"#   dx       = %lf\n",dx);
  fprintf(stderr,"#   space    = %d\n",space);
  fprintf(stderr,"#   space0   = %d\n",space0);
  for(i=0;i<icount;i++) {
    fprintf(stderr,"#   IVAL[%02d] = %d\n",i,ival[i]);
  }
  for(i=0;i<dcount;i++) {
    fprintf(stderr,"#   DVAL[%02d] = %g\n",i,dval[i]);
  }
  
  fprintf(stderr,"# changes:\n");
  if(remove_parameters == 0) {
    if(change_i>=0) {
      fprintf(stderr,"#   IVAL[%02d] = %d\n",change_i,new_ival);
    }
    if(change_d>=0) {
      fprintf(stderr,"#   DVAL[%02d] = %g\n",change_d,new_dval);
    }
    if(rescale>0) {
      fprintf(stderr,"#   RESCALE  = %g\n",densrescalefactor);
    }
  }else{
    fprintf(stderr,"#   REMOVE ALL PARAMETERS!\n");
  }
}



void cleanup() {
  free(dens);
  if(icount > 0)free(ival);
  if(dcount > 0)free(dval);
}

int main(int argn, char* argv[]) {
  
  
  parsecommandline(argn,argv);
  
  read_data();
  print_data();
  if(remove_parameters == 0) {
    if(change_i >= 0)adjust_i();
    if(change_d >= 0)adjust_d();
    if(rescale > 0)adjust_dens();
  }else{
    icount = 0;
    dcount = 0;
  }
  
  write_data();
  
  cleanup();


  return 0;
}

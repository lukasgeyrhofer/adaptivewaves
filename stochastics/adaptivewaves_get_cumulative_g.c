#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

char cinfilename[128];
char uinfilename[128];
char cgoutfilename[128];
int haveoutfile = 0;


double *c,*u,*g,*cg;
double dx;
int space,space0;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char* argv[]) {
  char c;
  int haveinfile = 0, haveufile = 0;
  while((c = getopt(argn,argv,"i:u:o:")) != -1) {
    switch(c) {
      case 'i':	strcpy(cinfilename,optarg);
		haveinfile = 1;
		break;
      case 'u':	strcpy(uinfilename,optarg);
		haveufile = 1;
		break;
      case 'o':	strcpy(cgoutfilename,optarg);
		haveoutfile = 1;
		break;
    }
  }
  if(haveinfile + haveufile < 2)print_error("need density- and u-files (-i & -u)");
}


void read_u() {
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fp;
  double u_dx;
  int u_space,u_space0;
  
  fp = fopen(uinfilename,"rb");
  if(fp==NULL)print_error("could not open u-infile");
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&u_dx,1,sizeof(double),fp);
  fread(&u_space,1,sizeof(int),fp);
  fread(&u_space0,1,sizeof(int),fp);
  
  if ((u_space!=space)||(u_space0 != space0)||(u_dx/dx > 1.1)||(u_dx/dx<0.9))print_error("lattice does not match");
  
  if(icount > 0) {
    tmpi=(int*)malloc(icount*sizeof(int));
    fread(&tmpi[0],icount,sizeof(int),fp);
    free(tmpi);
  }
  
  if(dcount > 0) {
    tmpd = (double*)malloc(dcount*sizeof(double));
    fread(&tmpd[0],dcount,sizeof(double),fp);
    free(tmpd);
  }
  u = (double*)malloc(space*sizeof(double));
  fread(&u[0],space,sizeof(double),fp);
  
  fclose(fp);
}


void read_c() {
  int i;
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fp;
  
  fp = fopen(cinfilename,"rb");
  if(fp==NULL)print_error("could not open c-infile");
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);
  
  if(icount > 0) {
    tmpi=(int*)malloc(icount*sizeof(int));
    fread(&tmpi[0],icount,sizeof(int),fp);
    free(tmpi);
  }
  
  if(dcount > 0) {
    tmpd = (double*)malloc(dcount*sizeof(double));
    fread(&tmpd[0],dcount,sizeof(double),fp);
    free(tmpd);
  }
  
  c = (double*)malloc(space*sizeof(double));
  fread(&c[0],space,sizeof(double),fp);
  fclose(fp);
}

void write_cg() {
  int icount = 0,dcount = 0;
  FILE *fpcg;
  fpcg = fopen(cgoutfilename,"wb");
  fwrite(&icount,1,sizeof(int),fpcg);
  fwrite(&dcount,1,sizeof(int),fpcg);
  fwrite(&dx,1,sizeof(double),fpcg);
  fwrite(&space,1,sizeof(int),fpcg);
  fwrite(&space0,1,sizeof(int),fpcg);
  fwrite(&cg[0],space,sizeof(double),fpcg);
  fclose(fpcg);
}

int main(int argn,char *argv[]) {
  int i,j;
  parsecommandline(argn,argv);
  read_c();
  read_u();
  
  g  = (double*)malloc(space*sizeof(double));
  cg = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++) {
    g[i] = u[i]*c[i];
  }
  for(i=0;i<space;i++) {
    cg[i] = 0.;
    for(j=i;j<space;j++) {
      cg[i] += g[j];
    }
    cg[i] *= dx;
  }
  
  if(haveoutfile == 1) {
    write_cg();
  }else{
    for(i=0;i<space;i++) {
      printf("%lf %e\n",(i-space0)*dx,cg[i]);
    }
  }
  return 0;
}
  
      
  
  
  
  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_blas.h>

int maxsteps   = 1000;
int outputstep = 100;
double epsilon = 1e-2;
int quiet = 0;

gsl_vector *c,*u,*x,*ones,*xx,*tmp,*prevc;
gsl_matrix *Lvstar;

gsl_vector_view c1,c2;
gsl_matrix *Lvstar1,*Lvstar2;

int haveufile = 0;
int havecfile = 0;
int havecoutfile = 0;
char uinfile[128];
char cinfile[128];
char coutfile[128];

int space;
int space0;
double dx;

double mutation_sigma;
double mutation_rate;
double wavespeed;

int override_mutation_sigma = 0;
int override_mutation_rate = 0;
int override_wavespeed = 0;

int updateconstraint = 0;

int stencil = 3;
int redblack = 0;

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"u:c:o:v:s:m:e:S:CqQO:FR")) != -1){
    switch(c) {
      case 'u':	strcpy(uinfile,optarg);
		haveufile = 1;
		break;
      case 'c':	strcpy(cinfile,optarg);
		havecfile = 1;
		break;
      case 'o':	strcpy(coutfile,optarg);
		havecoutfile = 1;
		break;
      case 'v':	wavespeed = atof(optarg);
		override_wavespeed = 1;
		break;
      case 's':	mutation_sigma = atof(optarg);
		override_mutation_sigma = 1;
		break;
      case 'm':	mutation_rate = atof(optarg);
		override_mutation_rate = 1;
		break;
      case 'S':	maxsteps = atoi(optarg);
		break;
      case 'O':	outputstep = atoi(optarg);
		break;
      case 'e':	epsilon = atof(optarg);
		break;
      case 'C':	updateconstraint = 1;
		break;
      case 'q': quiet = 2;
		break;
      case 'Q':	quiet = 1;
		break;
      case 'F':	stencil = 5;
		break;
      case 'R': redblack = 1;
		break;
    }
  }
  if(haveufile == 0) print_error("constraint needed, option -u FILENAME");
}


int read_c(int uselattice) {
  int icount,dcount;
  int *ival;
  double *dval;
  double cdx;
  int cspace,cspace0;
  FILE *fp;
  
  fp = fopen(cinfile,"rb");
  if (fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    fread(&cdx,sizeof(double),1,fp);
    fread(&cspace,sizeof(int),1,fp);
    fread(&cspace0,sizeof(int),1,fp);
    if (uselattice == 1) {
      space = cspace;
      space0 = cspace0;
      dx = cdx;
    }
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(double),dcount,fp);
      if (uselattice == 1) {
	if (override_mutation_rate == 0) mutation_rate = dval[0];
	if (override_wavespeed == 0) wavespeed = dval[1];
	if (override_mutation_sigma == 0) mutation_sigma = dval[2];
      }
      free(dval);
    }
    c = gsl_vector_alloc(space);
    fread(&c->data[0],sizeof(double),space,fp);
    fclose(fp);
  }
}


int read_u(int uselattice) {
  int icount,dcount;
  int *ival;
  double *dval;
  double udx;
  int uspace,uspace0;
  FILE *fp;
  
  fp = fopen(uinfile,"rb");
  if (fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    fread(&udx,sizeof(double),1,fp);
    fread(&uspace,sizeof(int),1,fp);
    fread(&uspace0,sizeof(int),1,fp);
    if (uselattice == 1) {
      space = uspace;
      space0 = uspace0;
      dx = udx;
    }
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(double),dcount,fp);
      if (uselattice == 1) {
	if (override_mutation_rate == 0) mutation_rate = dval[0];
	if (override_wavespeed == 0) wavespeed = dval[1];
	if (override_mutation_sigma == 0) mutation_sigma = dval[2];
      }
      free(dval);
    }
    u = gsl_vector_alloc(space);
    fread(&u->data[0],sizeof(double),space,fp);
    fclose(fp);
  }
}


int write_c() {
  int icount = 0,dcount = 3;
  FILE *fp;
  
  fp = fopen(coutfile,"wb");
  if (fp != NULL) {
    fwrite(&icount,sizeof(int),1,fp);
    fwrite(&dcount,sizeof(int),1,fp);
    fwrite(&dx,sizeof(double),1,fp);
    fwrite(&space,sizeof(int),1,fp);
    fwrite(&space0,sizeof(int),1,fp);
    fwrite(&mutation_rate,sizeof(double),1,fp);
    fwrite(&wavespeed,sizeof(double),1,fp);
    fwrite(&mutation_sigma,sizeof(double),1,fp);
    fwrite(&c->data[0],sizeof(double),space,fp);
    fclose(fp);
  }else{
    print_error("could not open cfile for writing");
  }
}
  

void update_constraint() {
  double tmp;
  gsl_blas_ddot(u,c,&tmp);
  gsl_blas_dscal(1./(tmp*dx),c);
}


int update_c_gsl(int timestep) {
  if(redblack == 0) {
    gsl_blas_dgemv(CblasNoTrans,epsilon,Lvstar,c,1,c);
  }else{
    gsl_blas_dgemv(CblasNoTrans,epsilon,Lvstar1,c,1,&c1.vector);
    gsl_blas_dgemv(CblasNoTrans,epsilon,Lvstar2,c,1,&c2.vector);
  }
}

int init_c_gaussian() {
  int i;
  double invvar2;
  invvar2 = dx*dx/(2*(wavespeed+mutation_rate*mutation_sigma));
  c = gsl_vector_alloc(space);
  for(i=0;i<space;i++) {
    gsl_vector_set(c,i,gsl_sf_exp(-(i-space0)*(i-space0)*invvar2));
  }
}
  


int initialize() {
  int i,j;
  gsl_vector *cpy;
  
  read_u(1);
  if(havecfile == 1) {
    read_c(0);
  }else{
    init_c_gaussian();
  }
  update_constraint();

  ones   = gsl_vector_alloc(space);
  x      = gsl_vector_alloc(space);
  xx     = gsl_vector_alloc(space);
  Lvstar = gsl_matrix_calloc(space,space);
  
  for(i=0;i<space;i++) {
    gsl_vector_set(ones,i,1.);
    gsl_vector_set(x,i,(i-space0)*dx);
    gsl_vector_set(xx,i,(i-space0)*(i-space0)*dx*dx);
    gsl_matrix_set(Lvstar,i,i,gsl_vector_get(x,i)-2*gsl_vector_get(u,i) - mutation_rate*gsl_sf_exp(-dx/mutation_sigma));
    
    if(stencil == 3) {
      if(i<space-1) gsl_matrix_set(Lvstar,i,i+1,gsl_matrix_get(Lvstar,i,i+1)+0.5*wavespeed/dx);
      if(i>0)       gsl_matrix_set(Lvstar,i,i-1,gsl_matrix_get(Lvstar,i,i-1)-0.5*wavespeed/dx);
    }else if(stencil == 5) {
      if(i<space-2) gsl_matrix_set(Lvstar,i,i+2,gsl_matrix_get(Lvstar,i,i+2)-1.*wavespeed/(12.*dx));
      if(i<space-1) gsl_matrix_set(Lvstar,i,i+1,gsl_matrix_get(Lvstar,i,i+1)+8.*wavespeed/(12.*dx));
      if(i>1)       gsl_matrix_set(Lvstar,i,i-1,gsl_matrix_get(Lvstar,i,i-1)-8.*wavespeed/(12.*dx));
      if(i>2)       gsl_matrix_set(Lvstar,i,i-2,gsl_matrix_get(Lvstar,i,i-2)-1.*wavespeed/(12.*dx));
    }else{
      print_error("stencil not defined");
    }
    
    for(j=1;j<=i;j++) {
      gsl_matrix_set(Lvstar,i,i-j,gsl_matrix_get(Lvstar,i,i-j) + mutation_rate*(1.-gsl_sf_exp(-dx/mutation_sigma))*gsl_sf_exp(-j*dx/mutation_sigma));
    }
  }
  
  if(redblack == 1) {
    if(space%2 != 0)print_error("need even 'space'!");
    c1 = gsl_vector_subvector_with_stride(c,0,2,space/2);
    c2 = gsl_vector_subvector_with_stride(c,1,2,space/2);
    cpy = gsl_vector_alloc(space);
    Lvstar1 = gsl_matrix_calloc(space/2,space);
    Lvstar2 = gsl_matrix_calloc(space/2,space);
    for(i=0;i<space;i+=2) {
      gsl_matrix_get_row(cpy,Lvstar ,i);
      gsl_matrix_set_row(    Lvstar1,i/2,cpy);
      gsl_matrix_get_row(cpy,Lvstar ,i+1);
      gsl_matrix_set_row(    Lvstar2,i/2,cpy);
    }
    gsl_vector_free(cpy);
  }
  

}




void print_c(int timestep) {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %lf %.10e\n",timestep*epsilon,gsl_vector_get(x,i),gsl_vector_get(c,i));
  }
  printf("\n");
}



void print_summary(int timestep) {
  int i;
  double xc,xxc,ps;
  gsl_blas_ddot(ones,c,&ps);
  gsl_blas_ddot(x,c,&xc);
  gsl_blas_ddot(xx,c,&xxc);
  fprintf(stderr,"%lf %.10e %.10e %.10e\n",timestep*epsilon,ps*dx,xc/ps,xxc/ps-(xc/ps)*(xc/ps));
}


void cleanup() {
  gsl_vector_free(c);
  gsl_vector_free(u);
  gsl_vector_free(x);
  gsl_vector_free(ones);
  gsl_vector_free(xx);
  gsl_vector_free(tmp);
  gsl_vector_free(prevc);
  gsl_matrix_free(Lvstar);
  if(redblack) {
    gsl_matrix_free(Lvstar1);
    gsl_matrix_free(Lvstar2);
  }
}


int main(int argn, char *argv[]) {
  int i;
  
  parsecommandline(argn,argv);
  initialize();
  
  if(quiet == 0)  print_c(0);
  if(quiet != 2)  print_summary(0);
  
  for(i=1;i<=maxsteps;i++) {
    update_c_gsl(i);
//     update_c_elements(i);
    if(updateconstraint == 1)update_constraint();
    if(i%outputstep == 0) {
      if(quiet == 0)  print_c(i);
      if(quiet != 2)  print_summary(i);
    }
  }
  
  if(havecoutfile == 1) {
    write_c();
  }
  
  cleanup();
  return 0;
}




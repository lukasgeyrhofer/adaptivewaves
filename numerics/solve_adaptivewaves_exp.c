#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_exp.h>

// lattice
int space,space0;
double dx;
int par_space=300,par_space0=100;
double par_dx=1e-3;

// mutation kernel
double mutation_rate = 1e-5;
double mutation_expdecay = 1e-2;


// wavespeed
double adaptationspeed = 1e-4;

// algorithm
double alpha = 1.;
int maxsteps = 100;
int quiet = 0;

// in- and output
int haveinfile = 0,haveoutfile = 0;
char u_infile[128],u_outfile[128];




// vectors & matrices

gsl_vector *u;
gsl_vector *fu;
gsl_matrix *jacobian;
gsl_matrix *invjacobian;
gsl_permutation *perm;


int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}



int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"s:z:d:D:M:v:i:o:a:S:qQ")) != -1){
    switch(c) {
      case 's':	par_space=atoi(optarg);
		break;
      case 'z':	par_space0=atoi(optarg);
		break;
      case 'd':	par_dx=atof(optarg);
		break;
      case 'i':	strcpy(u_infile,optarg);
		haveinfile = 1;
		break;
      case 'o':	strcpy(u_outfile,optarg);
		haveoutfile = 1;
		break;
      case 'D':	mutation_rate = atof(optarg);
		break;
      case 'M':	mutation_expdecay = atof(optarg);
		break;
      case 'v':	adaptationspeed = atof(optarg);
		break;
      case 'a':	alpha = atof(optarg);
		break;
      case 'S':	maxsteps = atoi(optarg);
		break;
      case 'q':	quiet = 2;
		break;
      case 'Q':	quiet = 1;
		break;
    }
  }
}

void read_u() {
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fpu;
  
  fpu = fopen(u_infile,"rb");
  if(fpu != NULL) {
    fread(&icount,1,sizeof(int),fpu);
    fread(&dcount,1,sizeof(int),fpu);
    fread(&dx,1,sizeof(double),fpu);
    fread(&space,1,sizeof(int),fpu);
    fread(&space0,1,sizeof(int),fpu);
    if(icount > 0) {
      tmpi = (int*)malloc(icount*sizeof(int));
      fread(&tmpi[0],icount,sizeof(int),fpu);
      free(tmpi);
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpu);
      free(tmpd);
    }
    
    u=gsl_vector_alloc(space);
    fread(&u->data[0],space,sizeof(double),fpu);
    fclose(fpu);
  }else{
    print_error("could not open u-file");
  }
}



void write_u() {
  int icount = 0,dcount = 0;
  int i;
  FILE *fp;
  
  fp = fopen(u_outfile,"wb");
  if(fp!=NULL) {
    fwrite(&icount,sizeof(int),1,fp);
    fwrite(&icount,sizeof(int),1,fp);
    fwrite(&dx,sizeof(double),1,fp);
    fwrite(&space,sizeof(int),1,fp);
    fwrite(&space0,sizeof(int),1,fp);
    fwrite(&u->data[0],sizeof(double),space,fp);
    fclose(fp);
  }
}

void print_u() {
  int i;
  for(i=0;i<space;i++) {
    printf("%10.6lf %.10e\n",(i-space0)*dx,gsl_vector_get(u,i));
  }
}

void approx_u() {
  int i;
  u=gsl_vector_alloc(space);
  for(i=0;i<=space0;i++)gsl_vector_set(u,i,1e-2*dx);
  for(i=space0+1;i<space;i++)gsl_vector_set(u,i,(i-space0)*dx*0.5);
}


void initialize() {
  int i;
  
  
  if(haveinfile) {
    read_u();
  }else{
    space=par_space;
    space0=par_space0;
    dx=par_dx;
    approx_u();
  }
  
  jacobian = gsl_matrix_calloc(space,space);
  invjacobian = gsl_matrix_calloc(space,space);
  perm = gsl_permutation_alloc(space);
  
  fu = gsl_vector_calloc(space);
}


void iterate_u_nonjacobian(int step) {
  int i,j;
  double f_tmp,fu_tmp;
  double lastu,nextu,curu;
  
  

  for(i=0;i<space;i++) {
    
//     printf("i=%d\n",i);
    if(i>0) {
      lastu = gsl_vector_get(u,i-1);
    }else{
      lastu = gsl_vector_get(u,i)*gsl_sf_exp(-dx/mutation_expdecay);
    }
    if(i<space-1) {
      nextu = gsl_vector_get(u,i+1);
    }else{
      nextu = gsl_vector_get(u,i)+0.5*dx;
    }
    curu = gsl_vector_get(u,i);

    f_tmp  = (adaptationspeed/(dx*dx)+0.5*((i-space0)*dx-mutation_rate+adaptationspeed/mutation_expdecay)/dx)*lastu - 2./dx*curu*lastu;
    f_tmp += (-2.*adaptationspeed/(dx*dx)+(i-space0)*dx/mutation_expdecay - 1. - 2./mutation_expdecay*curu)*curu;
    f_tmp += (adaptationspeed/(dx*dx)-0.5*((i-space0)*dx-mutation_rate+adaptationspeed/mutation_expdecay)/dx)*nextu + 2./dx*curu*nextu;
    fu_tmp = (-2*adaptationspeed/(dx*dx)+(i-space0)*dx/mutation_expdecay-1.-4/mutation_expdecay*curu+2./dx*nextu+2./dx*lastu);
    
    gsl_vector_set(fu,i,alpha*f_tmp/fu_tmp);
    
  }
  gsl_vector_add(u,fu);
}

void iterate_u(int step) {
  int i,j,s;
  double fu_tmp,nextu,j_tmp;
  
  fprintf(stderr,"%d\n",step);
  
  gsl_vector_set_zero(fu);
  gsl_matrix_set_zero(jacobian);
  
  for(i=0;i<space;i++) {
    fu_tmp = 0.;
    j_tmp = 0.;
  
    
    if(i>0) {
      fu_tmp =                     (adaptationspeed/(dx*dx)+0.5*((i-space0)*dx-mutation_rate+adaptationspeed/mutation_expdecay)/dx)*gsl_vector_get(u,i-1) - 2./dx*gsl_vector_get(u,i)*gsl_vector_get(u,i-1);
      gsl_matrix_set(jacobian,i,i-1,adaptationspeed/(dx*dx)+0.5*((i-space0)*dx-mutation_rate+adaptationspeed/mutation_expdecay)/dx                        - 2./dx*gsl_vector_get(u,i));
      j_tmp = -2./dx*gsl_vector_get(u,i-1);
    }
    fu_tmp += (-2.*adaptationspeed/(dx*dx)+(i-space0)*dx/mutation_expdecay - 1. - 2./mutation_expdecay*gsl_vector_get(u,i))*gsl_vector_get(u,i);
    j_tmp +=   -2.*adaptationspeed/(dx*dx)+(i-space0)*dx/mutation_expdecay - 1. - 4./mutation_expdecay*gsl_vector_get(u,i);
    if(i<space-1) {
      nextu = gsl_vector_get(u,i+1);
      gsl_matrix_set(jacobian,i,i+1,adaptationspeed/(dx*dx)-0.5*((i-space0)*dx-mutation_rate+adaptationspeed/mutation_expdecay)/dx                        + 2./dx*gsl_vector_get(u,i));
    }else{
      nextu = (space-1-space0)*dx*0.5;
    }

    fu_tmp +=                      (adaptationspeed/(dx*dx)-0.5*((i-space0)*dx-mutation_rate+adaptationspeed/mutation_expdecay)/dx)*nextu                 + 2./dx*gsl_vector_get(u,i)*nextu;
    j_tmp += 2./dx*nextu;
    
    gsl_vector_set(fu,i,fu_tmp);
    gsl_matrix_set(jacobian,i,i,j_tmp);
  }
  
  
  
  gsl_linalg_LU_decomp (jacobian, perm, &s);
  gsl_linalg_LU_invert (jacobian, perm, invjacobian);

  
  gsl_blas_dgemv(CblasNoTrans,-alpha,invjacobian,fu,1.0,u);
							// u = u - alpha*( J^-1 f_tmp)
  
}

void cleanup() {
  gsl_vector_free(u);
  gsl_vector_free(fu);
  gsl_matrix_free(jacobian);
  gsl_matrix_free(invjacobian);
  gsl_permutation_free(perm);
}


int main(int argn, char *argv[]) {
  int i;
  
  parsecommandline(argn,argv);
  initialize();
//   printf("in\n");
  for(i=1;i<=maxsteps;i++) {
    iterate_u_nonjacobian(i);
  }
  
  if(haveoutfile)write_u();
  if(quiet<2)print_u();
  
  cleanup();
  return 0;
}




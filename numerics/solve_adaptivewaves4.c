#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_exp.h>

// generic algorithm parameters
int maxsteps = 10000;
double alpha = 1.;

int adjust_step = 50;

int quiet = 0;


// switchc == 1 -> calculate density c(x)
// switchu == 1 -> calculate u*(x)
// switchx == 1 -> calculate spectral properties of linear operator for c(x)
int switchc = 0;
int switchu = 0;
int switchx = 0;

// names and flags for binary input- and output-file
char c_infile[128], c_outfile[128];
char u_infile[128], u_outfile[128];
int have_c_infile = 0, have_c_outfile = 0;
int have_u_infile = 0, have_u_outfile = 0;

// lattice
double dx = 1e-4;
int space = 10000,space0 = 5000;


// densities, etc...
// double *c,*cnew;
// double *u,*unew;

// double *c_coeff_p,*c_coeff_0,*c_coeff_m,*c_coeff_j;
// double *u_coeff_p,*u_coeff_0,*u_coeff_m,*u_coeff_j;


// model parameters
double wavespeed = 0.0001;


// mutations
int    mutation_type = 1; // 1 - diffusion, 2 - jumps, 3 - expdecay, 4 - external file
double mutation_rate = 0.00001;
int    mutation_jumpwidth = 20;
char   mutation_kernelfile[128];
double *mutation_kernel;
double mutation_expdecay = 0.001;
int    allowmutationoptions = 0;
double mutationalbias = 0.;
// initial conditions
int initialcond = 1; // 1 - known approximations, 2 - random, 3 - flat, 4 - gaussian (only c, unavailable for u)


// function pointers
// typedef int (*iteration_fp)(int timestep);
// iteration_fp iterate_u;
// iteration_fp iterate_c;


// matrix multiplication for iterative steps:

gsl_matrix *fu_jacobian;	// jacobian J_{ij} = \left(\frac{\partial F_i}{\partial u_j}\right)
gsl_matrix *const_jacobian;	// only diagonal elements depend on u[i] and have to be updated, in const_jacobian the constant part is stored
gsl_matrix *inverse_jacobian;	// J^-1
gsl_matrix *fu;			// matrix Fu in equation: F_i = (Fu_ij*u_j) - 2 u_i*u_i

gsl_vector *u, *unew;
gsl_vector *u2;
gsl_vector *f;
gsl_vector *f_u;
// gsl_vector *u_coeff_m, *u_coeff_0, *u_coeff_p;
gsl_vector *f_tmp;

gsl_vector *mutation_bc_correction;

gsl_vector *c, *cnew;
gsl_vector *df_dc;
gsl_vector *df_dc_inv;
gsl_matrix *fc;


gsl_permutation *perm;

// ===========================================================================================
// print_error
// ===========================================================================================

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}



int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"UCXi:o:u:I:S:a:s:d:z:m:M:v:D:QqA:")) != -1){
    switch(c) {
      case 'C':	if(switchu+switchx >= 1) {
		  print_error("can only use one of the options -U -C -X");
		}else{
		  switchc = 1;
		}
		break;
      case 'U':	if(switchc+switchx >= 1) {
		  print_error("can only use one of the options -U -C -X");
		}else{
		  switchu = 1;
		}
		break;
      case 'X': if(switchu+switchc >= 1) {
		  print_error("can only use one of the options -U -C -X");
		}else{
		  switchx = 1;
		}
		break;
      case 'i':	if(switchc == 1) {
		  strcpy(c_infile,optarg);
		  have_c_infile = 1;
		}else if(switchu == 1) {
		  strcpy(u_infile,optarg);
		  have_u_infile = 1;
		}else{
		  print_error("option -C or -U has to be specified first");
		}
		break;
      case 'o':	if(switchc == 1) {
		  strcpy(c_outfile,optarg);
		  have_c_outfile = 1;
		}else if(switchu == 1) {
		  strcpy(u_outfile,optarg);
		  have_u_outfile = 1;
		}else{
		  print_error("option -C or -U has to be specified first");
		}
		break;
      case 'u':	if(switchc+switchx == 1) {
		  strcpy(u_infile,optarg);
		  have_u_infile = 1;
		}else{
		  print_error("option -u requires option -C/-X");
		}
		break;
		
      case 'm':	if(strcmp(optarg,"diffusion") == 0) {
		  mutation_type = 1;
		  allowmutationoptions = 1;
		}else if(strcmp(optarg,"jump") == 0) {
		  mutation_type = 2;
		  allowmutationoptions = 1;
		}else if(strcmp(optarg,"expdecay") == 0) {
		  mutation_type = 3;
		  allowmutationoptions = 1;
		}else if(strcmp(optarg,"file") == 0) {
		  mutation_type = 4;
		  allowmutationoptions = 1;
		}else{
		  print_error("argument provided to -m option is not valid: choose from 'diffusion', 'jump', 'expdecay' or 'file'");
		}
		break;
      case 'M':	if(allowmutationoptions == 1) {
		  switch(mutation_type) {
		    case 1:	print_error("'diffusion' mutation type does not have additional parameters");
				break;
		    case 2:	mutation_jumpwidth = atoi(optarg);
				break;
		    case 3:	mutation_expdecay = atof(optarg);
				break;
		    case 4:	strcpy(mutation_kernelfile,optarg);
				break;
		  }
		}else{
		  print_error("mutation type is not specified yet, use option -m");
		}
		break;
		
      case 'I':	if(strcmp(optarg,"approx") == 0) {
		  initialcond = 1;
		}else if(strcmp(optarg,"random") == 0) {
		  initialcond = 2;
		}else if(strcmp(optarg,"flat") == 0) {
		  initialcond = 3;
		}else if(strcmp(optarg,"gaussian") == 0) {
		  if(switchc==1) {
		    initialcond = 4;
		  }else{
		    print_error("'gaussian' initial conditions only available for density c(x), use option -C");
		  }
		}else if(strcmp(optarg,"zero") == 0) {
		  if(switchu==1) {
		    initialcond = 4;
		  }else{
		    print_error("'zero' initial conditions only available for u(x), use option -U");
		  }
		}else{
		  print_error("argument provided to option -I is not valid: choose from 'approx', 'random', 'flat' or 'gaussian' (only mode -C)");
		}
		break;
		
      case 'S':	maxsteps = atoi(optarg);
		break;
      case 'a':	alpha = atof(optarg);
		break;

      case 'd':	dx = atof(optarg);
		break;
      case 's':	space = atoi(optarg);
		break;
      case 'z':	space0 = atoi(optarg);
		break;
		
      case 'v':	wavespeed = atof(optarg);
		break;
      case 'D':	mutation_rate = atof(optarg);
		break;
      case 'Q':	quiet = 1;
		break;
      case 'q':	quiet = 2;
		break;
      case 'A':	adjust_step = atoi(optarg);
		break;
    }
  }
  if(switchc + switchu +switchx == 0) {
    print_error("have to use at least one of the options, -C -U -X");
  }
  if(switchc+switchx == 1) {
    if(have_u_infile == 0)print_error("need u* file to calculate density profile c(x)");
  }
  return 0;
}


// ===========================================================================================
// adjust the density to fulfil the constraint \int dx u(x) c(x) = 1
// ===========================================================================================

void adjust_constraint(int timestep) {
  int i;
  double sum = 0.,invsum;
  gsl_blas_ddot(u,c,&sum);
  invsum = 1./(sum * dx);
  gsl_blas_dscal(invsum,c);
  printf("adjusting, timestep = %d, invsum = %17.10e\n",timestep,invsum);
}



// ===========================================================================================
// initial conditions
// ===========================================================================================

double get_gaussian(double pos, double sigma2) {
  return exp(-pos*pos/(2.*sigma2));
}

void initial_c_gaussian() {
  int i;
  double sigma2 = wavespeed - mutationalbias;
  for(i=0;i<space;i++){
    gsl_vector_set(c,i,get_gaussian((i-space0)*dx,sigma2));
  }
  adjust_constraint(0);
}

void initial_c_approx() {
//   int i;
//   if(mutation_type == 1) { // this approximation is only valid for diffusion type mutations ...
//     for(i=0;i<space;i++)c[i] = exp(-wavespeed * (i-space0) *dx / mutation_rate)*u->data[i];
//   }else{
//     initial_c_gaussian();
//   }
}

void initial_c_random() {
//   int i;
//   for(i=0;i<space;i++)c[i] = (dx*rand())/(1.*RAND_MAX);
}

void initial_c_flat() {
//   int i;
//   for(i=0;i<space;i++)c[i] = dx;
}

void initial_u_approx() {
  int i;
  double x;
  for(i=0;i<=space0;i++)u->data[i] = dx;
  for(i=space0+1;i<space;i++)u->data[i] = 0.5*(i-space0)*dx;
}

void initial_u_flat() {
  int i;
  for(i=0;i<space;i++)u->data[i] = dx;
}

void initial_u_random() {
  int i;
  for(i=0;i<space;i++)u->data[i] = (dx*rand())/(1.*RAND_MAX);
}


void initial_u_zero() {
  int i;
  for(i=0;i<space0;i++)u->data[i] = 0.;
  for(i=space0;i<space;i++)u->data[i] = 0.5*(i-space0)*dx;
}
    


// ===========================================================================================
// read external files
// ===========================================================================================

void read_u_file(int extractparameters) {
  FILE *fpu;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(u_infile,"rb");
  if(fpu != NULL) {
    fread(&icount,1,sizeof(int),fpu);
    fread(&dcount,1,sizeof(int),fpu);
    fread(&u_dx,1,sizeof(double),fpu);
    fread(&u_space,1,sizeof(int),fpu);
    fread(&u_space0,1,sizeof(int),fpu);
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
    if(extractparameters == 1) {
      dx = u_dx;
      space = u_space;
      space0 = u_space0;
    }else if (( space != u_space) || (space0 != u_space0)) {
      print_error("lattice does not match");
    }
    u=gsl_vector_alloc(space);
    fread(&u->data[0],space,sizeof(double),fpu);
    fclose(fpu);
  }else{
    print_error("could not open u-file");
  }
}
  
  
void read_c_file(int extractparameters) {
  FILE *fpc;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int c_space,c_space0;
  double c_dx;
  
  fpc = fopen(c_infile,"rb");
  if(fpc != NULL) {
    fread(&icount,1,sizeof(int),fpc);
    fread(&dcount,1,sizeof(int),fpc);
    fread(&c_dx,1,sizeof(double),fpc);
    fread(&c_space,1,sizeof(int),fpc);
    fread(&c_space0,1,sizeof(int),fpc);
    if(icount > 0) {
      tmpi = (int*)malloc(icount*sizeof(int));
      fread(&tmpi[0],icount,sizeof(int),fpc);
      free(tmpi);
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpc);
      free(tmpd);
    }
    if(extractparameters == 1) {
      dx = c_dx;
      space = c_space;
      space0 = c_space0;
    }else if (( space != c_space) || (space0 != c_space0)) {
      print_error("lattice does not match");
    }
    c=gsl_vector_alloc(space);
    fread(&c->data[0],space,sizeof(double),fpc);
    fclose(fpc);
  }else{
    print_error("could not open c-file");
  }
}




// ===========================================================================================
// write output files
// ===========================================================================================

void write_c_file() {
  FILE *fpc;
  int icount = 0,dcount = 4;
  double tmp_starttime = 0.,tmp_rho = 0.;
  fpc = fopen(c_outfile,"wb");
  fwrite(&icount,1,sizeof(int),fpc);
  fwrite(&dcount,1,sizeof(int),fpc);
  fwrite(&dx,1,sizeof(double),fpc);
  fwrite(&space,1,sizeof(int),fpc);
  fwrite(&space0,1,sizeof(int),fpc);
  fwrite(&tmp_rho,1,sizeof(double),fpc);
  fwrite(&tmp_starttime,1,sizeof(double),fpc);
  fwrite(&mutation_rate,1,sizeof(double),fpc);
  fwrite(&wavespeed,1,sizeof(double),fpc);
  fwrite(&c[0],space,sizeof(double),fpc);
  fclose(fpc);
}

void write_u_file() {
  FILE *fpu;
  int icount = 0,dcount = 0;
  fpu = fopen(u_outfile,"wb");
  fwrite(&icount,1,sizeof(int),fpu);
  fwrite(&dcount,1,sizeof(int),fpu);
  fwrite(&dx,1,sizeof(double),fpu);
  fwrite(&space,1,sizeof(int),fpu);
  fwrite(&space0,1,sizeof(int),fpu);
  fwrite(&u->data[0],space,sizeof(double),fpu);
  fclose(fpu);
}



// ===========================================================================================
// iteration procedures
// ===========================================================================================

int iterate_c(int timestep) {
  int i;
  // c'_i = c - alpha (df/dc)^-1_i (F c)_i
  // where fc_ij = (df/dc)^1_i F_ij
  gsl_blas_dgemv(CblasNoTrans,-alpha,fc,c,1,c);
}

int iterate_u(int timestep) {
  int i,j,s;
  double tmp;
  
  gsl_vector_memcpy(f_tmp,u);				// f_tmp[i] = u[i]
  gsl_vector_mul(f_tmp,u);				// f_tmp[i] = f_tmp[i]*u[i]
  gsl_blas_dgemv(CblasNoTrans,1.,fu,u,-2.,f_tmp);	// f_tmp = 1. * (SUM_j fu[ij]*u[j]) -2. *f_tmp[i]
  
  // special boundary conditions: mostly mutation correction
  gsl_vector_add(f_tmp,mutation_bc_correction);
  
  // u[space] does not exist, but we can approximate: u[space] = 0.5*(space-space0)*dx
  // u(x) ~ x/2 for x>xc>0
  
  // add diagonal elements to jacobian:
  gsl_matrix_memcpy(fu_jacobian,const_jacobian);
  for(i=0;i<space;i++) {
    tmp = gsl_matrix_get(fu_jacobian,i,i);
    gsl_matrix_set(fu_jacobian,i,i,tmp-4.*u->data[i]);	// J[ii] = -mutrate + x - 4*u[i]
  }
  // now jacobian can be inverted, with LU decomposition
  gsl_linalg_LU_decomp (fu_jacobian, perm, &s);
  gsl_linalg_LU_invert (fu_jacobian, perm, inverse_jacobian);

  
  gsl_blas_dgemv(CblasNoTrans,-alpha,inverse_jacobian,f_tmp,1.0,u);
							// u = u - alpha*( J^-1 f_tmp)
  
}

// ===========================================================================================
// 
// gives the stationary solution to our problem
// ===========================================================================================

void solve_spectrum() {
  int i,j;
  int index_0,index_l;
  double cur_evl = 1.e40,cur_ev0 = -1.;
  double cur_abs = 0.;
  
  printf("# solver starting...\n");
  
  gsl_vector_complex *eigenvalues;
  gsl_eigen_nonsymmv_workspace *w;
  gsl_complex ev0,evl;
  gsl_matrix_complex *evec;
  
  gsl_vector_view lr,li,zr,zi,evr,evi;
  
  gsl_vector_complex *largestev,*zeroev;
  largestev = gsl_vector_complex_alloc(space);

  w = gsl_eigen_nonsymmv_alloc(space);
  eigenvalues = gsl_vector_complex_alloc(space);
  evec = gsl_matrix_complex_alloc(space,space);
  
  printf("# allocation done...\n");
  
  gsl_eigen_nonsymmv(fc,eigenvalues,evec,w);
  
  printf("# solved for eigenvalues and eigenvectors...\n");
  
  evr=gsl_vector_complex_real(eigenvalues);
  evi=gsl_vector_complex_imag(eigenvalues);
  
  for(i=0;i<space;i++) {
    fprintf(stdout,"%4d\t%16.10e\t%16.10e\n",i,gsl_vector_get(&evr.vector,i),gsl_vector_get(&evi.vector,i));
    gsl_matrix_complex_get_row(largestev,evec,i);
    zr=gsl_vector_complex_real(largestev);
    zi=gsl_vector_complex_imag(largestev);
    
    fprintf(stderr,"%14.10lf",(i-space0)*dx);
    for(j=0;j<space;j++) {
      fprintf(stderr,"\t%16.10e\t%16.10e",gsl_vector_get(&zr.vector,j),gsl_vector_get(&zi.vector,j));
    }
    fprintf(stderr,"\n");
  }
  fflush(stdout);
  fflush(stderr);
  
  gsl_vector_complex_free(eigenvalues);
  gsl_eigen_nonsymmv_free(w);
  gsl_matrix_complex_free(evec);
  
  
  gsl_vector_complex_free(largestev);
  
  /*
  
  gsl_matrix *VV;
  gsl_vector *SS;
  
  VV = gsl_matrix_alloc(space,space);
  SS = gsl_vector_alloc(space);
  
  gsl_linalg_SV_decomp_jacobi(fc,VV,SS);
  
  
  free(VV);
  free(SS);
  
  exit(1);
  */
  
}



// ===========================================================================================
// initialize all variables
// ===========================================================================================

void initialize() {
  int i,j;
  double x,y,tmp;
  
  // ==================================================
  // initilization for C mode
  // ==================================================
  
  if (switchc+switchx == 1) {
    if(have_c_infile == 1) {
      read_c_file(1);
      read_u_file(0);
    }else{
      read_u_file(1);
      c = gsl_vector_alloc(space);
      switch(initialcond) {
	case 1:	initial_c_approx();
		break;
	case 2:	initial_c_random();
		break;
	case 3:	initial_c_flat();
		break;
	case 4:	initial_c_gaussian();
		break;
	default:print_error("something went wrong, option for initial condition not available");
      }
    }
    
    cnew = gsl_vector_alloc(space);
    fc = gsl_matrix_calloc(space,space);
    df_dc = gsl_vector_alloc(space);
    df_dc_inv = gsl_vector_alloc(space);
    
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      if(i>0) {
	gsl_matrix_set(fc,i,i-1,-0.5*wavespeed/dx);
      }
      gsl_matrix_set(fc,i,i,x-2.*u->data[i]);
      gsl_vector_set(df_dc,i,x-2.*u->data[i]);
      if(i<space-1) {
	gsl_matrix_set(fc,i,i+1,0.5*wavespeed/dx);
      }
      
      switch(mutation_type) {
	case 1:	if(i>0) {
		  gsl_matrix_set(fc,i,i-1,gsl_matrix_get(fc,i,i-1)+mutation_rate/(dx*dx));
		}
		tmp = gsl_matrix_get(fc,i,i);
		gsl_matrix_set(fc,i,i,tmp-2.*mutation_rate/(dx*dx));
		gsl_vector_set(df_dc,i,tmp-2*mutation_rate/(dx*dx));
		if(i<space-1) {
		  gsl_matrix_set(fc,i,i+1,gsl_matrix_get(fc,i,i+1)+mutation_rate/(dx*dx));
		}
		mutationalbias = 0.;
		break;
	case 2:	tmp = gsl_matrix_get(fc,i,i);
		gsl_matrix_set(fc,i,i,tmp-mutation_rate);
		gsl_vector_set(df_dc,i,tmp-mutation_rate);
		if(i>=mutation_jumpwidth) {
		  gsl_matrix_set(fc,i,i-mutation_jumpwidth,gsl_matrix_get(fc,i,i-mutation_jumpwidth)+mutation_rate);
		}
		mutationalbias = mutation_rate * mutation_jumpwidth * dx;
		break;
	case 3:	gsl_matrix_set(fc,i,i,gsl_matrix_get(fc,i,i)+mutation_rate*(dx/mutation_expdecay-1.));
		for(j=i+1;j<space;j++) {
		  gsl_matrix_set(fc,i,j,gsl_matrix_get(fc,i,j) + mutation_rate*dx/mutation_expdecay*gsl_sf_exp(-(j-i)*dx/mutation_expdecay));
		}
		mutationalbias = mutation_rate * mutation_expdecay;
		break;
	default:print_error("other mutation models not yet implemented");
      }
    }
  }
  // multiply the derivation (df/dc)^-1 directly into the matrix, as it is constant for every i
  // c profile is then calculated as
  // c'_i = c_i - alpha (df/dc)_i^-1 (F c)_i
  // with F including the mutational scheme (diffusion, jump, expdecay)
  // has to be done elementwise, as the (Fc) depends only linearly on c, so the inverse jacobian would be just F^-1, resulting in c' = c - alpha F^-1 F c = 0
  if(switchc == 1) {
    for(i=0;i<space;i++) {
      for(j=0;j<space;j++) {
	gsl_matrix_set(fc,i,j,gsl_matrix_get(fc,i,j)/gsl_vector_get(df_dc,i));
      }
    }
  }

  // ==================================================
  // initilization for U mode
  // ==================================================
  
  if (switchu == 1) {
    if (have_u_infile == 1) {
      read_u_file(1);
    }else{
      u = gsl_vector_alloc(space);
      switch(initialcond) {
	case 1:	initial_u_approx();
		break;
	case 2:	initial_u_random();
		break;
	case 3:	initial_u_flat();
		break;
	case 4:	initial_u_zero();
		break;
	default:print_error("something went wrong, option for initial condition not available");
      }
    }
    unew = gsl_vector_alloc(space);
    u2 = gsl_vector_alloc(space);
    f_tmp = gsl_vector_alloc(space);
    const_jacobian = gsl_matrix_calloc(space,space);
    fu_jacobian = gsl_matrix_calloc(space,space);
    inverse_jacobian = gsl_matrix_alloc(space,space);
    fu = gsl_matrix_calloc(space,space);
    perm = gsl_permutation_alloc(space);
    mutation_bc_correction = gsl_vector_calloc(space);
    
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      if(i>0) {
	gsl_matrix_set(const_jacobian,i,i-1,0.5*wavespeed/dx);
	gsl_matrix_set(fu,i,i-1,0.5*wavespeed/dx);
      }
      gsl_matrix_set(const_jacobian,i,i,x);
      gsl_matrix_set(fu,i,i,x);
      if(i<space-1) {
	gsl_matrix_set(const_jacobian,i,i+1,-0.5*wavespeed/dx);
	gsl_matrix_set(fu,i,i+1,-0.5*wavespeed/dx);
      }

      switch(mutation_type) {
	case 1:	// diffusion mutation scheme
		if(i>0) {
		  tmp = gsl_matrix_get(const_jacobian,i,i-1);
		  gsl_matrix_set(const_jacobian,i,i-1,tmp + mutation_rate/(dx*dx));
		  tmp = gsl_matrix_get(fu,i,i-1);
		  gsl_matrix_set(fu,i,i-1,tmp+mutation_rate/(dx*dx));
		}
		tmp = gsl_matrix_get(const_jacobian,i,i);
		gsl_matrix_set(const_jacobian,i,i,tmp - 2.*mutation_rate/(dx*dx));
		tmp = gsl_matrix_get(fu,i,i);
		gsl_matrix_set(fu,i,i,tmp - 2.*mutation_rate/(dx*dx));
		if(i<space-1) {
		  tmp = gsl_matrix_get(const_jacobian,i,i+1);
		  gsl_matrix_set(const_jacobian,i,i+1,tmp + mutation_rate/(dx*dx));
		  tmp = gsl_matrix_get(fu,i,i+1);
		  gsl_matrix_set(fu,i,i+1,tmp + mutation_rate/(dx*dx));
		}
		
		if(i==space-1) {
		  gsl_vector_set(mutation_bc_correction,space-1,0.5*(space-space0)*dx*(mutation_rate/(dx*dx) - 0.5*wavespeed/dx));
		}
		
		
		break;
		
	case 2:	// close to DFm model, delta-mutation scheme
		tmp = gsl_matrix_get(const_jacobian,i,i);
		gsl_matrix_set(const_jacobian,i,i, tmp - mutation_rate/dx);
		gsl_matrix_set(fu,i,i,tmp - mutation_rate/dx);
		
		if(i<space-mutation_jumpwidth) {
		  tmp = gsl_matrix_get(const_jacobian,i,i+mutation_jumpwidth);
		  gsl_matrix_set(const_jacobian,i,i+mutation_jumpwidth, tmp + mutation_rate/dx);
		  gsl_matrix_set(fu,i,i+mutation_jumpwidth, tmp + mutation_rate/dx);
		}else{
		  gsl_vector_set(mutation_bc_correction,i,0.5*(i-space0+mutation_jumpwidth)*mutation_rate);
		}
		
		if(i==space-1) {
		  gsl_vector_set(mutation_bc_correction,space-1,gsl_vector_get(mutation_bc_correction,space-1) - 0.25*(space-space0)*wavespeed);
		}
		break;
		
	case 3:	// positive exponential mutation scheme
		// i == j
		y=(space-i)*dx;
		gsl_vector_set(mutation_bc_correction,i,0.5*mutation_rate*(mutation_expdecay+(space-space0)*dx)*gsl_sf_exp(-y/mutation_expdecay));
		gsl_matrix_set(fu,i,i,gsl_matrix_get(fu,i,i)+mutation_rate*(dx/mutation_expdecay-1.));
		gsl_matrix_set(const_jacobian,i,i,gsl_matrix_get(const_jacobian,i,i)+mutation_rate*(dx/mutation_expdecay-1.));
		
		// i < j
		for(j=i+1;j<space;j++) {
		  gsl_matrix_set(fu,i,j,gsl_matrix_get(fu,i,j)+dx*mutation_rate/mutation_expdecay*gsl_sf_exp(-(j-i)*dx/mutation_expdecay));
		  gsl_matrix_set(const_jacobian,i,j,gsl_matrix_get(const_jacobian,i,j)+dx*mutation_rate/mutation_expdecay*gsl_sf_exp(-(j-i)*dx/mutation_expdecay));
		}
		
		if(i==space-1) {
		  gsl_vector_set(mutation_bc_correction,space-1,gsl_vector_get(mutation_bc_correction,space-1) - 0.25*(space-space0)*wavespeed);
		}
		break;

	default:print_error("other mutation models not yet implemented");
      }
    }
  }
}



void print_options() {
  printf("###################################################\n");
  printf("# solve profiles and weight-functions numerically #\n");
  printf("###################################################\n");  
  printf("#\n");
  printf("# mutationmodel     = %d\n",mutation_type);
  printf("# wavespeed         = %lf\n",wavespeed);
  printf("# mutationrate      = %lf\n",mutation_rate);
  if(mutation_type == 2) {
    printf("# option            = %d\n",mutation_jumpwidth);
  }
}

void print_u() {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %e\n",(i-space0)*dx,u->data[i]);
  }
}

void print_c() {
  int i;
  for(i=0;i<space;i++) {
    fprintf(stderr,"%lf %e\n",(i-space0)*dx,gsl_vector_get(c,i));
  }
}



// ===========================================================================================
// cleanup
// ===========================================================================================

void cleanup() {
  int i;
  if(switchc == 1) {
    gsl_vector_free(u);
  }else if (switchu == 1){
    gsl_vector_free(u);
    gsl_vector_free(unew);
    gsl_vector_free(u2);
    gsl_vector_free(f_tmp);
    gsl_vector_free(mutation_bc_correction);
    gsl_matrix_free(const_jacobian);
    gsl_matrix_free(fu_jacobian);
    gsl_matrix_free(inverse_jacobian);
    gsl_permutation_free(perm);
  }else if(switchx==1) {
    gsl_vector_free(c);
    gsl_vector_free(cnew);
    gsl_matrix_free(fc);
    gsl_vector_free(df_dc);
    gsl_vector_free(df_dc_inv);
    gsl_vector_free(u);
  }
  
  if(mutation_type == 4) {
    free(mutation_kernel);
  }
  
  
}



// ===========================================================================================
// main
// ===========================================================================================

int main(int argn, char *argv[]) {
  int i;
  
  
  parsecommandline(argn,argv);
  initialize();
  
  
  if(quiet<1)print_options();
  if(switchu == 1) {
    for(i=1;i<=maxsteps;i++) {
      iterate_u(i);
    }
    if(have_u_outfile == 1)write_u_file();
    if(quiet<2)print_u();
  }else if(switchx == 1) {
    solve_spectrum();
  }else if(switchc == 1) {
    for(i=1;i<=maxsteps;i++) {
      iterate_c(i);
      if(i%adjust_step == 0)adjust_constraint(i);
    }
    print_c();
    if(have_c_outfile == 1)write_c_file();
  }
  cleanup();
  return 0;
}



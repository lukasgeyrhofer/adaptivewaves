// solve_adaptivewaves
// numerical solution of the differential equations of u* and c in CRBWs
// lukas.geyrhofer@ds.mpg.de

// using newton's method at each lattice point to find solution

// v0.1		120604	first version
// v0.2		120619  minor updates
// v0.3		120704	flag for periodic/reflecting bc


// v0.4 	120820	modified code from oasis-problem (S = a \delta(x) -r) to adaptive waves (S= x -r)
//			using reduced variables

// v0.5		120906	new fileformat for binary files
// 			include output for densities

// v0.6		120910	can now read previous configurations,
//			-u -c output
//			-U -C input
//			option -m MODE: (MODE&1>0 - calc u), (MODE&2>0 - calc c)


// v0.6.01	120917	option -N, take no parameters, read all of them from previous c-file


// v0.7		130131	different mutation models
//			

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



double r = 0.0;
double v = 0.0;
double diffusionconstant = 1.0;
double lattice = 0.1;
double xmax = 500.;
int itersteps = 10000;
int verbose = 0;
int outputstep = 1000;
int startwithgaussian = 0;

double initialscale = 1.0;

int set_c_to_zero = 0;

int ccheck_interval = 1000;

int periodic_bc = 0;

int write_u_to_file = 0;
char ufilename[128];
int u_outputstep = 1;

int write_c_to_file = 0;
char cfilename[128];

int calc_u = 1;
int calc_c = 1;

double r_sqrt;

int space, space0;
double zeropos = 0.5;

int sizecontrainedby_space = 0;


double *c,*u,*unew,*cnew;

double *c_coeff_m1, *c_coeff_0, *c_coeff_p1;
double *u_coeff_m1, *u_coeff_0, *u_coeff_p1, *u2_coeff;
double *uder_coeff_0,*cder_const;

int randominitialcond = 0;
const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
unsigned long randseed = 0;


int noparameters = 0;

int read_u_from_file = 0;
char uinfilename[128];

int read_c_from_file = 0;
char cinfilename[128];

int initialcondition_c_from_approx_uEXP = 0;

// ===========================================================================================
// print_error
// ===========================================================================================

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


// ===========================================================================================
// parsecommandline
// ===========================================================================================

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn, argv,"r:V:d:M:z:S:O:U:u:o:C:Rs:D:Zc:A:I:m:NE")) != -1){
    switch(c) {
      case 'r':	r = atof(optarg);		// (global) death rate
		break;
      case 'V':	v = atof(optarg);		// convection velocity
		break;
      case 'D': diffusionconstant = atof(optarg);
		break;
      case 'z': space0 = atoi(optarg);
		break;
      case 'd':	lattice = atof(optarg);		// lattice constant (distance between two points in the discretized lattice)
		break;
      case 'M':	xmax = atof(optarg);		// size of the lattice in real units. symmetric around oasis, so lattice from -xmax to xmax
		break;
      case 'S':	itersteps = atoi(optarg);	// number of iterations to solve c and u*
		break;
      case 'O':	outputstep = atoi(optarg);	// output of convergence criteria every "outputstep" steps (criteria usually intergral over u/c)
		break;
      case 'u': strcpy(ufilename,optarg);	// write the solution u(x) to a file (name as parameter to "-u")
		write_u_to_file =1;
		break;
      case 'o': u_outputstep = atoi(optarg);	// only used with option "-u FILENAME"
		break;				// output only every U_OUTPUTSTEP lattice point for the u-file
						// useful when the numerical solution is on a finer grid than the stochastic simulation
      case 'I': ccheck_interval = atoi(optarg);
		break;
      case 'R': randominitialcond = 1;
		break;
      case 's': space = atoi(optarg);
		sizecontrainedby_space = 1;
		break;
      case 'Z': set_c_to_zero = 1;
		break;
      case 'c': strcpy(cfilename,optarg);
		write_c_to_file = 1;
		break;
      case 'A':	initialscale = atof(optarg);
		break;
      case 'm': if(strcmp(optarg,"u") == 0) {
		  calc_u = 1;
		  calc_c = 0;
		}else if(strcmp(optarg,"c") == 0) {
		  calc_u = 0;
		  calc_c = 1;
		}else{
		  calc_u = 1;
		  calc_c = 1;
		}
		break;
      case 'C':	read_c_from_file = 1;
		strcpy(cinfilename,optarg);
		break;
      case 'U':	read_u_from_file = 1;
		strcpy(uinfilename,optarg);
		break;
      case 'N':	noparameters = 1;
		break;
      case 'E': initialcondition_c_from_approx_uEXP = 1;
		break;
    }
  }
  if((noparameters==1)&&(read_c_from_file!=1))print_error("need c-file, option -C");
  randseed = time(NULL);
}


// ===========================================================================================
// get_gaussian
// ===========================================================================================

double get_gaussian(double x, double sigma2) {
  return exp(-x*x/sigma2);
}



// ===========================================================================================
// read_u
// ===========================================================================================

void read_u() {
  FILE *fpu;
  double ux;
  int i,tmp_count_params[2];
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(uinfilename,"rb");
  if(fpu != NULL) {
    fread(&tmp_count_params[0],2,sizeof(int),fpu);
    
    fread(&u_dx,1,sizeof(double),fpu);
    fread(&u_space,1,sizeof(int),fpu);
    fread(&u_space0,1,sizeof(int),fpu);
    
    if(tmp_count_params[0] > 0) {
      tmpi = (int*)malloc(tmp_count_params[0]*sizeof(int));
      fread(&tmpi[0],tmp_count_params[0],sizeof(int),fpu);
      free(tmpi);
    }
    if(tmp_count_params[1] > 0) {
      tmpd = (double*)malloc(tmp_count_params[1]*sizeof(double));
      fread(&tmpd[0],tmp_count_params[1],sizeof(double),fpu);
      free(tmpd);
    }
    
    if((space!= u_space) || (space0 != u_space0)) {
      print_error("reading u: lattice does not match");
    }
    u=(double*)malloc(space*sizeof(double));
    fread(&u[0],space,sizeof(double),fpu);
    
    fclose(fpu);
  }
}



// ===========================================================================================
// read_c
// ===========================================================================================

void read_c() {
  FILE *fpc;
  int i,tmp_count_params[2];
  int *tmpi;
  double *tmpd;
  int c_space,c_space0;
  double c_dx;
  
  fpc = fopen(cinfilename,"rb");
  if(fpc != NULL) {
    fread(&tmp_count_params[0],2,sizeof(int),fpc);
    
    fread(&c_dx,1,sizeof(double),fpc);
    fread(&c_space,1,sizeof(int),fpc);
    fread(&c_space0,1,sizeof(int),fpc);
    
    if(tmp_count_params[0] > 0) {
      tmpi = (int*)malloc(tmp_count_params[0]*sizeof(int));
      fread(&tmpi[0],tmp_count_params[0],sizeof(int),fpc);
      free(tmpi);
    }
    if(noparameters == 0) {
      if(tmp_count_params[1] > 0) {
	tmpd = (double*)malloc(tmp_count_params[1]*sizeof(double));
	fread(&tmpd[0],tmp_count_params[1],sizeof(double),fpc);
	free(tmpd);
      }
      if(sizecontrainedby_space == 0) {
	space = 2*(int)(xmax/lattice)+1;
      }else{
	xmax = (1.-zeropos)*lattice*space;
      }
      space0 = zeropos*space;
      
    }else{
      tmpd = (double*)malloc((tmp_count_params[1]-2)*sizeof(double));
      fread(&tmpd[0],2,sizeof(double),fpc);
      fread(&diffusionconstant,1,sizeof(double),fpc);
      fread(&v,1,sizeof(double),fpc);
      if(tmp_count_params[1]>4)fread(&tmpd[2],tmp_count_params[1]-4,sizeof(double),fpc);
      free(tmpd);
      
      space = c_space;
      space0 = c_space0;
      lattice = c_dx;
    }
    
    if((space!= c_space) || (space0 != c_space0)) {
      print_error("reading c: lattice does not match");
    }
    c=(double*)malloc(space*sizeof(double));
    fread(&c[0],space,sizeof(double),fpc);
    
    fclose(fpc);
  }
}

// ===========================================================================================
// initialize
// ===========================================================================================

void initialize() {
  int i,j;
  double x;
  double u2int_inv;
  double sigma2;
  double N_expectednovelocity,rho;

  double Dl2i;
  double v2li;

  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);


  
  
  
  if(read_c_from_file == 1) {
    read_c();
  }

  if(noparameters == 0) {
  }

  
  if(read_c_from_file == 0) { 
    c = (double*)malloc(space*sizeof(double));
    if(randominitialcond == 0) {
      for(i=0;i<space;i++)c[i] = exp(-v/diffusionconstant*fabs(i-space0)*lattice);
    }else{
      for(i=0;i<space;i++)c[i] = 0.001*gsl_rng_uniform(rg);
    }
  }
  
  printf("test\n");
  if(read_u_from_file == 1) {
    read_u();
  }else{
    u = (double*)malloc(space*sizeof(double));
    if(randominitialcond == 0) {
      for(i=0;i<space0;i++)u[i] = 0.0001;
      for(i=space0;i<space;i++)u[i]=0.5*(i-space0)*lattice;
    }else{
      for(i=0;i<space;i++)u[i] = 0.001*gsl_rng_uniform(rg);
    }
  }
  printf("testintermed\n");
  if(initialcondition_c_from_approx_uEXP == 1) {
    for(i=0;i<space;i++) {
      c[i] = u[i]*exp(-(i-space0)*lattice*v/diffusionconstant);
    }
  }
  
  printf("test2\n");
  

  cnew = (double*)malloc(space*sizeof(double));
  unew = (double*)malloc(space*sizeof(double));
  
  c_coeff_0  = (double*)malloc(space*sizeof(double));	// coefficients in the discretized ODE for c[i], c[i-1] and c[i+1]
  c_coeff_m1 = (double*)malloc(space*sizeof(double));
  c_coeff_p1 = (double*)malloc(space*sizeof(double));

  u_coeff_0  = (double*)malloc(space*sizeof(double));	// coefficients in the discretized ODE for u[i], u[i-1] and u[i+1]
  u_coeff_m1 = (double*)malloc(space*sizeof(double));
  u_coeff_p1 = (double*)malloc(space*sizeof(double));
  u2_coeff   = (double*)malloc(space*sizeof(double));
  
  uder_coeff_0 = (double*)malloc(space*sizeof(double));	// ODE differentiated with respect to u, for newton iteration
  cder_const = (double*)malloc(space*sizeof(double));	//				      c
  
  
  
//   sigma2 = xmax*xmax*initialscale*initialscale/(8.);					// (characteristic lenght)^2, only used with gaussian initial conditions
  
  // set initial values for u,c and all the coefficients
  
  Dl2i = diffusionconstant/(lattice*lattice);
  v2li = 0.5*v/lattice;
  
  for(i=0;i<space;i++) {
    x = (i-space0)*lattice;

    u_coeff_m1[i] = Dl2i + v2li;
    u_coeff_0[i] = -2.*Dl2i - r + x;
    u_coeff_p1[i] = Dl2i - v2li;
    
    if(x>-1.) {
      u2_coeff[i] = 2. + x;
    }else{
      u2_coeff[i] = 1.;
    }
    
    c_coeff_m1[i] = Dl2i - v2li;
    c_coeff_0[i] = -2.*Dl2i - r + x;
    c_coeff_p1[i] = Dl2i + v2li;
    
    uder_coeff_0[i] = -2.*Dl2i-r + x;
    
    cder_const[i] = uder_coeff_0[i];
    
    // periodic_bc == 0 means reflecting boundary condition,
    // no dependence on convection anymore, coefficient at boundaries is doubled or zero.
//     if(i==0) {
//       u_coeff_m1[i] = 0.;
//       u_coeff_p1[i] = 2.*Dl2i;
//       c_coeff_m1[i] = 0.;
//       c_coeff_p1[i] = 2.*Dl2i;
//     }
//     if(i==space-1) {
//       u_coeff_m1[i] = 2.*Dl2i;
//       u_coeff_p1[i] = 0.;
//       c_coeff_m1[i] = 2.*Dl2i;
//       c_coeff_p1[i] = 0.;
//     }
      
    
  }
}

// ===========================================================================================
// cleanup
// ===========================================================================================

void cleanup() {
  // free all the allocated memory again
  free(c);
  free(u);
  free(unew);
  free(cnew);
  free(c_coeff_0);
  free(c_coeff_m1);
  free(c_coeff_p1);
  free(u_coeff_0);
  free(u_coeff_m1);
  free(u_coeff_p1);
  free(u2_coeff);
  free(uder_coeff_0);
  free(cder_const);
}



// ===========================================================================================
// iterate_u
// ===========================================================================================

int iterate_u(int timestep) {
  int i,j;
  double f,fu;
  double u_int = 0.0;

  // periodic boundary condition, hence u[0] and u[space-1] (first and last points) have to be treated separately (to use correct indices)
  // by using the predefined coefficients "u_coeff_?" and not doing calculations every step, the speedup is something between 15-20%
  f =u_coeff_0[0]*u[0] + u_coeff_p1[0]*u[1] - 2.*u[0]*u[0];
  fu = uder_coeff_0[0]-4.*u[0];
  unew[0] = u[0] - f/fu;
#pragma omp parallel for 
  for(i=1;i<space-1;i++) {
    f = u_coeff_m1[i]*u[i-1] + u_coeff_0[i]*u[i] + u_coeff_p1[i]*u[i+1]-2.*u[i]*u[i];
    fu = uder_coeff_0[i]-4.*u[i];
    unew[i] = u[i] - f/fu;
  }
#pragma omp barrier
  
  f = u_coeff_m1[space-1] * u[space-2] + u_coeff_0[space-1]*u[space-1] + u_coeff_p1[space-1]*0.5*(space-space0)*lattice - 2.*u[space-1]*u[space-1];
  fu = uder_coeff_0[space-1]-4.*u[space-1];
  unew[space-1] = u[space-1] - f/fu;
  
  // after updating unew, everything copied to u
  memcpy(&u[0],&unew[0],space*sizeof(double));
  
  
  // output integral over u to see convergence of algorithm
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)u_int += u[i];
    u_int *= lattice;
    printf("%8d %.10e nan\n",timestep,u_int);
  }
  return 0;
}


// ===========================================================================================
// iterate_c
// ===========================================================================================

int iterate_c(int timestep) {
  int i,j;
  double f,fc;
  double c_int = 0.0;

  
  // same as for u, first and last treated separately
  
  if(set_c_to_zero == 1) {
    cnew[0] = 0.;
    cnew[space-1] = 0.;
  }else{
    f = c_coeff_m1[0] * c[space-1] + c_coeff_0[0]*c[0] + c_coeff_p1[0]*c[1];
    cnew[0] = c[0] - f/cder_const[0];
    
    f = c_coeff_m1[space-1] * c[space-2] + c_coeff_0[space-1]*c[space-1] + c_coeff_p1[space-1]*c[0];
    cnew[space-1] = c[space-1] - f/cder_const[space-1];
  }

  for(i=1;i<space-1;i++) {
    f = c_coeff_m1[i]*c[i-1] + c_coeff_0[i]*c[i] + c_coeff_p1[i]*c[i+1];
    cnew[i] = c[i] - f/cder_const[i];
  }
  
  memcpy(&c[0],&cnew[0],space*sizeof(double));

  // output population size to see convergence
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)c_int += c[i];
    c_int *= lattice;
    printf("%8d nan %.10e\n",timestep,c_int);
  }
  
  return 0;
}


// ===========================================================================================
// print_output
// ===========================================================================================

void print_output() {
  int i,j;
  double x;
  
  double popsize = 0.0;
  for(i=0;i<space;i++) {
    popsize += c[i];
  }
  popsize *= lattice;
  
  fprintf(stderr,"# N = %lf\n",popsize);
  
  
  if((calc_u == 1)&&(calc_c == 0)) {
    fprintf(stderr,"# x\t\t--\t\tu(x)\n");
  }else{
    fprintf(stderr,"# x\t\tc(x)\t\t\tu(x)\t\t\tc(x)*u(x)\n");
  }
  for(i=0;i<space;i++) {
    x = (i-space0)*lattice;
    if((calc_u == 1)&&(calc_c == 0)) {	// no solution for c present with option "-U", so second column is just "nan"
					// (second column is still written, so that u(x) is ALWAYS in third column)
      fprintf(stderr,"%lf\tnan\t%.10e\n",x,u[i]);
    }else{						// second column now c, fourth u*c
      fprintf(stderr,"%lf\t%.10e\t%.10e\t%.10e\n",x,c[i],u[i],c[i]*u[i]);
    }
  }
}



// ===========================================================================================
// write_u
// ===========================================================================================

void write_u() {
  // write u(x) solution to an additional file
  // which can then be read into the stochastic simulation
  // to have an arbitrary constraint
  
  FILE* fp;
  int i,j;
  
  int p_int    = 0; // number of integer parameters
  int p_double = 0; // number of double parameters
  
  double *u_write;
  int space_write, space0_write;
  double dx_write;
  
  int offset = space0 % u_outputstep;
  
  space_write = space/u_outputstep;
  dx_write = lattice*u_outputstep;
  u_write = (double*)malloc(space_write*sizeof(double));
  
  j=0;
  for(i=0;i<space;i++) {
    if(i%u_outputstep == offset) {
      u_write[j] = u[i];
      if(i==space0)space0_write = j;
      j++;
    }
  }
  
  fp = fopen(ufilename,"wb");
  if(fp!=NULL) {
    fwrite(&p_int,1,sizeof(int),fp);
    fwrite(&p_double,1,sizeof(int),fp);
    fwrite(&dx_write,1,sizeof(double),fp);
    fwrite(&space_write,1,sizeof(int),fp);
    fwrite(&space0_write,1,sizeof(int),fp);
    fwrite(&u_write[0],space,sizeof(double),fp);
  }
  fclose(fp);
  free(u_write);
}



// ===========================================================================================
// write_c
// ===========================================================================================
  
void write_c() {
  FILE* fp;
  int i,j;
  
  int p_int    = 0; // number of integer parameters
  int p_double = 4; // number of double parameters
  
  double tmp_rho = 0.0;
  double tmp_starttime = 0.0;
  
  double *c_write;
  int space_write, space0_write;
  double dx_write;
  
  int offset = space0 % u_outputstep;
  
  space_write = space/u_outputstep;
  dx_write = lattice*u_outputstep;
  c_write = (double*)malloc(space_write*sizeof(double));
  
  j=0;
  for(i=0;i<space;i++) {
    if(i%u_outputstep == offset) {
      c_write[j] = c[i];
      if(i==space0)space0_write = j;
      j++;
    }
  }
  
  fp = fopen(cfilename,"wb");
  if(fp!=NULL) {
    fwrite(&p_int,1,sizeof(int),fp);
    fwrite(&p_double,1,sizeof(int),fp);
    
    fwrite(&dx_write,1,sizeof(double),fp);
    fwrite(&space_write,1,sizeof(int),fp);
    fwrite(&space0_write,1,sizeof(int),fp);
    
    fwrite(&tmp_rho,1,sizeof(double),fp);
    fwrite(&tmp_starttime,1,sizeof(double),fp);
    fwrite(&diffusionconstant,1,sizeof(double),fp);
    fwrite(&v,1,sizeof(double),fp);
    
    fwrite(&c_write[0],space,sizeof(double),fp);
  }
  fclose(fp);
  free(c_write);
}
  
// ===========================================================================================
// consistency_check
// ===========================================================================================
  
void consistency_check(int timestep) {
  
  // integral (uc)dx should be 1
  
  int i,j;
  double popsize = 0.;
  double sum = 0.;
  
  for(i=0;i<space;i++) {
    popsize += c[i];
    sum += u[i]*c[i];
  }
  popsize *= lattice;
  sum *= lattice;
  
  // ODE for c(x) is only linear, so doesnt give correct prefactor
  // hence rescale c(x) to comply constraint:
  for(i=0;i<space;i++) {
    c[i] /= sum;
  }
  popsize /= sum; // population size has to be rescaled, too
}





// ===========================================================================================
// MAIN
// ===========================================================================================

int main(int argn, char *argv[]) {
  int i,j;
  parsecommandline(argn,argv);
  initialize();
  
  if(calc_u == 1) {
    for(i=0;i<itersteps;i++) iterate_u(i);
  }
  if(write_u_to_file == 1) {
    write_u();
  }
  
  if(calc_c == 1) {
    for(i=0;i<space;i++) {
      c_coeff_0[i] -= 2.*u[i];
      cder_const[i] -= 2.*u[i];
    }
    printf("\n");
    for(i=0;i<itersteps;i++) {
      iterate_c(i);
      if(i%ccheck_interval == 0)consistency_check(i);
    }
  }
  consistency_check(itersteps);
  print_output();

  if(write_c_to_file == 1) {
    write_c();
  }
  
  
  cleanup();
  return 0;
}


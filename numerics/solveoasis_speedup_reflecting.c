// solveoasis
// numerical solution of the differential equations of u* and c in CRBWs
// lukas.geyrhofer@ds.mpg.de

// using newton's method at each lattice point to find solution

// v0.1		120604	first version
// v0.2		120619  minor updates
// v0.3		120704	flag for periodic/reflecting bc


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>


double alpha = 0.1;
double a_effective = 0.0;
double r = 0.002;
double v = 0.0;
double diff = 1.0;
double lattice = 0.1;
double xmax = 500.;
int itersteps = 10000;
int verbose = 0;
int outputstep = 1000;
int startwithgaussian = 0;

int periodic_bc = 0;

int write_u_to_file = 0;
char ufilename[128];
int u_outputstep = 10;

int onlyU = 0;

double r_sqrt;

int space, space0;
double *c,*u,*unew,*cnew;

double *c_coeff_m1, *c_coeff_0, *c_coeff_p1;
double *u_coeff_m1, *u_coeff_0, *u_coeff_p1;
double *uder_coeff_0,*cder_const;



int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn, argv,"a:r:v:D:d:M:S:gO:Uu:o:P")) != -1){
    // parameter parsing, from library "unistd.h"
    // letters with a semicolon after them need an additional value
    // without semicolon just used as a flag (here: -V -U -g)
    // options are used as e.g. "-a 0.1" directly in the command-line
    // if option is not provided in command-line, default value is used (see definitions of variables above for default values)
    switch(c) {
      case 'a':	alpha = atof(optarg);		// growth rate alpha at oasis
		break;
      case 'r':	r = atof(optarg);		// (global) death rate
		break;
      case 'v':	v = atof(optarg);		// convection velocity
		break;
      case 'D':	diff = atof(optarg);		// diffusion coefficient
		break;
      case 'd':	lattice = atof(optarg);		// lattice constant (distance between two points in the discretized lattice)
		break;
      case 'M':	xmax = atof(optarg);		// size of the lattice in real units. symmetric around oasis, so lattice from -xmax to xmax
		break;
      case 'S':	itersteps = atoi(optarg);	// number of iterations to solve c and u*
		break;
      case 'g':	startwithgaussian = 1;		// special initial conditions
		break;
      case 'O':	outputstep = atoi(optarg);	// output of convergence criteria every "outputstep" steps (criteria usually intergral over u/c)
		break;
      case 'U': onlyU = 1;			// flag "-U", only the u* ODE is solved, not c
		break;
      case 'u': strcpy(ufilename,optarg);	// write the solution u(x) to a file (name as parameter to "-u")
		write_u_to_file =1;
		break;
      case 'o': u_outputstep = atoi(optarg);	// only used with option "-u FILENAME"
		break;				// output only every U_OUTPUTSTEP lattice point for the u-file
						// useful when the numerical solution is on a finer grid than the stochastic simulation
      case 'P': periodic_bc = 1;		// use periodic boundary conditions. default is reflecting boundaries.
		break;
    }
  }
}

double get_popsize(double r) {
  // convection-less case
  double rho = 4.*diff*r/(alpha*alpha);
  return diff/alpha/alpha * 16./(1+sqrt(rho)-2.*rho);
}

double get_u(double x) {
  // small velocity approximation
  double r_effective = r + v*v/4.;
  double aa = (alpha-2.*sqrt(diff*r_effective))/(alpha+2.*sqrt(diff*r_effective));
  double sqrt_ro2 = sqrt(r_effective/diff)/2.;
  double correctforpopsize = get_popsize(r)/get_popsize(r_effective);
  double tmp_numerator = exp(sqrt_ro2*fabs(x)) - aa * exp(-sqrt_ro2*fabs(x));
  return 3.*r_effective*aa*correctforpopsize/(tmp_numerator*tmp_numerator)*exp(0.5*v*x/diff);
}

double get_c(double x) {
  // small velocity approximation
  double r_effective = r + v*v/4.;
  double aa = (alpha-2.*sqrt(diff*r_effective))/(alpha+2.*sqrt(diff*r_effective));
  double sqrt_ro2 = sqrt(r_effective/diff)/2.;
  double correctforpopsize = get_popsize(r)/get_popsize(r_effective);
  double u2int = 3./32.*(alpha-2.*sqrt(diff*r_effective))*(alpha-2.*sqrt(diff*r_effective))*(alpha+4.*sqrt(diff*r_effective));
  double tmp_numerator = exp(sqrt_ro2*fabs(x)) - aa * exp(-sqrt_ro2*fabs(x));
  return 3.*r_effective*aa*correctforpopsize/(tmp_numerator*tmp_numerator)*exp(-0.5*v*x/diff)/u2int;
}

double get_gaussian(double x, double sigma2) {
  return exp(-x*x/sigma2);
}


void initialize() {
  int i,j;
  double x;
  double u2int_inv;
  double sigma2;
  double N_expectednovelocity,rho;

  double Dl2i = diff/(lattice*lattice);
  double v2li = 0.5*v/lattice;
  
  
  a_effective = alpha/lattice;
  
  space = 2*(int)(xmax/lattice)+1;
  space0 = (space-1)/2;
  
  r_sqrt = sqrt(r);
  
  rho = 4.*diff*r/(alpha*alpha);
  N_expectednovelocity = diff/alpha/alpha * 16./(1+sqrt(rho)-2.*rho);
  
  fprintf(stderr,"# parameters:\n");
  fprintf(stderr,"#    alpha      = %lf\n",alpha);
  fprintf(stderr,"#    r          = %lf\n",r);
  fprintf(stderr,"#    v          = %lf\n#\n",v);
  fprintf(stderr,"#    xmax       = %lf\n",xmax);
  fprintf(stderr,"#    lattice    = %lf\n",lattice);
  fprintf(stderr,"#    iterations = %d\n",itersteps);
  fprintf(stderr,"#    space      = %d\n",space);
  fprintf(stderr,"#    space0     = %d\n#\n",space0);
  fprintf(stderr,"#    N_expectv0 = %lf\n",N_expectednovelocity);
  
  c = (double*)malloc(space*sizeof(double));
  u = (double*)malloc(space*sizeof(double));
  cnew = (double*)malloc(space*sizeof(double));
  unew = (double*)malloc(space*sizeof(double));
  
  c_coeff_0  = (double*)malloc(space*sizeof(double));	// coefficients in the discretized ODE for c[i], c[i-1] and c[i+1]
  c_coeff_m1 = (double*)malloc(space*sizeof(double));
  c_coeff_p1 = (double*)malloc(space*sizeof(double));

  u_coeff_0  = (double*)malloc(space*sizeof(double));	// coefficients in the discretized ODE for u[i], u[i-1] and u[i+1]
  u_coeff_m1 = (double*)malloc(space*sizeof(double));
  u_coeff_p1 = (double*)malloc(space*sizeof(double));
  
  uder_coeff_0 = (double*)malloc(space*sizeof(double));	// ODE differentiated with respect to u, for newton iteration
  cder_const = (double*)malloc(space*sizeof(double));	//				      c
  
  sigma2 = r/(4.*diff);					// (characteristic lenght)^2, only used with gaussian initial conditions
  
  
  // set initial values for u,c and all the coefficients
  for(i=0;i<space;i++) {
    x = (i-space0)*lattice;
    if(startwithgaussian == 1) {
      u[i] = get_gaussian(x,sigma2);
      c[i] = get_gaussian(x,sigma2);
    }else{
      u[i] = get_u(x);
      c[i] = get_c(x);
    }
    u_coeff_m1[i] = Dl2i + v2li;
    u_coeff_0[i] = -2.*Dl2i - r;
    u_coeff_p1[i] = Dl2i - v2li;
    if(i==space0)u_coeff_0[i] += a_effective;
    
    c_coeff_m1[i] = Dl2i - v2li;
    c_coeff_0[i] = -2.*Dl2i - r;
    c_coeff_p1[i] = Dl2i + v2li;
    if(i==space0)c_coeff_0[i] += a_effective;
    
    uder_coeff_0[i] = -2.*Dl2i-r;
    if(i==space0)uder_coeff_0[i] += a_effective;
    
    cder_const[i] = uder_coeff_0[i];
    
    
    
    // periodic_bc == 0 means reflecting boundary condition,
    // no dependence on convection anymore, coefficient at boundaries is doubled or zero.
    if((periodic_bc == 0)&&(i==0)) {
      u_coeff_m1[i] = 0.;
      u_coeff_p1[i] = 2.*Dl2i;
      c_coeff_m1[i] = 0.;
      c_coeff_p1[i] = 2.*Dl2i;
    }
    if((periodic_bc == 0)&&(i==space-1)) {
      u_coeff_m1[i] = 2.*Dl2i;
      u_coeff_p1[i] = 0.;
      c_coeff_m1[i] = 2.*Dl2i;
      c_coeff_p1[i] = 0.;
    }
      
    
  }
}

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
  free(uder_coeff_0);
  free(cder_const);
}

int iterate_u(int timestep) {
  int i,j;
  double f,fu;
  double u_int = 0.0;

  // periodic boundary condition, hence u[0] and u[space-1] (first and last points) have to be treated separately (to use correct indices)
  // by using the predefined coefficients "u_coeff_?" and not doing calculations every step, the speedup is something between 15-20%
  f = u_coeff_m1[0] * u[space-1] + u_coeff_0[0]*u[0] + u_coeff_p1[0]*u[1] - 2.*u[0]*u[0];
  fu = uder_coeff_0[0]-4.*u[0];
  unew[0] = u[0] - f/fu;

  for(i=1;i<space-1;i++) {
    f = u_coeff_m1[i]*u[i-1] + u_coeff_0[i]*u[i] + u_coeff_p1[i]*u[i+1]-2.*u[i]*u[i];
    fu = uder_coeff_0[i]-4.*u[i];
    unew[i] = u[i] - f/fu;
  }
  
  f = u_coeff_m1[space-1] * u[space-2] + u_coeff_0[space-1]*u[space-1] + u_coeff_p1[space-1]*u[0] - 2.*u[space-1]*u[space-1];
  fu = uder_coeff_0[space-1]-4.*u[space-1];
  unew[space-1] = u[space-1] - f/fu;
  
  // after updating unew, everything copied to u
  memcpy(&u[0],&unew[0],space*sizeof(double));
  
  
  // output integral over u to see convergence of algorithm
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)u_int += u[i];
    u_int *= lattice;
    printf("%8d %.10e\n",timestep,u_int);
  }
  return 0;
}

int iterate_c(int timestep) {
  int i,j;
  double f,fc;
  double c_int = 0.0;

  
  // same as for u, first and last treated separately
  f = c_coeff_m1[0] * c[space-1] + c_coeff_0[0]*c[0] + c_coeff_p1[0]*c[1];
  cnew[0] = c[0] - f/cder_const[0];

  for(i=1;i<space-1;i++) {
    f = c_coeff_m1[i]*c[i-1] + c_coeff_0[i]*c[i] + c_coeff_p1[i]*c[i+1];
    cnew[i] = c[i] - f/cder_const[i];
  }
  f = c_coeff_m1[space-1] * c[space-2] + c_coeff_0[space-1]*c[space-1] + c_coeff_p1[space-1]*c[0];
  cnew[space-1] = c[space-1] - f/cder_const[space-1];
  
  memcpy(&c[0],&cnew[0],space*sizeof(double));

  // output population size to see convergence
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)c_int += c[i];
    c_int *= lattice;
    printf("%8d %.10e\n",timestep,c_int);
  }
  
  return 0;
}

void print_output() {
  int i,j;
  double x;
  if(onlyU == 1) {
    fprintf(stderr,"# x\t\t--\t\tu(x)\n");
  }else{
    fprintf(stderr,"# x\t\tc(x)\t\t\tu(x)\t\t\tc(x)*u(x)");
  }
  for(i=0;i<space;i++) {
    x = (i-space0)*lattice;
    if(onlyU == 1) {			// no solution for c present with option "-U", so second column is just "nan"
					// (second column is still written, so that u(x) is ALWAYS in third column)
      fprintf(stderr,"%lf\tnan\t%.10e\n",x,u[i]);
    }else{						// second column now c, fourth u*c
      fprintf(stderr,"%lf\t%.10e\t%.10e\t%.10e\n",x,c[i],u[i],c[i]*u[i]);
    }
  }
}




void write_u() {
  // write u(x) solution to an additional file
  // which can then be read into the stochastic simulation
  // to have an arbitrary constraint
  
  FILE* fp;
  int i;
  double x;
  
  fp = fopen(ufilename,"w");
  if(fp!=NULL) {
    for(i=0;i<space;i+=u_outputstep) {
      x = (i-space0)*lattice;
      fprintf(fp,"%.3lf %.10e\n",x,u[i]);
    }
    fclose(fp);
  }
}
  

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
  
  fprintf(stderr,"#    time       = %d\n",timestep);
  fprintf(stderr,"#    N          = %10.6lf\n",popsize);
}


int main(int argn, char *argv[]) {
  int i,j;
  parsecommandline(argn,argv);
  initialize();
  
  
  for(i=0;i<itersteps;i++) iterate_u(i);
  
  
  if(write_u_to_file == 1) {
    write_u();
  }
  
  // with flag "-U" (onlyU == 1) only calculate the solution u*, do not iterate for c
  if(onlyU != 1) {
    
    for(i=0;i<space;i++) {
      // update coefficients for c with solution from u
      c_coeff_0[i] -= 2.*u[i];
      cder_const[i] -= 2.*u[i];
    }
    printf("\n");
    
    for(i=0;i<itersteps;i++) {
      iterate_c(i);
    }
  }

  consistency_check(itersteps);
  
  
  print_output();

  
  
  cleanup();
  return 0;
}


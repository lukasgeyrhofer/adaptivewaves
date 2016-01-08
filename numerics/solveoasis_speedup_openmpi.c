#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

double alpha = 0.1;
double a_effective = 0.0;
double r = 0.002;
double v = 0.0;
double diff = 1.0;
double lattice = .1;
double xmax = 500.;
int periods = 1;
int itersteps = 10000;
int verbose = 0;
int outputstep = 1000;
int startwithgaussian = 0;

double r_sqrt;

int space, space0;
double *c,*u,*unew,*cnew;

double *c_coeff_m1, *c_coeff_0, *c_coeff_p1;
double *u_coeff_m1, *u_coeff_0, *u_coeff_p1;
double *uder_coeff_0,*cder_const;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn, argv,"a:r:v:D:x:M:S:P:gVO:")) != -1){
    switch(c) {
      case 'a':	alpha = atof(optarg);
		break;
      case 'r':	r = atof(optarg);
		break;
      case 'v':	v = atof(optarg);
		break;
      case 'D':	diff = atof(optarg);
		break;
      case 'x':	lattice = atof(optarg);
		break;
      case 'M':	xmax = atof(optarg);
		break;
      case 'S':	itersteps = atoi(optarg);
		break;
      case 'P':	periods = atoi(optarg);
		break;
      case 'g':	startwithgaussian = 1;
		break;
      case 'V': verbose = 1;
		break;
      case 'O':	outputstep = atoi(optarg);
		break;
    }
  }
}

double get_popsize(double r) {
  double rho = 4.*diff*r/(alpha*alpha);
  return diff/alpha/alpha * 16./(1+sqrt(rho)-2.*rho);
}

double get_u(double x) {
  double r_effective = r + v*v/4.;
  double aa = (alpha-2.*sqrt(diff*r_effective))/(alpha+2.*sqrt(diff*r_effective));
  double sqrt_ro2 = sqrt(r_effective/diff)/2.;
  double correctforpopsize = get_popsize(r)/get_popsize(r_effective);
  double tmp_numerator = exp(sqrt_ro2*fabs(x)) - aa * exp(-sqrt_ro2*fabs(x));
  return 3.*r_effective*aa*correctforpopsize/(tmp_numerator*tmp_numerator)*exp(0.5*v*x/diff);
}

double get_c(double x) {
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


void initialize(int rank, int numprocs) {
  int i,j;
  double x;
  double u2int_inv;
  double sigma2;
  double N_expectednovelocity,rho;

  double Dl2i = diff/(lattice*lattice);
  double v2li = 0.5*v/lattice;
  
  
  a_effective = alpha/lattice;
  
  space = 2*(int)(xmax/lattice);
  space0 = (space)/2;
  
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
  
  c_coeff_0  = (double*)malloc(space*sizeof(double));
  c_coeff_m1 = (double*)malloc(space*sizeof(double));
  c_coeff_p1 = (double*)malloc(space*sizeof(double));

  u_coeff_0  = (double*)malloc(space*sizeof(double));
  u_coeff_m1 = (double*)malloc(space*sizeof(double));
  u_coeff_p1 = (double*)malloc(space*sizeof(double));
  
  uder_coeff_0 = (double*)malloc(space*sizeof(double));
  cder_const = (double*)malloc(space*sizeof(double));
  
  sigma2 = r/(4.*diff);
  
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
    
  }
  //dummy
}

void cleanup() {
  free(c);
  free(u);
  free(unew);
  free(cnew);
}

int iterate_u(int timestep) {
  int i,j;
  double f,fu;
  double Dl2i = diff/(lattice*lattice);
  double v2li = 0.5*v/lattice;
  double t[3];
  double u_int = 0.0;

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
  
  memcpy(&u[0],&unew[0],space*sizeof(double));
  
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
  double Dl2i = diff/(lattice*lattice);
  double v2li = 0.5*v/lattice;

  f = c_coeff_m1[0] * c[space-1] + c_coeff_0[0]*c[0] + c_coeff_p1[0]*c[1];
  cnew[0] = c[0] - f/cder_const[0];

  for(i=1;i<space-1;i++) {
    f = c_coeff_m1[i]*c[i-1] + c_coeff_0[i]*c[i] + c_coeff_p1[i]*c[i+1];
    cnew[i] = c[i] - f/cder_const[i];
  }
  f = c_coeff_m1[space-1] * c[space-2] + c_coeff_0[space-1]*c[space-1] + c_coeff_p1[space-1]*c[0];
  cnew[space-1] = c[space-1] - f/cder_const[space-1];
  
  memcpy(&c[0],&cnew[0],space*sizeof(double));
  
  return 0;
}

void print_output() {
  int i,j;
  double x;
  fprintf(stderr,"# x\t\tc(x)\t\t\tu(x)\t\t\tc(x)*u(x)\n");
  for(j=0;j<periods;j++) {
    for(i=0;i<space;i++) {
      x = (i-space0)*lattice+ 2*j*xmax;
      fprintf(stderr,"%lf\t%.10e\t%.10e\t%.10e\n",x,c[i],u[i],c[i]*u[i]);
    }
  }
}

void consistency_check(int timestep) {
  int i,j;
  double popsize = 0.;
  double sum = 0.;
  
  for(i=0;i<space;i++) {
    popsize += c[i];
    sum += u[i]*c[i];
  }
  popsize *= lattice;
  sum *= lattice;
  
  fprintf(stderr,"# results");
  if(verbose) {
    if(timestep>0) {
      fprintf(stderr," (at timestep %d)",timestep);
    }else{
      fprintf(stderr," (at end)");
    }
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"#    N          = %10.6lf\n",popsize);
  fprintf(stderr,"#    <u|c>      = %.10e\n",sum);
  fprintf(stderr,"#    1 - <u|c>  = %.10e\n",1.-sum);
}


int main(int argn, char *argv[]) {
  int i,j;
  int numprocs,rank;

  parsecommandline(argn,argv);
  
  
  MPI_Init(&argn, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  initialize(rank,numprocs);

  for(i=0;i<itersteps;i++) iterate_u(i);
  
  for(i=0;i<space;i++) {
    c_coeff_0[i] -= 2.*u[i];
    cder_const[i] -= 2.*u[i];
  }
  
  for(i=0;i<itersteps;i++) {
    iterate_c(i);
//     if(verbose)if(i%10000==0)consistency_check(i);
  }
  
  
  print_output();

  consistency_check(-1);
  
  
  cleanup();
  MPI_Finalize();
  return 0;
}


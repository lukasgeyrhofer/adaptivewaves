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
double lattice = 1.0;
double xmax = 500.;
double tstep = 1.0;
int periods = 1;
int itersteps = 10000;
int verbose = 0;
int outputstep = 1000;
int startwithgaussian = 0;

double r_sqrt;

int space, space0;
double *c,*u,*unew,*cnew;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn, argv,"a:r:v:D:x:M:S:T:P:gVO:")) != -1){
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
      case 'T': tstep = atof(optarg);
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


double get_u(double x) {
  double tmp_a,tmp_numerator;
  tmp_a = (alpha - 2.*r_sqrt)/(alpha + 2.*r_sqrt);
  tmp_numerator = exp(r_sqrt*fabs(x)*0.5) - tmp_a * exp(-r_sqrt*fabs(x)*0.5);
  return 3.*r*tmp_a/(tmp_numerator*tmp_numerator);
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
  
  a_effective = alpha/lattice;
  
  space = 2*(int)(xmax/lattice);
  space0 = (space)/2;
  
  r_sqrt = sqrt(r);
  u2int_inv = 32./(3.*(alpha*alpha*alpha-12.*alpha*r+16.*r*r_sqrt));
  
  rho = 4.*diff*r/(alpha*alpha);
  N_expectednovelocity = diff/alpha/alpha * 16./(1+sqrt(rho)-2.*rho);
  
  fprintf(stderr,"# parameters:\n");
  fprintf(stderr,"#    alpha      = %lf\n",alpha);
  fprintf(stderr,"#    r          = %lf\n",r);
  fprintf(stderr,"#    v          = %lf\n#\n",v);
  fprintf(stderr,"#    xmax       = %lf\n",xmax);
  fprintf(stderr,"#    lattice    = %lf\n",lattice);
  fprintf(stderr,"#    tstep      = %lf\n",tstep);
  fprintf(stderr,"#    iterations = %d\n",itersteps);
  fprintf(stderr,"#    space      = %d\n",space);
  fprintf(stderr,"#    space0     = %d\n#\n",space0);
  fprintf(stderr,"#    N_expectv0 = %lf\n",N_expectednovelocity);
  
  c = (double*)malloc(space*sizeof(double));
  u = (double*)malloc(space*sizeof(double));
  cnew = (double*)malloc(space*sizeof(double));
  unew = (double*)malloc(space*sizeof(double));
  
  sigma2 = r/(4.*diff);
  
  for(i=0;i<space;i++) {
    x = (i-space0)*lattice;
    if(startwithgaussian == 1) {
      u[i] = get_gaussian(x,sigma2);
      c[i] = get_gaussian(x,sigma2);
    }else{
      u[i] = get_u(x);
      c[i] = get_u(x)*u2int_inv;
    }
    
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

  f = Dl2i*(u[space-1]-2.*u[0]+u[1])-v2li*(u[1]-u[space-1])-r*u[0]-2.*u[0]*u[0];
  fu = -2.*Dl2i-r-4.*u[0];
  unew[0] = u[0] - tstep * f/fu;

  for(i=1;i<space-1;i++) {
    f = Dl2i*(u[i-1]-2.*u[i]+u[i+1])-v2li*(u[i+1]-u[i-1])-r*u[i]-2.*u[i]*u[i];
    fu = -2.*Dl2i-r-4.*u[i];
    if(i==space0) {
      f += a_effective * u[i];
      fu += a_effective;
    }
    unew[i] = u[i] - tstep * f/fu;
  }
  f = Dl2i*(u[space-2]-2.*u[space-1]+u[0])-v2li*(u[0]-u[space-2])-r*u[space-1]-2.*u[space-1]*u[space-1];
  fu = -2.*Dl2i-r-4.*u[space-1];
  unew[space-1] = u[space-1] - tstep * f/fu;
  
  memcpy(&u[0],&unew[0],space*sizeof(double));
  
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)u_int += u[i];
    u_int *= lattice;
    printf("%6d %.10e\n",timestep,u_int);
  }
  return 0;
}

int iterate_c(int timestep) {
  int i,j;
  double f,fc;
  double Dl2i = diff/(lattice*lattice);
  double v2li = 0.5*v/lattice;

  f = Dl2i*(c[space-1]-2.*c[0]+c[1])+v2li*(c[1]-c[space-1])-r*c[0]-2.*u[0]*c[0];
  fc = -2.*Dl2i-r-2.*u[0];
  cnew[0] = c[0] - tstep * f/fc;

  for(i=1;i<space-1;i++) {
    f = Dl2i*(c[i-1]-2.*c[i]+c[i+1])+v2li*(c[i+1]-c[i-1])-r*c[i]-2.*u[i]*c[i];
    fc = -2.*Dl2i-r-2.*u[i];
    if(i==space0) {
      f += a_effective * c[i];
      fc += a_effective;
    }
    cnew[i] = c[i] - tstep * f/fc;
  }
  f = Dl2i*(c[space-2]-2.*c[space-1]+c[0])+v2li*(c[0]-c[space-2])-r*c[space-1]-2.*u[space-1]*c[space-1];
  fc = -2.*Dl2i-r-2.*u[space-1];
  cnew[space-1] = c[space-1] - tstep * f/fc;
  
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
  parsecommandline(argn,argv);
  initialize();
  
  
  for(i=0;i<itersteps;i++) iterate_u(i);
  
  for(i=0;i<itersteps;i++) {
    iterate_c(i);
//     if(verbose)if(i%10000==0)consistency_check(i);
  }
  
  
  print_output();

  consistency_check(-1);
  
  
  cleanup();
  return 0;
}


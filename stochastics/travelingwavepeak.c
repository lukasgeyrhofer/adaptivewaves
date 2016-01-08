#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>

int space = 100;
int space0 = 50;
double dx = 1e-2;
int maxSteps = 1000;

int outputstep = 100;
int quiet = 0;

double epsilon = 1e-2, twoepssqrt;
double mutationrate = 1e-5;

double populationsize = 1e0;

int shift_threshold = 1;
int xmean_lastshift_bin;

gsl_vector *nn;
gsl_vector *binfitness;
gsl_vector *tmp;
gsl_vector *ones;

const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
unsigned long int randseed = 0;

void parsecomamndline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn,argv,"s:z:d:S:e:D:O:qQT:R:N:")) != -1) {
    switch(c) {
      case 's':	space = atoi(optarg);
		break;
      case 'z':	space0 = atoi(optarg);
		break;
      case 'd':	dx = atof(optarg);
		break;
      case 'S':	maxSteps = atoi(optarg);
		break;
      case 'e': epsilon = atof(optarg);
		break;
      case 'D':	mutationrate = atof(optarg);
		break;
      case 'O':	outputstep = atoi(optarg);
		break;
      case 'q':	quiet = 2;
		break;
      case 'Q':	quiet = 1;
		break;
      case 'T':	shift_threshold = atoi(optarg);
		break;
      case 'R':	randseed = atoi(optarg);
		break;
      case 'N':	populationsize = atof(optarg);
		break;
    }
  }
  if(randseed==0)randseed = time(NULL);
}



  
double get_mean_fit() {
  double retval;
  gsl_blas_ddot(nn,binfitness,&retval);
  return retval;
}


double get_mean_fit2() {
  double retval;
  gsl_vector_memcpy(tmp,nn);
  gsl_vector_mul(tmp,binfitness);
  gsl_blas_ddot(tmp,binfitness,&retval);
  return retval;
}




void shift_forward(int step) {
  int i;
  for(i=0;i<space-step;i++) {
    gsl_vector_set(nn,i,gsl_vector_get(nn,i+step));
  }
  for(i=space-step;i<space;i++) {
    gsl_vector_set(nn,i,0.);
  }
}

void shift_backward(int step) {
  int i;
  for(i=space-1;i>step;i--) {
    gsl_vector_set(nn,i,gsl_vector_get(nn,i-step));
  }
  for(i=step-1;i>=0;i--) {
    gsl_vector_set(nn,i,0.);
  }
}


double shift_population(double xmean_diff,int timestep) {
  int shift = (int)floor(xmean_diff/dx);
  xmean_lastshift_bin += shift;
//   gsl_vector_fprintf(stdout,nn,"## %14.10e");
//   fprintf(stdout,"## --\n");
  if(shift > 0) {
    fprintf(stdout,"#shift forward: %d [time %d]\n",shift,timestep);
    shift_forward(shift);
    return -shift * dx;
  }
  if(shift < 0) {
    fprintf(stdout,"#shift backward: %d [time %d]\n",shift,timestep);
    shift_backward(shift);
    return shift * dx;
  }

//   gsl_vector_fprintf(stdout,nn,"## %14.10e");
//   exit(0);
  
}

  
double reproduce(int timestep) {
  double xmean = get_mean_fit();
  double nntmp;
  int i;
  
  if((xmean > shift_threshold*dx) || (xmean < -shift_threshold*dx))xmean += shift_population(xmean,timestep);
  
  nntmp = gsl_vector_get(nn,0)*(1.+epsilon*(gsl_vector_get(binfitness,0)-xmean-mutationrate));
//   nntmp += twoepssqrt*(gsl_ran_poisson(rg,nntmp)-nntmp);
  gsl_vector_set(tmp,0,nntmp);
  for(i=1;i<space;i++) {
    nntmp = gsl_vector_get(nn,i)*(1.+epsilon*(gsl_vector_get(binfitness,i)-xmean-mutationrate))+gsl_vector_get(nn,i-1)*epsilon*mutationrate;
//     nntmp += twoepssqrt*(gsl_ran_poisson(rg,nntmp)-nntmp);
    gsl_vector_set(tmp,i,nntmp);
  }
  gsl_vector_memcpy(nn,tmp);
  
  return xmean;
}


void populationconstraint() {
  double currentpopsize = 0.;
  int i;
  for(i=0;i<space;i++) {
    currentpopsize += gsl_vector_get(nn,i);
  }
  gsl_vector_scale(nn,populationsize/currentpopsize);
}

 
void print_populationdensity(int timestep,double meanfitness) {
  int i;
  for(i=0;i<space;i++) {
    fprintf(stderr,"%d %14.10lf %20.10e\n",timestep,gsl_vector_get(binfitness,i)+xmean_lastshift_bin*dx,gsl_vector_get(nn,i));
  }
  fprintf(stderr,"\n");
}



int initialize() {
  int i;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);

  
  twoepssqrt = 2.*sqrt(epsilon);
  
  nn = gsl_vector_calloc(space);
  gsl_vector_set(nn,space0,populationsize);
  
  binfitness = gsl_vector_calloc(space);
  for(i=0;i<space;i++) {
    gsl_vector_set(binfitness,i,(i-space0)*dx);
  }
  
  tmp = gsl_vector_calloc(space);
  ones = gsl_vector_alloc(space);
  for(i=0;i<space;i++)gsl_vector_set(ones,i,1.);
  
  xmean_lastshift_bin = 0;
  
}
  


  
void cleanup() {
  gsl_vector_free(nn);
  gsl_vector_free(binfitness);
  gsl_vector_free(tmp);
  gsl_vector_free(ones);
}


int main(int argn, char *argv[]) {
  int i;
  double x,xx;
  parsecomamndline(argn,argv);
  
  initialize();

  for(i=0;i<maxSteps;i++) {
    x = reproduce(i);
    populationconstraint();
    if(i%outputstep == 0) {
      if(quiet<2) {
	xx = get_mean_fit2();
	fprintf(stdout,"%5d %14.10lf %14.10lf %14.10lf\n",i,x,xx-x*x,xmean_lastshift_bin*dx);
      }
      if(quiet==0)print_populationdensity(i,x);
    }
  }

  cleanup();

  return 0;
}


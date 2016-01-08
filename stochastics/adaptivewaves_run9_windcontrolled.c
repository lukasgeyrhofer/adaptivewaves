// Stochastic simulation of adaptation in the high-mutation, strong-selection regime
// based on travelling waves
// Lukas Geyrhofer, 2012-2013

// version 0.1 (first named version)
//	130130		included different mutation models

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
  

// commandline options

char cinfilename[128];
char uinfilename[128];
char coutfilename[128];
int write_conf_to_file = 0;

double epsilon = 0.05;
int maxSteps = 10000;
int outputstep = 100;

int quiet = 0;


// from c file
double dx;
int space;
int space0;

double rho;
double starttime;
double driftvelocity;


// densities, profiles, coefficients

double *n, *n_new;		// dynamics based on occupancies n[i],
double *c;			// in- and output is still based on densities c[i], (n[i] = c[i]*dx)
double *u;
double *coeff_m,*coeff_0,*coeff_p;
double *zeros;

// time averages

int timeavg = 0;
char timeavg_filename[128];
double **timeavg_moments;
double *timeavg_ones;
double *timeavg_tmp;
double timeavg_count= 0.;
int timeavg_countmoments = 2;


// covariance

int covar = 0;
char covar_filename[128];
double **covariance;
double covar_count = 0.;
int covar_outputstep=10;
int covar_offset = 0;
int covar_space;


// integrated moments of density

double dens_mom0,dens_mom1,dens_mom2;
double dv = 0.;
double uc = 0.;

// random number generator variables

const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
unsigned long int randseed = 0;


// other constants and variables

double lambda = 1.0;
double twoepssqrt;
double twoeps;

// mutation models

int mutationtype = 1;	// 1 - diffusion
			//      -m diffusion
			// 2 - single (beneficial) jumps
			//      -m jumps
			// 3 - exponential kernel, only beneficial
			// 4 - exponential kernel, beneficial AND deleterious
			// 5 - kernel from external file
double diffusionconstant;



int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char* argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn,argv,"c:C:u:e:S:O:R:t:M:qQv:V:m:")) != -1) {
    switch(c) {
      case 'c':	strcpy(cinfilename,optarg);
		haveinfile++;
		break;
      case 'C': strcpy(coutfilename,optarg);
		write_conf_to_file = 1;
		break;
      case 'u': strcpy(uinfilename,optarg);
		haveinfile++;
		break;
      case 'e':	epsilon=atof(optarg);
		break;
      case 'S': maxSteps = atoi(optarg);
		break;
      case 'O':	outputstep = atoi(optarg);
		break;
      case 'R':	randseed = atoi(optarg);
		break;
      case 't':	strcpy(timeavg_filename,optarg);
		timeavg = 1;
		break;
      case 'M':	timeavg_countmoments = atoi(optarg);
		break;
      case 'Q': quiet = 1;
		break;
      case 'q': quiet = 2;
		break;
      case 'v':	strcpy(covar_filename,optarg);
		covar = 1;
		break;
      case 'V':	covar_outputstep = atoi(optarg);
		break;
    }
  }
  if(haveinfile<2)print_error("need both options -c and -u");
  if(randseed==0)randseed=time(NULL);
}




void read_c() {
  int i;
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fp;
  
  // file format of binary files:
  // 2 integers, specifying the number of integer and double parameters, respectively
  // 1 double, 2 integers, specifying the lattice: dx, space (number of lattice point), space0 (index of x=0)
  // all integer parameters
  // all double parameters: for adaptive waves: rho, starttime, 
  
  
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
  
  fread(&rho,1,sizeof(double),fp);
  fread(&starttime,1,sizeof(double),fp);
  fread(&diffusionconstant,1,sizeof(double),fp);
  fread(&driftvelocity,1,sizeof(double),fp);
  
  if(dcount > 4) {
    tmpd = (double*)malloc((dcount-4)*sizeof(double));
    fread(&tmpd[0],dcount-4,sizeof(double),fp);
    free(tmpd);
  }
  
  c = (double*)malloc(space*sizeof(double));
  n = (double*)malloc(space*sizeof(double));
  fread(&c[0],space,sizeof(double),fp);
  for(i=0;i<space;i++)n[i] = c[i]*dx;
  fclose(fp);
}

void write_c() {
  int i;
  FILE *fp;
  int icount = 0,dcount = 4;
  double endtime = starttime + epsilon*maxSteps;
  
  fp = fopen(coutfilename,"wb");
  if(fp==NULL)print_error("count not open file for writing");
  
  fwrite(&icount,1,sizeof(int),fp);
  fwrite(&dcount,1,sizeof(int),fp);
  
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  fwrite(&rho,1,sizeof(double),fp);
  fwrite(&endtime,1,sizeof(double),fp);
  fwrite(&diffusionconstant,1,sizeof(double),fp);
  fwrite(&driftvelocity,1,sizeof(double),fp);
  
  
  for(i=0;i<space;i++)c[i]=n[i]/dx;
  fwrite(&c[0],space,sizeof(double),fp);
  
  fclose(fp);
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
  
  if ((u_space!=space)||(u_space0 != space0))print_error("lattice of u and c not identical!");
  
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
 
 
void calc_timeavg() {
  int i,j;
  memcpy(&timeavg_tmp[0],&timeavg_ones[0],space*sizeof(double));
  for(j=0;j<timeavg_countmoments;j++) {
    for(i=0;i<space;i++) {
      timeavg_tmp[i] *= n[i]/dx;
      timeavg_moments[j][i] += timeavg_tmp[i];
    }
  }
  timeavg_count +=1.;
}


void write_timeavg() {
  FILE *fp;
  int i,j;

  fprintf(stdout,"# writing time-averaged moments to '%s'\n",timeavg_filename);
  
  fp = fopen(timeavg_filename,"w");
  fprintf(fp,"# simulationtime = %lf\n",maxSteps*epsilon);
  fprintf(fp,"# datapoints = %.0lf\n",timeavg_count);
  for(i=0;i<space;i++) {
    fprintf(fp,"%lf",(i-space0)*dx);
    for(j=0;j<timeavg_countmoments;j++) {
      fprintf(fp,"\t%e",timeavg_moments[j][i]/timeavg_count);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
 
 
 
void calc_covar() {
  int i,j;
  int ii,jj;
  for(i=0;i<covar_space;i++) {
    ii = i*covar_outputstep+covar_offset;
    for(j=0;j<=i;j++) {
      jj = j*covar_outputstep+covar_offset;
      covariance[i][j] += n[ii]*n[jj];
    }
  }
  covar_count += 1.;
}


void write_covar() {
  int i,j;
  FILE *fp;
  int covar_index;
  
  fprintf(stdout,"# writing covariance matrix to '%s'\n",covar_filename);
  
  fp = fopen(covar_filename,"w");
  if(fp!=NULL) {
    for(i=0;i<covar_space;i++) {
      for(j=0;j<covar_space;j++) {
	if(j<=i) {
	  fprintf(fp,"%lf %lf %20.10e\n",(i*covar_outputstep+covar_offset-space0)*dx,(j*covar_outputstep+covar_offset-space0)*dx, \
	  covariance[i][j]/(dx*dx*covar_count));
	}else{
	  fprintf(fp,"%lf %lf %20.10e\n",(i*covar_outputstep+covar_offset-space0)*dx,(j*covar_outputstep+covar_offset-space0)*dx,covariance[j][i]/(dx*dx*covar_count));
	}
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
}
 
void reproduction() {
  int i;
  int last = space-1;
  
  n_new[0] = 0.;
  for(i=1;i<space-1;i++) {
    n_new[i] =  coeff_m[i]*n[i-1] + coeff_0[i]*n[i] + coeff_p[i]*n[i+1];
  }
  memcpy(n,n_new,space*sizeof(double));
  for(i=1;i<space-1;i++)if(n[i]<1e9)n[i] += twoepssqrt*(gsl_ran_poisson(rg,n[i])-n[i]);

//   n_new[last] = 0.;
  
//   memcpy(&n[0],&n_new[0],space*sizeof(double));
}


double get_int_uc() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++)sum+=u[i]*n[i];
  return sum;
}

double get_int_ucp() {
  int i;
  double sum = 0.;
  for(i=1;i<space-1;i++)sum += u[i]*(n[i+1]-n[i-1]);
  return sum/(2.*dx*dx);
}

double get_dv() {
  int i;
  double sum1 = 0., sum2 = 0.;
  uc = 0.;
  for(i=1;i<space-1;i++) {
    uc += u[i]*n[i];
    sum2 += u[i]*(n[i+1]-n[i-1]);
  }
  return (1.-1./uc)/sum2; // *2dx    // this cancels out anyway below!
}

void update_integrated_moments() {
  int i;
  double x;
  dens_mom0 = 0.;
  dens_mom1 = 0.;
  dens_mom2 = 0.;
  for(i=0;i<space;i++) {
    x = (i-space0)*dx;
    dens_mom0 += n[i];
    dens_mom1 += x*n[i];
    dens_mom2 += x*x*n[i];
  }
}

void populationcontrol() {
  int i;
//   double dv;
  memcpy(n_new,n,space*sizeof(double));
  dv = get_dv();
  for(i=1;i<space-1;i++) {
    n_new[i] -= dv*(n[i+1]-n[i-1]);  // actually dv/2dx, but 2dx is canceled by the computation above
    if(n_new[i]<0)n_new[i]=0;
  }
  memcpy(n,n_new,space*sizeof(double));
}


void write_output(int step) {
  int i;
  double a;
  if(quiet < 2) {
    a=dens_mom1/dens_mom0;
    fprintf(stdout,"%10.3lf %17.10e %17.10e %17.10e %17.10e %17.10e\n",epsilon*step,a,dens_mom2/dens_mom0-a*a,dens_mom0,dv,uc);
    if(quiet < 1) {
      for(i=0;i<space;i++) {
	fprintf(stderr,"%lf %lf %g\n",epsilon*step,(i-space0)*dx,n[i]/dx);
      }
      fprintf(stderr,"\n");
    }
  }
}



void print_options() {
  if(quiet < 2) {
    fprintf(stdout,"# stochastic simulation of adaptive waves...\n");
    fprintf(stdout,"# options provided by commandline:\n");
    fprintf(stdout,"#     epsilon        = %g\n",epsilon);
    fprintf(stdout,"#     maxSteps       = %d\n",maxSteps);
    fprintf(stdout,"#     outputstep     = %d\n",outputstep);
    fprintf(stdout,"#     cfile          = %s\n",cinfilename);
    fprintf(stdout,"#     ufile          = %s\n",uinfilename);
    fprintf(stdout,"#     RANDSEED       = %d\n",randseed);
    
    fprintf(stdout,"# options read from density-file:\n");
    fprintf(stdout,"#     dx             = %g\n",dx);
    fprintf(stdout,"#     space          = %d\n",space);
    fprintf(stdout,"#     space0         = %d\n",space0);
    fprintf(stdout,"#     rho            = %g\n",rho);
    fprintf(stdout,"#     starttime      = %g\n",starttime);
    fprintf(stdout,"#     diffusionconst = %g\n",diffusionconstant);
    fprintf(stdout,"#     driftvelocity  = %g\n",driftvelocity);
  }  
}  



void initialize() {
  int i,j;
  double x;
  
  read_c();
  read_u();
  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  n_new   = (double*)malloc(space*sizeof(double));
  
  coeff_m = (double*)malloc(space*sizeof(double));
  coeff_0 = (double*)malloc(space*sizeof(double));
  coeff_p = (double*)malloc(space*sizeof(double));
  
  zeros = (double*)malloc(space*sizeof(double));
  
  
  coeff_m[0] = 0.;
  coeff_0[0] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (0-space0)*dx*epsilon - rho*epsilon;
  coeff_p[0] = epsilon*diffusionconstant/(dx*dx) + driftvelocity*epsilon/dx;
  zeros[0]   = 0.;
  
  for(i=1;i<space-1;i++) {
    coeff_m[i] = epsilon*diffusionconstant/(dx*dx) - 0.5*driftvelocity*epsilon/dx;
    coeff_0[i] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (i-space0)*dx*epsilon - rho*epsilon;
    coeff_p[i] = epsilon*diffusionconstant/(dx*dx) + 0.5*driftvelocity*epsilon/dx;
    zeros[i]   = 0.;
  }
  coeff_m[space-1] = epsilon*diffusionconstant/(dx*dx) - driftvelocity*epsilon/dx;
  coeff_0[space-1] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (space-1-space0)*dx*epsilon - rho*epsilon;
  coeff_p[space-1] = 0.;
  zeros[space-1]   = 0.;
  
  
  twoepssqrt = sqrt(2.*epsilon);
  twoeps = 2.*epsilon;
  
  if(timeavg == 1) {
    timeavg_moments = (double**)malloc(timeavg_countmoments*sizeof(double*));
    timeavg_ones = (double*)malloc(space*sizeof(double));
    timeavg_tmp  = (double*)malloc(space*sizeof(double));
    for(j=0;j<timeavg_countmoments;j++) {
      timeavg_moments[j] = (double*)malloc(space*sizeof(double));
      memcpy(&timeavg_moments[j][0],&zeros[0],space*sizeof(double));
    }
    for(i=0;i<space;i++)timeavg_ones[i] = 1.;
  }
  
  
  if(covar == 1) {
    covar_space = space/covar_outputstep;
    covar_offset = space0%covar_outputstep;
    if(covar_offset==0)covar_space++;
    covariance = (double**)malloc(covar_space*sizeof(double*));
    for(j=0;j<covar_space;j++) {
      covariance[j] = (double*)malloc((j+1)*sizeof(double));
      memcpy(&covariance[j][0],&zeros[0],(j+1)*sizeof(double));
    }
  }
  
}



void cleanup() {
  int i;
  free(c);
  free(n);
  free(n_new);
  free(u);
  free(coeff_m);
  free(coeff_0);
  free(coeff_p);  
  free(zeros);
  if(timeavg==1) {
    for(i=0;i<timeavg_countmoments;i++)free(timeavg_moments[i]);
    free(timeavg_moments);
    free(timeavg_ones);
    free(timeavg_tmp);
  }
  if(covar==1) {
    for(i=0;i<covar_space;i++)free(covariance[i]);
    free(covariance);
  }
}


int main(int argn, char* argv[]) {
  int i;
  
  parsecommandline(argn,argv);
  
  initialize();
  print_options();
  
  update_integrated_moments();
  write_output(0);
  
  for(i=1;i<=maxSteps;i++) {
    reproduction();
    populationcontrol();
    if(i%outputstep==0) {
      update_integrated_moments();
      write_output(i);
      if(timeavg==1)calc_timeavg();
      if(covar==1)calc_covar();
    }
  }
  
//   printf("1\n");
  
  if(timeavg==1)write_timeavg();
//   printf("2\n");
  if(covar==1)write_covar();
//   printf("3\n");
  if(write_conf_to_file == 1)write_c();
//   printf("4\n");
  
  cleanup();
  
  
  return 0;
}
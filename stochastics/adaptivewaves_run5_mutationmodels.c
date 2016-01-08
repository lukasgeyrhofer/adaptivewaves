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


// definitions for function pointers

typedef void (*repr_func)();
typedef void (*popctrl_func)();
// typedef void (*avg_func)();

repr_func reproduction;
popctrl_func populationcontrol;
// avg_func calc_covar;
// avg_func calc_timeavg;
// avg_func update_integrated_moments;


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
double *coeff_m,*coeff_0,*coeff_p, *coeff_j,*coeff_kernel;
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
			//      -m exp
			// 4 - kernel from external file
			//      -m file
double diffusionconstant;
int allowmutationoption = 0;
double mutationrate = 0.001;
double mutation_exprate = 0.001;
int mutation_jumpwidth = 1;
char mutation_kernelfile[128];
double *mutationkernel;
int space_m,space0_m;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char* argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn,argv,"c:C:u:e:S:O:R:t:T:qQv:V:m:M:U:")) != -1) {
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
      case 'T':	timeavg_countmoments = atoi(optarg);
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
      case 'm':	if(strcmp(optarg,"diffusion") == 0) {
		  mutationtype = 1;
		}else if(strcmp(optarg,"jumps") == 0) {
		  mutationtype = 2;
		}else if(strcmp(optarg,"exp") == 0) {
		  mutationtype = 3;
		}else if(strcmp(optarg,"file") == 0) {
		  mutationtype = 4;
		}else{
		  print_error("could not find mutationtype! possible options: 'diffusion', 'jumps', 'exp', 'file'");
		}
		allowmutationoption = 1;
		break;
      case 'M': if(allowmutationoption == 1) {
		  switch(mutationtype) {
		    case 1:	diffusionconstant = atof(optarg);
				break;
		    case 2:	mutation_jumpwidth = atoi(optarg);
				break;
		    case 3:	mutation_exprate = atof(optarg);
				break;
		    case 4:	strcpy(mutation_kernelfile,optarg);
				break;
		  }
		}else{
		  print_error("mutationtype not specified!");
		}
		break;
      case 'U':	if(allowmutationoption == 1) {
		  mutationrate = atof(optarg);
		}else{
		  print_error("mutationtype not specified!");
		}
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



void read_mutationkernel() {
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  FILE *fp;
  double m_dx;
  
  fp = fopen(mutation_kernelfile,"rb");
  fread(&icount,1,sizeof(int),fp);
  fread(&dcount,1,sizeof(int),fp);
  
  fread(&m_dx,1,sizeof(double),fp);
  fread(&space_m,1,sizeof(int),fp);
  fread(&space0_m,1,sizeof(int),fp);
  
  if((space_m != 2*space+1) && (space0_m != space))print_error("mutation kernel has to be twice as large as population lattice");
  
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
  
  mutationkernel = (double*)malloc(space_m*sizeof(double));
  fread(&mutationkernel[0],space_m,sizeof(double),fp);
  
  fclose(fp);
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


void write_output(int step) {
  int i;
  if(quiet < 2) {
    fprintf(stdout,"%10.4lf %16.3lf %16.10lf %16.10lf %16.10lf\n",epsilon*step,dens_mom0,dens_mom1/dens_mom0,dens_mom2/dens_mom0,lambda);
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
    fprintf(stdout,"# stochasitic simulation of adaptive waves...\n");
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




void generate_exp_mutationkernel() {
  int i;
  space_m = 2*space+1;
  space0_m = space;
  mutationkernel = (double*)malloc(space_m*sizeof(double));
  for(i=0;i<space0_m;i++)mutationkernel[i] = 0.;
  for(i=space0_m;i<space_m;i++)mutationkernel[i] = 1./mutation_exprate*exp((i-space0_m)*dx/mutation_exprate);
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

 

double get_int_uc() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++)sum+=u[i]*n[i];
  return sum;
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



// =====================================================================================================
// ==  POPULATIONCONTROL ==
// =====================================================================================================



void populationcontrol_tunedmodel() {
  int i;
  double lambda_inv;
  lambda = get_int_uc();
  lambda_inv = 1./lambda;
  for(i=0;i<space;i++)n[i] *= lambda_inv;
}

void populationcontrol_fixedNmodel() {
  int i;
  double popsize_cur = 0., ninv;
  for(i=0;i<space;i++)popsize_cur += n[i];
  ninv = 1./popsize_cur;
  for(i=0;i<space;i++)n[i] *= ninv;
}


// =====================================================================================================
// ==  REPRODUCTION ==
// =====================================================================================================

void reproduction_diffusion() {
  int i;
  int last = space-1;
  
  n_new[0] = 0.;
  for(i=1;i<space-1;i++) {
    n_new[i] =  coeff_m[i]*n[i-1] + coeff_0[i]*n[i] + coeff_p[i] * n[i+1];
    n_new[i] += twoepssqrt*(gsl_ran_poisson(rg,n[i])-n[i]);
  }
  n_new[last] = 0.;
  
  memcpy(&n[0],&n_new[0],space*sizeof(double));
}

void reproduction_jumps() {
  int i;
  int last = space-1;

  n_new[0] = 0.;
  for(i=1;i<mutation_jumpwidth;i++) {
    n_new[i] = coeff_m[i]*n[i-1] + coeff_0[i]*n[i] + coeff_p[i] * n[i+1];
    n_new[i] += twoepssqrt*(gsl_ran_poisson(rg,n[i])-n[i]);
  }
  for(i=mutation_jumpwidth;i<space-1;i++) {
    n_new[i] = coeff_m[i]*n[i-1] + coeff_0[i]*n[i] + coeff_p[i] * n[i+1] + coeff_j[i] * n[i-mutation_jumpwidth];
    n_new[i] += twoepssqrt*(gsl_ran_poisson(rg,n[i])-n[i]);
  }
  n_new[last] = 0.;
  
  memcpy(&n[0],&n_new[0],space*sizeof(double));
}


void reproduction_expkernel() {
  int i,j;
  int last = space-1;
  
  n_new[0] = 0.;
  for(i=1;i<space;i++) {
    n_new[i] = coeff_m[i] * n[i-1] + coeff_0[i] * n[i] + coeff_p[i] * n[i+1];
    for(j=0;j<i;j++) n_new[i] += coeff_kernel[j]*n[i-j];
    n_new[i] += twoepssqrt*(gsl_ran_poisson(rg,n[i])-n[i]);
  }
  n_new[last] = 0.;
  
  memcpy(&n[0],&n_new[0],space*sizeof(double));
}

void reproduction_kernel() {
  int i,j;
  int last = space-1;
  
  n_new[0] = 0.;
  for(i=1;i<space-1;i++) {
    n_new[i] =  coeff_m[i]*n[i-1] + coeff_0[i]*n[i] + coeff_p[i] * n[i+1];
    for(j=0;j<space;j++)n_new[j] += coeff_kernel[j-i+space0_m]*n[j];
    n_new[i] += twoepssqrt*(gsl_ran_poisson(rg,n[i])-n[i]);
  }
  n_new[last] = 0.;
  
  memcpy(&n[0],&n_new[0],space*sizeof(double));
}




void initialize() {
  int i,j;
  double x;
  double integrated_mutationkernel;
  
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
  coeff_j = (double*)malloc(space*sizeof(double));
  
  zeros = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++)zeros[i] = 0.;
  
  
  populationcontrol = populationcontrol_tunedmodel;
  
  switch(mutationtype) {
    case 1:	reproduction = reproduction_diffusion;
		coeff_m[0] = 0.;
		coeff_0[0] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (0-space0)*dx*epsilon - rho*epsilon;
		coeff_p[0] = epsilon*diffusionconstant/(dx*dx) + driftvelocity*epsilon/dx;
		
		for(i=1;i<space-1;i++) {
		  coeff_m[i] = epsilon*diffusionconstant/(dx*dx) - 0.5*driftvelocity*epsilon/dx;
		  coeff_0[i] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (i-space0)*dx*epsilon - rho*epsilon;
		  coeff_p[i] = epsilon*diffusionconstant/(dx*dx) + 0.5*driftvelocity*epsilon/dx;
		}
		coeff_m[space-1] = epsilon*diffusionconstant/(dx*dx) - driftvelocity*epsilon/dx;
		coeff_0[space-1] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (space-1-space0)*dx*epsilon - rho*epsilon;
		coeff_p[space-1] = 0.;
		break;
    case 2:	reproduction = reproduction_jumps;
		coeff_m[0] = 0.;
		coeff_0[0] = 1. + (0-space0)*dx*epsilon - rho*epsilon - mutationrate/dx*epsilon;
		coeff_p[0] = driftvelocity*epsilon/dx;
		coeff_j[0] = mutationrate/dx;
		
		for(i=1;i<space-1;i++) {
		  coeff_m[i] = 0.5*driftvelocity*epsilon/dx;
		  coeff_0[i] = 1. + (i-space0)*dx*epsilon - rho*epsilon - mutationrate/dx*epsilon;
		  coeff_p[i] = 0.5*driftvelocity*epsilon/dx;
		  coeff_j[i] = epsilon*mutationrate/dx;
		}
		coeff_m[space-1] = -driftvelocity*epsilon/dx;
		coeff_0[space-1] = 1. + (space-1-space0)*dx*epsilon - rho*epsilon - mutationrate/dx*epsilon;
		coeff_p[space-1] = 0.;
		coeff_j[space-1] = mutationrate/dx*epsilon;
		break;
    case 3:	reproduction = reproduction_expkernel;
		generate_exp_mutationkernel();
		integrated_mutationkernel = 0.;
		for(i=0;i<space;i++)integrated_mutationkernel += mutationkernel[i+space];
		integrated_mutationkernel *= dx*mutationrate*epsilon;
		coeff_m[0] = 0.;
		coeff_0[0] = 1. - space0*dx*epsilon - rho*epsilon - integrated_mutationkernel;
		coeff_p[0] = driftvelocity*epsilon/dx;
		for(i=1;i<space-1;i++) {
		  coeff_m[i] = -driftvelocity*epsilon/dx;
		  coeff_0[i] = 1. - (i-space0)*dx*epsilon - rho*epsilon - integrated_mutationkernel;
		  coeff_p[i] = driftvelocity*epsilon/dx;
		}
		coeff_m[space-1] = -driftvelocity*epsilon/dx;
		coeff_0[space-1] = 1. + (space-1-space0)*dx*epsilon - integrated_mutationkernel;
		coeff_p[space-1] = driftvelocity*epsilon/dx;
		break;
    case 4:	reproduction = reproduction_kernel;
		read_mutationkernel();
		integrated_mutationkernel = 0.;
		for(i=0;i<space_m;i++)integrated_mutationkernel += mutationkernel[i];
// 		integrated_mutationkernel -= mutationkernel[space0_m];
		integrated_mutationkernel *= dx*mutationrate*epsilon;
		coeff_m[0] = 0.;
		coeff_0[0] = 1. - space0*dx*epsilon - rho*epsilon - integrated_mutationkernel;
		coeff_p[0] = driftvelocity*epsilon/dx;
		for(i=0;i<space;i++) {
		  coeff_m[i] = -driftvelocity*epsilon/dx;
		  coeff_0[i] = 1. + (i-space0)*dx*epsilon - rho*epsilon - integrated_mutationkernel;
		  coeff_p[i] = driftvelocity*epsilon/dx;
		}
		coeff_m[space-1] = -driftvelocity*epsilon/dx;
		coeff_0[space-1] = 1. + (space-1-space0)*dx*epsilon - rho*epsilon - integrated_mutationkernel;
		coeff_p[space-1] = 0.;
		break;
  }
      
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
  if((mutationtype==3) || (mutationtype==4)) {
    free(mutationkernel);
    free(coeff_kernel);
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
    if(i%outputstep == 0) {
      update_integrated_moments();
      write_output(i);
      if(timeavg == 1)calc_timeavg();
      if(covar == 1)calc_covar();
    }
  }
  
  if(timeavg == 1)write_timeavg();
  if(covar == 1)write_covar();
  if(write_conf_to_file == 1)write_c();
  
  cleanup();
  
  
  return 0;
}
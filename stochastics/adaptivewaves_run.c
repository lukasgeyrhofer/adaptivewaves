// adaptivewaves_run.c
// Lukas Geyrhofer

// v0.01, 111006
// v0.02, 111017, switched from chi-square-noise to poisson noise

// v0.03, 111101, added periodic boundary condidtions

// v0.04, 120104, + drift velocity
// v0.04.1, 120110, +option -K (keep), add +1.0 to c[x=0] to keep population in oasis
// v0.04.2, 120217, added output of central peak and +- 5 bins

// v0.05, 120330, added new u=delta option, "-m delta -b BETA"

// v0.06, 120605, added new mode to load u from file, "-m file -u FILENAME"
// 		  added paramater to option "-K KEEP" which adds KEEP to the density at the oasis in one generation

// v0.07, 120619, minor updates to comply with documentation
//		  added some comments

// v0.08, 120730, implemented moment-closure at higher levels, parameter -n (at least for v=0)
//		  time-averages can now be calculated up to moment M, option -M
// v0.08.1, 120806, output "covariance" matrix <c(x)c(y)>, option -C FILENAME

// v1.00, 120827, switched to adaptive waves
//		  version now included in config file

// v1.01, 120906, binary file structure changed. can now include arbitrary number of integer and double parameters in file,
//		  first two entries give the numbers of parameters,
//		  then lattice is defined (dx, space, space0),
// 		  then int/double parameters,
//		  finally density.
//		  u* read from file now also has this file-format (0,0)
//		  c has format (0,2) - rho, starttime

// v1.02, 120907, c format (0,4): added diffusionconstant, driftvelocity in binary file
//		  -D and -V options not available anymore

#define CONFVERSION 1
#define sqrt2 1.4142135623

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double r;
double dx;
int space,space0;
int version = CONFVERSION;


double beta = 1.;

double driftvelocity = 0.;
double diffusionconstant = 1.0;

double starttime;

int mode = 4;
int moment_closure = 1;
int calc_moments = 2;
double eps = 0.2;
int maxsteps = 10000;
unsigned long randseed = 0;
int outputinterval = 100;

int quiet = 0;


int distance = 5;

int offset_output = 0;
int outputinterval_dx = 1;

char in_filename[128];
char out_filename[128];
int outputtofile = 0;

char ufilename[128];

double *c,*c_new;
double *u;
double *xpos_eps;


double **c_timeavg_sum;
double *c_timeavg_tmp,*c_timeavg_ones;
int timeavg = 0;
char timeavg_filename[128];
double n_timeavg = 0.;

double **covariance;
char covar_filename[128];
double n_covar = 0.;
int covar = 0;

double r_sqrt;
double N_avg_theoretical;
double N;
double lambda = 1.;
double xcae;



double N_constraint = 0.;

const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


// *******************************
// get all options from the command line.
// provived with flags, i.e. "-i FILENAME". see list of letters in the
// GETOPT arguments. letters with a semicolon after them require an
// additional value, see "-i FILENAME", or "-S MAXSTEPS". if flag is not given,
// standard value of variables is taken, see declaration above.
// *******************************


void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn, argv,"i:o:e:S:R:m:O:d:qQt:N:b:u:M:C:")) != -1){
    switch(c) {
      case 'i':	strcpy(in_filename,optarg);
		haveinfile = 1;
		break;
      case 'o': strcpy(out_filename,optarg);
		if(strcmp(out_filename,"-") == 0) {
		  outputtofile = 2;
		}else{
		  outputtofile = 1;
		}
		break;
      case 'e': eps = atof(optarg);
		break;
      case 'S': maxsteps = atoi(optarg);
		break;
      case 'R': randseed = atoi(optarg);
		break;
      case 'm': if(strcmp(optarg,"flat") == 0) {
		  mode = 2;
		}else if(strcmp(optarg,"delta") == 0) {
		  mode = 3;
		}else if(strcmp(optarg,"file") == 0) {
		  mode = 4;
		}else print_error("mode not known!");
		break;
      case 'b':	if(mode == 3) {
		  beta = atof(optarg);
		}else{
		  print_error("option b only in mode -m \"delta\"");
		}
		break;
      case 'u':	if(mode==4) {
		  strcpy(ufilename,optarg);
		}else{
		  print_error("option u only in mode -m \"file\"");
		}
		break;
      case 'N': if(mode!=2) {
		  print_error("need mode '-m flat' to use option N for populations size control");
		}else{
		  N_constraint = atof(optarg);
		}
		break;
      case 'O':	outputinterval = atoi(optarg);
		break;
      case 'd': outputinterval_dx = atoi(optarg);
		break;
      case 'q':	quiet = 1;
		break;
      case 'Q': quiet = 2;
		break;
      case 't': strcpy(timeavg_filename,optarg);
		timeavg = 1;
		break;
      case 'M': if(timeavg == 1) {
		  calc_moments = atoi(optarg);
		}else{
		  print_error("option M only with time-averages, option -t FILENAME");
		}
		break;
//       case 'V': driftvelocity = atof(optarg);
// 		break;
//       case 'D':	diffusionconstant = atof(optarg);
// 		break;
      case 'C': strcpy(covar_filename,optarg);
		covar = 1;
		break;
    }
    
  }
  if(haveinfile == 0)print_error("infile necessary! use option -i");
  if(randseed <= 0)randseed = time(NULL);
  if(mode == 2) {
    if(N_constraint <= 0.) {
      print_error("need to set population size for mode '-m flat'");
    }
  }else if(mode == 4) {
    if(strcmp(ufilename,"") == 0) {
      print_error("need filename of u*-file in mode '-m file'");
    }
  }
  
}

// *******************************
// write configuration file at end of run
// *******************************

void write_data() {
  FILE *fp;
  int tmp_count_params[2];
  double curtime = starttime + eps * maxsteps;
  tmp_count_params[0] = 0;
  tmp_count_params[1] = 4;
  if(outputtofile == 2) {
    fp = stdout;
  }else{
    fp = fopen(out_filename,"wb");
  }
  fwrite(&tmp_count_params[0],2,sizeof(int),fp);
  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  fwrite(&r,1,sizeof(double),fp);
  fwrite(&curtime,1,sizeof(double),fp);
  fwrite(&diffusionconstant,1,sizeof(double),fp);
  fwrite(&driftvelocity,1,sizeof(double),fp);
  fwrite(&c[0],space,sizeof(double),fp);
  if(outputtofile == 1)fclose(fp);
}


// *******************************
// read configuration file
// *******************************

int read_data() {
  FILE *fp;
  int usefile = 0;
  int tmp_count_params[2];
  int i,tmpi;
  double tmpd;
  if(strcmp(in_filename,"-") == 0) { // input comes from STDIN, i.e. directly from the console.
				     // usage: "./oasis_run -i - [OTHEROPTIONS] < start.conf"
    fp = stdin;
  }else{
    fp=fopen(in_filename,"rb");
    if(fp == NULL)print_error("input file not found!");
    usefile = 1;
  }
  fread(&tmp_count_params[0],2,sizeof(int),fp);
  if(tmp_count_params[1]<4)print_error("configuration file does not contain necessary data.");
    
  fread(&dx,1,sizeof(double),fp);
  fread(&space,1,sizeof(int),fp);
  fread(&space0,1,sizeof(int),fp);

  for(i=0;i<tmp_count_params[0];i++)fread(&tmpi,1,sizeof(int),fp); // forward compatibility
  fread(&r,1,sizeof(double),fp);
  fread(&starttime,1,sizeof(double),fp);
  fread(&diffusionconstant,1,sizeof(double),fp);
  fread(&driftvelocity,1,sizeof(double),fp);
  for(i=4;i<tmp_count_params[1];i++)fread(&tmpd,1,sizeof(double),fp); // forward compatibility

  c=(double*)malloc(space*sizeof(double));
  fread(&c[0],space,sizeof(double),fp);
  if(usefile == 1)fclose(fp);
}

// *******************************
// read u(x) from external file
// note that lattice has to be the same, otherwise the simulation runs on random u(x)
// *******************************

int get_u_from_file() {
  FILE *fpu;
  double ux;
  int i,tmp_count_params[2];
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(ufilename,"rb");
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
    fread(&u[0],space,sizeof(double),fpu);
    
    fclose(fpu);
  }
  return 0;
}



// *******************************
// calculate COM
// *******************************

double get_mean_pos(double popsize) {
  int i;
  double meanpos = 0.0;
  for(i=0;i<space;i++) {
    meanpos += (i-space0)*c[i];
  }
  return meanpos*dx*dx/popsize;
}


// *******************************
// \int dx u(x)c(x,t)
// *******************************

double get_lambda() {
  int i;
  double sum=0.;
  for(i=0;i<space;i++) {
    sum += c[i]*u[i];
  }
  return sum*dx;
}


// *******************************
// initialization
// *******************************

void initialize() {
  int i,j;
  double x;
  
  // setup random number generator
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  
  u = (double*)malloc(space*sizeof(double));
  c_new = (double*)malloc(space*sizeof(double));
  xpos_eps = (double*)malloc(space*sizeof(double));
  
  for(i=0;i<space;i++) {
    xpos_eps[i] = (i-space0)*dx*eps;
//     printf("xpos_eps[%06d] = %lf\n",i,xpos_eps[i]);
  }
  
//   exit(1);
  
  if(timeavg == 1) {
    c_timeavg_sum = (double**)malloc(calc_moments*sizeof(double*));
    for(j=0;j<calc_moments;j++) {
      c_timeavg_sum[j] = (double*)malloc(space*sizeof(double));
      if(j==0) {
	c_timeavg_ones = (double*)malloc(space*sizeof(double));
	c_timeavg_tmp = (double*)malloc(space*sizeof(double));
      }
      for(i=0;i<space;i++) {
	c_timeavg_ones[i] = 1.;
	c_timeavg_sum[j][i]  = 0.;
      }
    }
    n_timeavg = 0;
  }
  
  if(covar == 1) {
    covariance = (double**)malloc(space*sizeof(double*));
    for(i=0;i<space;i++) {
      covariance[i] = (double*)malloc(space*sizeof(double));
      for(j=0;j<space;j++) {
	covariance[i][j] = 0.;
      }
    }
  }
    
  
  if(mode == 4) {
    // mode 4 loads u(x) from external file
    get_u_from_file();
  }
  
//   double sum =0.;
//   for(i=0;i<space;i++) {
//     sum+=c[i]*u[i]*dx;
//     printf("%lf %g %g %g\n",(i-space0)*dx,c[i],u[i],sum);
//   }
//   exit(1);
  
  
  offset_output = space0 % outputinterval_dx;
  
}

// *******************************
// reproduction step
// *******************************

void reproduction() {
  int i;
  double dw;
  
  // preset all coefficients
  double epsdx2i = diffusionconstant*eps/dx/dx;
  double determchangefactor = 1. - 2.*epsdx2i -eps*r;
  double eps_sqrt_2 = sqrt2*sqrt(eps)*sqrt(dx);
  double vepshalfdxi = 0.5*eps*driftvelocity/dx;
  
  
  
  printf("vepshalfdxi = %g, determchangefactor= %g, epsdx2i = %g\n",vepshalfdxi,determchangefactor,epsdx2i);
  exit(1);
  
//   printf("determchangefactor = %14.10lf, eps_sqrt_2 = %14.10lf, vepshalfdxi = %14.10lf, epsdx2i = %14.10lf, r=%g\n",determchangefactor,eps_sqrt_2,vepshalfdxi,epsdx2i,r);
//   exit(1);
  
  
//   printf("%lf\n",vepshalfdxi);
  
  int last = space-1;
  
  // update occupancies for all positions, first in variable c_new which is then copied to c
  // first and last index are treated separately to implement periodic boundaries

//   c_new[0] = c[0] * (determchangefactor + xpos_eps[0]);
//   c_new[0] += 2.*c[1]*epsdx2i;
//   dw = gsl_ran_poisson(rg,c_new[0])-c_new[0];
//   c_new[0] += eps_sqrt_2*dw;
//   
  c_new[0] = 0.;
  

  for(i=1;i<space-1;i++) {
    c_new[i]  = c[i] * (determchangefactor + xpos_eps[i]);
    c_new[i] += (c[i+1]+c[i-1])*epsdx2i;
    c_new[i] += (c[i+1]-c[i-1])*vepshalfdxi;
//     dw = gsl_ran_poisson(rg,c_new[i])-c_new[i];
//     c_new[i] += eps_sqrt_2*dw;
  }
  c_new[last] = 0.;
//   c_new[last] = c[last] * (determchangefactor + xpos_eps[last]);
//   c_new[last] += 2.*c[last-1]*epsdx2i;
//   dw = gsl_ran_poisson(rg,c_new[last])-c_new[last];
//   c_new[last] += eps_sqrt_2*dw;
  
  // memcpy copies memory blocks. fastest way to copy
  memcpy(&c[0],&c_new[0],space*sizeof(double));
}


// *******************************
// \int dx c(x,t)
// *******************************

double get_N() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++)sum += c[i];
  return sum*dx;
}


// *******************************
// different modes for population control
// *******************************

void population_control() {
  double l,li;
  int i;
  if(mode == 2) {
    // fixed N
    lambda = get_N();
    lambda /= N_constraint;
  }else if(mode == 3) {
    // only occupancy at oasis is kept fixed
    lambda = c[space0]/beta;
  }else if(mode == 4) {
    // still contraint <u|c>=1, but u(x) is from external file
    lambda = get_lambda();
  }
  li = 1./lambda;
  // each occupancy is globally rescaled
  for(i=0;i<space;i++)c[i] *= li;
}



void calc_timeavg() {
  int i,j;
  memcpy(&c_timeavg_tmp[0],&c_timeavg_ones[0],space*sizeof(double));
  for(j=0;j<calc_moments;j++) {
    for(i=0;i<space;i++) {
      c_timeavg_tmp[i] *= c[i];
      c_timeavg_sum[j][i] += c_timeavg_tmp[i];
    }
  }
  n_timeavg +=1.;
}

void write_timeavg() {
  FILE *fp;
  int i,j;
  fp = fopen(timeavg_filename,"w");
  fprintf(fp,"# simulationtime = %lf\n",maxsteps*eps);
  fprintf(fp,"# datapoints = %d\n",n_timeavg);
  for(i=0;i<space;i++) {
    fprintf(fp,"%lf",(i-space0)*dx);
    for(j=0;j<calc_moments;j++) {
      fprintf(fp,"\t%e",c_timeavg_sum[j][i]/n_timeavg);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}


void calc_covariance() {
  int i,j;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      covariance[i][j] += c[i]*c[j];
    }
  }
  n_covar += 1.;
}


void write_covariance() {
  FILE *fp;
  int i,j;
  fp=fopen(covar_filename,"w");
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      fprintf(fp,"%lf %lf %e",(i-space0)*dx,(j-space0)*dx,covariance[i][j]/n_covar);
      if(timeavg == 1) {
	fprintf(fp," %e",covariance[i][j]/n_covar-c_timeavg_sum[0][i]*c_timeavg_sum[0][j]/(n_timeavg*n_timeavg));
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void print_data(int step) {
  int i;
  double curtime = starttime + eps*step;
  double meanpos;
  if(outputtofile != 2) {
    N = get_N();
    meanpos = get_mean_pos(N);
    fprintf(stdout,"%lf %lf %lf %lf %14.12lf\n",curtime,N,lambda,meanpos,1.-1./lambda);
  }
  if(quiet!=2) {
    for(i=0;i<space;i++) {
      if( (i) %outputinterval_dx == 0)fprintf(stderr,"%lf %lf %e\n",curtime,(i-space0)*dx,c[i]);
    }
    fprintf(stderr,"\n");
  }
}

void print_options() {
  fprintf(stdout,"# running stochastic simulations of adaptive waves...\n");
  fprintf(stdout,"# options read from infile:\n");
  fprintf(stdout,"#   rho       = %6.4e\n",r);
  fprintf(stdout,"#   dx        = %lf\n",dx);
  fprintf(stdout,"#   starttime = %lf\n",starttime);
  fprintf(stdout,"#   space     = %d\n",space);
  fprintf(stdout,"#   space0    = %d\n",space0);
  fprintf(stdout,"# options provided for running stochastic simulation:\n");
  fprintf(stdout,"#   epsilon   = %lf\n",eps);
  fprintf(stdout,"#   v         = %e\n",driftvelocity);
  fprintf(stdout,"#   DiffConst = %e\n",diffusionconstant);
  fprintf(stdout,"#   maxsteps  = %d\n",maxsteps);
  fprintf(stdout,"#   output_t  = %d\n",outputinterval);
  fprintf(stdout,"#   output_x  = %d\n",outputinterval_dx);
  fprintf(stdout,"#   mode      = ");
  switch(mode) {
    case 1:fprintf(stdout,"u\n");break;
    case 2:fprintf(stdout,"flat\n");break;
    case 3:fprintf(stdout,"delta\n");break;
    case 4:fprintf(stdout,"file\n");break;
  }
}



void cleanup() {
  int i,j;
  if(timeavg == 1) {
    for(j=0;j<calc_moments;j++) {
      free(c_timeavg_sum[j]);
    }
    free(c_timeavg_ones);
    free(c_timeavg_sum);
    free(c_timeavg_tmp);
  }
  free(c);
  free(c_new);
  free(u);
  free(xpos_eps);
}

int main(int argn, char* argv[]) {
  int i,j;
  
  parsecommandline(argn, argv);

  read_data();
  initialize();
  
//   double popsizetmp = get_N();
//   double tmp = get_mean_pos(popsizetmp);
//   printf("com = %lf\n",tmp);
//   exit(1);
  
  if(quiet!=1)print_options();
  if(quiet!=1) {
    fprintf(stdout,"\n\n# time N lambda com r\n");
    print_data(0);
  }
  for(i=1;i<=maxsteps;i++) {
//     fprintf(stderr,"before reproduction!\n");
    reproduction();
//     fprintf(stderr,"before control!\n");
    population_control();
//     fprintf(stderr,"before output!\n");
    if(i%outputinterval == 0) {
      if(quiet!=1)print_data(i);
      if(timeavg==1)calc_timeavg();
      if(covar==1)calc_covariance();
    }
  }
  
  if(outputtofile >= 1)write_data();
  if(timeavg == 1)write_timeavg();
  if(covar == 1)write_covariance();
  cleanup();
  return 0;
}

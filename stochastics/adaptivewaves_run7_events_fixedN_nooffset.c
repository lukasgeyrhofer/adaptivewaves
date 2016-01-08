// Stochastic simulation of adaptation in the high-mutation, strong-selection regime
// based on travelling waves
// Lukas Geyrhofer, 2012-2013

// version 0.1 (first named version)
//	130130		included different mutation models

//	130513		including subpopulation w
//	130521		stop-condition based on #events (fixation/extinction)


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
int maxEvents = 1000;
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
double *n_start;		// fallback density after fixation/extinction event
double *n_tmp;			// for shift_windows procedure
double *c;			// in- and output is still based on densities c[i], (n[i] = c[i]*dx)
double *u;
double *coeff_m,*coeff_0,*coeff_p;
double *zeros;

double populationsize;


// labelled subpopulations

double *w;
double *nw,*nw_new, *nw_start;
double *nv,*nv_new, *nv_start;
int subpop_inject = 0;
double subpop_fixationthreshold = 1e-3;
double subpop_extinctionthreshold = 1e-8;
double subpop_labelthreshold = 0.;
double subpop_expectedfixationprobability = 0.;

int subpop_identicalstartpopulation = 1;


char winfilename[128];
char woutfilename[128];
int subpop_active = 0;
int subpop_fromfile = 0;
int subpop_label = 0;
int subpop_haveoutfile = 0;
double subpop_fraction = 0.;

int subpop_extinctions = 0;
int subpop_fixations = 0;

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

int mutation_type = 1;	// 1 - diffusion
			//      -m diffusion
			// 2 - single (beneficial) jumps
			//      -m jumps
			// 3 - exponential kernel, only beneficial
			// 4 - exponential kernel, beneficial AND deleterious
			// 5 - kernel from external file
double diffusionconstant;


// shift window...

int shifthreshold = 10;
int currentshift = 0;




// function pointers
typedef void (*iteration_fp)();
iteration_fp reproduction;
iteration_fp populationcontrol;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

void parsecommandline(int argn, char* argv[]) {
  char c;
  int havecinfile = 0, haveuinfile = 0, havepopulationsize = 0;
  while((c = getopt(argn,argv,"c:C:u:e:S:O:R:t:M:qQv:V:m:L:N:T:")) != -1) {
    switch(c) {
      case 'c':	strcpy(cinfilename,optarg);
		havecinfile = 1;
		break;
      case 'C': strcpy(coutfilename,optarg);
		write_conf_to_file = 1;
		break;
      case 'u': strcpy(uinfilename,optarg);
		haveuinfile = 1;
		break;
      case 'e':	epsilon=atof(optarg);
            	break;
      case 'S': maxEvents = atoi(optarg);
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
      case 'm':	if(strcmp(optarg,"diffusion") == 0) {
                  mutation_type = 1;
                }else if(strcmp(optarg,"jumps") == 0) {
                  mutation_type = 2;
                }else{
                  print_error("mutation type not known!");
                }
                break;
		
      case 'L':	subpop_expectedfixationprobability = atof(optarg);
		subpop_active = 1;
		break;
      case 'N':	populationsize = atof(optarg);
		havepopulationsize = 1;
		break;
      case 'T':	shifthreshold = atoi(optarg);
		break;
    }
  }
  if(havecinfile + haveuinfile + havepopulationsize < 3)print_error("need all options -c, -u and -N");
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

void write_c(int step) {
  int i;
  FILE *fp;
  int icount = 0,dcount = 4;
  double endtime = starttime + epsilon*step;
  
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



void write_w() {
  int icount = 0,dcount = 0;
  int i;
  FILE *fpw;
  fpw = fopen(woutfilename,"wb");
  if(fpw==NULL)print_error("could not open subpop-outfile");
  fwrite(&icount,1,sizeof(int),fpw);
  fwrite(&dcount,1,sizeof(int),fpw);
  fwrite(&dx,1,sizeof(double),fpw);
  fwrite(&space,1,sizeof(int),fpw);
  fwrite(&space0,1,sizeof(int),fpw);
  for(i=0;i<space;i++)w[i] = nw[i]/dx;
  fwrite(&w[0],space,sizeof(double),fpw);
}


void read_w(int allocatemem) {
  int i;
  int icount,dcount;
  int *ival;
  double *dval;
  double w_dx;
  int w_space,w_space0;
  FILE *fpw;
  fpw = fopen(winfilename,"rb");
  if(fpw==NULL)print_error("could not open subpop-infile");
  fread(&icount,1,sizeof(int),fpw);
  fread(&dcount,1,sizeof(int),fpw);
  fread(&w_dx,1,sizeof(double),fpw);
  fread(&w_space,1,sizeof(int),fpw);
  fread(&w_space0,1,sizeof(int),fpw);
  if(icount > 0) {
    ival = (int*)malloc(icount*sizeof(int));
    fread(&ival[0],icount,sizeof(int),fpw);
    free(ival);
  }
  if(dcount > 0) {
    dval = (double*)malloc(dcount*sizeof(double));
    fread(&dval[0],dcount,sizeof(double),fpw);
    free(dval);
  }
  if((space != w_space) || (space0 != w_space0) || (w_dx/dx > 1.1) || (w_dx/dx < 0.9))print_error("lattice does not match (subpop-infile)");
  if(allocatemem == 1) {
    w  = (double*)malloc(space*sizeof(double));
    nw = (double*)malloc(space*sizeof(double));
  }
  fread(&w[0],space,sizeof(double),fpw);
  for(i=0;i<space;i++)nw[i] = w[i]*dx;
  fclose(fpw);
}


void check_w() {
  int i;
  for(i=0;i<space;i++) {
    if(w[i] > c[i])print_error("subpopulation larger than normal population");
  }
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


void write_timeavg(int step) {
  FILE *fp;
  int i,j;

  fprintf(stdout,"# writing time-averaged moments to '%s'\n",timeavg_filename);
  
  fp = fopen(timeavg_filename,"w");
  fprintf(fp,"# simulationtime = %lf\n",step*epsilon);
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



int get_offset() {
  int i;
  int currentoffset = 0;
  double int_uc = 0., last_int_uc;
  
  for(i=0;i<space;i++)int_uc += u[i]*n[i];
  
//   printf("first try: int_uc = %g, offset = %d\n",int_uc,currentoffset);
  
  if(int_uc < 1.) {
    while(int_uc < 1.) {
      currentoffset++;
      last_int_uc = int_uc;
      int_uc = 0.;
      for(i=currentoffset;i<space-currentoffset;i++)int_uc += u[i+currentoffset]*n[i];
//       printf("int_uc = %g, offset = %d\n",int_uc,currentoffset);
    }
    // going up in int_uc, so last one (int_uc) would be larger than 1., last_int_uc would be smaller than 1.
    // check (absolute) differences, take offset with smallest difference to 1.
    if(int_uc - 1. > 1. - last_int_uc)currentoffset--;
  }else if(int_uc > 1.) {
    while(int_uc > 1.) {
      currentoffset--;
      last_int_uc = int_uc;
      int_uc = 0.;
      for(i=currentoffset;i<space-currentoffset;i++)int_uc += u[i+currentoffset]*n[i];
//       printf("int_uc = %g, offset = %d\n",int_uc,currentoffset);
    }
    // going down in int_uc, int_uc<1, last_int_uc>1
    if(1. - int_uc > last_int_uc - 1.)currentoffset++;
  }
  fprintf(stdout,"# offset = %d\n",currentoffset);
  return currentoffset;
}



 
int get_fitness_bin() {
  int i;
  double cumulative_g = 0.,tmp;
  int offset = 0;
  
//   offset = get_offset();
  
  for(i=space-offset-1;i>offset;i--) {
    cumulative_g += u[i+offset]*n[i];
    if(cumulative_g > subpop_expectedfixationprobability)break;
  }
  tmp = cumulative_g - u[i+offset]*n[i];
  fprintf(stdout,"# parameter fixation threshold: %lf\n# real fitness threshold: %lf %lf\n",subpop_expectedfixationprobability,cumulative_g,tmp);
  return i;
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
 
void reproduction_nosubpop() {
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


void reproduction_subpop() {
  int i;
  int last = space-1;
  double rhoepsilon = rho*epsilon;
  
  for(i=1;i<space-1;i++) {
    nv_new[i] = coeff_m[i]*nv[i-1] + (coeff_0[i] - rhoepsilon)*nv[i] + coeff_p[i] * nv[i+1];
    nv_new[i] += twoepssqrt*(gsl_ran_poisson(rg,nv[i])-nv[i]);
    nw_new[i] = coeff_m[i]*nw[i-1] + (coeff_0[i] - rhoepsilon)*nw[i] + coeff_p[i] * nw[i+1];
    nw_new[i] += twoepssqrt*(gsl_ran_poisson(rg,nw[i])-nw[i]);
  }
  memcpy(&nv[1],&nv_new[1],(space-2)*sizeof(double));
  memcpy(&nw[1],&nw_new[1],(space-2)*sizeof(double));
}

double get_int_c_subpop() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++)sum+=(nv[i]+nw[i]);
  return sum;
}

double update_population() {
  int i;
  for(i=1;i<space-1;i++)n[i] = nv[i] + nw[i];
}
  

double get_int_c_nosubpop() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++)sum+=u[i]*n[i];
  return sum;
}


void update_integrated_moments() {
  int i;
  double x;
  if(subpop_active==1)update_population();
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

void populationcontrol_nosubpop() {
  int i;
  double lambda_inv;
  lambda = get_int_c_nosubpop()/populationsize;
  lambda_inv = 1./lambda;
  for(i=0;i<space;i++)n[i] *= lambda_inv;
}

void populationcontrol_subpop() {
  int i;
  double lambda_inv;
  lambda = get_int_c_subpop()/populationsize;
  lambda_inv = 1./lambda;
  for(i=1;i<space-1;i++) {
    nv[i] *= lambda_inv;
    nw[i] *= lambda_inv;
  }
}

double get_subpop_fraction() {
  int i;
  double nsum = 0.,nwsum = 0.;
  for(i=0;i<space;i++) {
    nsum += nw[i]+nv[i];
    nwsum += nw[i];
  }
  return nwsum/nsum;
}




void write_output(int step) {
  int i;
  if(quiet < 2) {
    fprintf(stdout,"%10.4lf %16.3lf %16.10lf %16.10lf %16.10lf",epsilon*step,populationsize,dens_mom1/dens_mom0,dens_mom2/dens_mom0,lambda);
    if(subpop_active == 1) {
      subpop_fraction = get_subpop_fraction();
      fprintf(stdout,"%16.10lf",subpop_fraction);
    }
    fprintf(stdout,"\n");
    if(quiet < 1) {
      for(i=0;i<space;i++) {
	fprintf(stderr,"%lf %lf %g %g\n",epsilon*step,(i-space0)*dx,n[i]/dx,nw[i]/dx);
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
    fprintf(stdout,"#     maxEvents      = %d\n",maxEvents);
    fprintf(stdout,"#     outputstep     = %d\n",outputstep);
    fprintf(stdout,"#     cfile          = %s\n",cinfilename);
    fprintf(stdout,"#     ufile          = %s\n",uinfilename);
    fprintf(stdout,"#     populationsize = %g\n",populationsize);
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

void relabel_population(int step, int ef) {
  int i,j;
  int label_startindex;
  
  // update complete profile
  if(ef == 0) {
    fprintf(stdout,"# exctintion time %lf\n",step*epsilon);
  }else if(ef == 1) {
    fprintf(stdout,"# fixation time %lf\n",step*epsilon);
  }
    
  memcpy(&nw[0],&nw_start[0],space*sizeof(double));
  memcpy(&nv[0],&nv_start[0],space*sizeof(double));
  memcpy(&n[0], &n_start[0], space*sizeof(double));
  
  subpop_fraction = get_subpop_fraction();
}    



void create_subpopulation() {
  int i,j;
  int label_startindex;
  
  printf("# create subpopulation\n");
  
  
  if(subpop_active == 1) {
    w =      (double*)malloc(space*sizeof(double));
    nw =     (double*)malloc(space*sizeof(double));
    nw_new = (double*)malloc(space*sizeof(double));
    nv =     (double*)malloc(space*sizeof(double));
    nv_new = (double*)malloc(space*sizeof(double));
    
    label_startindex = get_fitness_bin();
    if(label_startindex > space)print_error("subpop label threshold larger than simulationbox");
    memcpy(&nw[0],&zeros[0],label_startindex*sizeof(double));
    memcpy(&nw[label_startindex],&n[label_startindex],(space-label_startindex)*sizeof(double));
    memcpy(&nv[0],&n[0],label_startindex*sizeof(double));
    memcpy(&nv[label_startindex],&zeros[0],(space-label_startindex)*sizeof(double));
  }

  n[0]  = 0.;  n[space-1]  = 0.;
  nw[0] = 0.;  nw[space-1] = 0.;
  nv[0] = 0.;  nv[space-1] = 0.;

  nw_start = (double*)malloc(space*sizeof(double));
  nv_start = (double*)malloc(space*sizeof(double));
  n_start  = (double*)malloc(space*sizeof(double));
  memcpy(&nw_start[0],&nw[0],space*sizeof(double));
  memcpy(&nv_start[0],&nv[0],space*sizeof(double));
  memcpy(&n_start[0], &n[0], space*sizeof(double));
  
  subpop_fraction = get_subpop_fraction();
}



void shift_window() {
  double meanfit = dens_mom1/dens_mom0/dx;
  int actualshift;
  if(meanfit > shifthreshold) {
    actualshift = (int)floor(meanfit);

    memcpy(&n_tmp[0],&nw[actualshift],(space-actualshift)*sizeof(double));
    memcpy(&n_tmp[space-actualshift],&zeros[0],actualshift*sizeof(double));
    memcpy(&nw[0],&n_tmp[0],space*sizeof(double));

    memcpy(&n_tmp[0],&nv[actualshift],(space-actualshift)*sizeof(double));
    memcpy(&n_tmp[space-actualshift],&zeros[0],actualshift*sizeof(double));
    memcpy(&nv[0],&n_tmp[0],space*sizeof(double));
    
    currentshift += actualshift;
    rho -= actualshift*dx;
  }else if(meanfit < -shifthreshold) {
    actualshift = -(int)ceil(meanfit);

    memcpy(&n_tmp[actualshift],&nw[0],(space-actualshift)*sizeof(double));
    memcpy(&n_tmp[0],&zeros[0],actualshift*sizeof(double));
    memcpy(&nw[0],&n_tmp[0],space*sizeof(double));

    memcpy(&n_tmp[actualshift],&nv[0],(space-actualshift)*sizeof(double));
    memcpy(&n_tmp[0],&zeros[0],actualshift*sizeof(double));
    memcpy(&nv[0],&n_tmp[0],space*sizeof(double));

    currentshift -= actualshift;
    rho += actualshift*dx;
  }
}





void initialize() {
  int i,j;
  double x;
  
  read_c();
  read_u();
  
  print_options();
  
  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  n_new   = (double*)malloc(space*sizeof(double));
  n_tmp   = (double*)malloc(space*sizeof(double));
  
  coeff_m = (double*)malloc(space*sizeof(double));
  coeff_0 = (double*)malloc(space*sizeof(double));
  coeff_p = (double*)malloc(space*sizeof(double));
  
  zeros = (double*)malloc(space*sizeof(double));
  
  
  coeff_m[0] = 0.;
  coeff_0[0] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (0-space0)*dx*epsilon - rho*epsilon;
  coeff_p[0] = epsilon*diffusionconstant/(dx*dx);// + driftvelocity*epsilon/dx;
  zeros[0]   = 0.;
  
  for(i=1;i<space-1;i++) {
    coeff_m[i] = epsilon*diffusionconstant/(dx*dx); // - 0.5*driftvelocity*epsilon/dx;
    coeff_0[i] = 1. - 2.*epsilon*diffusionconstant/(dx*dx) + (i-space0)*dx*epsilon - rho*epsilon;
    coeff_p[i] = epsilon*diffusionconstant/(dx*dx);// + 0.5*driftvelocity*epsilon/dx;
    zeros[i]   = 0.;
  }
  coeff_m[space-1] = epsilon*diffusionconstant/(dx*dx);// - driftvelocity*epsilon/dx;
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
  
  if(subpop_active == 1) {
    create_subpopulation();
    reproduction = reproduction_subpop;
    populationcontrol = populationcontrol_subpop;
  }else{
    print_error("need fixation threshold");
    reproduction = reproduction_nosubpop;
    populationcontrol = populationcontrol_nosubpop;
  }
}



void cleanup() {
  int i;
  free(c);
  free(n);
  free(n_new);
  free(n_tmp);
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
  if(subpop_active==1) {
    free(w);
    free(nw);
    free(nv);
    free(nv_new);
    free(nw_new);
    free(n_start);
    free(nw_start);
    free(nv_start);
  }
}


int main(int argn, char* argv[]) {
  int eventcount=0,i=0;
  
  parsecommandline(argn,argv);
  
  initialize();
  
  update_integrated_moments();
  write_output(0);
  
  while(eventcount < maxEvents) {
    i++;
    reproduction();
    update_integrated_moments();
    rho = dens_mom1/dens_mom0;
    populationcontrol();
    shift_window();

    if(i%outputstep==0) {
      write_output(i);
      if(timeavg==1)calc_timeavg();
      if(covar==1)calc_covar();
    }
    if(subpop_fraction < subpop_extinctionthreshold) {
      subpop_extinctions++;
      eventcount++;
      relabel_population(i,0);
    }
    if(1. - subpop_fraction < subpop_fixationthreshold) {
      subpop_fixations++;
      eventcount++;
      relabel_population(i,1);
    }
    
  }

  if(timeavg==1)write_timeavg(i);
  if(covar==1)write_covar();
  if(write_conf_to_file == 1)write_c(i);
  if(subpop_haveoutfile == 1)write_w();
  
  printf("# extinctions: %d\n# fixations: %d\n",subpop_extinctions,subpop_fixations);
  
  cleanup();
  
  return 0;
}

/* #################################################
 * # stochastic simulation of adapting populations #
 * # with variing constraint type                  #
 * # and diffusion mutation kernel                 #
 * #################################################
 * 
 * 2012-2015, Lukas Geyrhofer
 * 
 * 
 * #################################################
 * 
 * simplest usage (population size 1e7):
 * ./travelingwavepeak_diff -N 1e7 1> out.txt 2> conf.txt
 * 
 * ################################################# */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>


int space = 300;
int space0 = 100;
double dx = 1e-2;
int maxSteps = 1000;

int outputstep = 100;
int quiet = 0;

int noise = 0;
double populationsize = 1.;
double populationvariance;
double xg;
int correctformeanfitness =0;
int printhistotype = 0;

int windcontrolled = 0;
double dv = 0.;

int allshifts = 0;
int shiftthreshold = 1;
double current_mean_fitness = 0.;

int read_from_file = 0;
int write_to_file = 0;
char c_infile[128],c_outfile[128];

int have_u_infile = 0;
char u_infile[128];
double *u_read,*u;
int dens_ustar_latticeratio = 1;
int u_space,u_space0;
double u_dx;
double wavespeed = 0.;
double speedprefactor;


double epsilon = 1e-2;
double twoepssqrt;
double *nn;
double *tmp;
double *x;

double mutationrate = 1e-5;
double *mutation_inflow;
double mutation_outflow;
// double mutation_sigma = 1e-2;

double mutation_jump = 0.;

const gsl_rng* rg; /* gsl, global generator */
const gsl_rng_type* T;
unsigned long int randseed = 0;


double popdens_0thmom;
double popdens_1stmom;
double popdens_2ndmom;
double ucp;


int averagepopdens = 0;
int averagepopdens_center = 1;
int averagepopdens_resolution = 1;
int averagepopdens_havefile = 0;
int averagepopdens_havecovarfile = 0;
double *averagepopdens_dens;
double *averagepopdens_densnose;
double **averagepopdens_covar;
double **averagepopdens_covarnose;
char averagepopdens_outputfile[128];
char averagepopdens_covarfile[128];
double averagepopdens_count = 0.;
double averagepopdens_covarcount = 0.;
double averagepopdens_dx;
int averagepopdens_space;
int averagepopdens_space0;
int averagepopdens_lower;
int averagepopdens_higher;

double *timetrace_popsize;
char timetrace_filename[128];
int timetrace_ps_interval = 100;
int timetrace_length;


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

// ************************************************************
// **   parameters
// ************************************************************


void parsecomamndline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn,argv,"s:z:d:S:e:D:O:qQi:o:R:N:Cu:U:T:a:A:PWn:g:")) != -1) {
    switch(c) {
      case 's':	space = atoi(optarg);
		break;
      case 'z':	space0 = atoi(optarg);
		break;
      case 'd':	dx = atof(optarg);
		break;
      case 'i':	strcpy(c_infile,optarg);
		read_from_file = 1;
		break;
      case 'o':	strcpy(c_outfile,optarg);
		write_to_file = 1;
		break;
      case 'u':	strcpy(u_infile,optarg);
		if(noise > 0)
		  print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME or -n FILENAME_TIMETRACE)");
		noise = 2;
		break;
      case 'U':	if(noise == 2) {
		  dens_ustar_latticeratio = atoi(optarg);
		}else{
		  print_error("option -u FILENAME needed before option -U RATIO");
		}
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
      case 'R':	randseed =atoi(optarg);
		break;
      case 'N': populationsize = atof(optarg);
		if(noise > 0)
		  print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME or -n FILENAME_TIMETRACE)");
		noise = 1;
		break;
      case 'C':	correctformeanfitness = 1;
		break;
//       case 'M':	mutation_sigma = atof(optarg);
// 		break;
      case 'T':	shiftthreshold = atoi(optarg);
		break;
      case 'a':	strcpy(averagepopdens_outputfile,optarg);
		averagepopdens_havefile = 1;
		averagepopdens = 1;
		break;
      case 'A': strcpy(averagepopdens_covarfile,optarg);
		averagepopdens_havecovarfile = 1;
		averagepopdens = 1;
		break;
      case 'P':	printhistotype = 1;
		break;
      case 'W':	windcontrolled = 1;
		break;
      case 'n':	strcpy(timetrace_filename,optarg);
		if(noise > 0)print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME or -n FILENAME_TIMETRACE)");
		noise = 3;
		break;
      case 'g':	timetrace_ps_interval = atoi(optarg);
		break;
      }
  }
  if(randseed==0)randseed=time(NULL);
}

// ************************************************************
// **   input and output for configuration files
// ************************************************************

void read_popdens(int importparameters) {
  int i;
  int icount,dcount;
  int *ival;
  double *dval;
  int c_space,c_space0;
  double c_dx;
  FILE *fp;
  
  fp = fopen(c_infile,"rb");
  if(fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    fread(&c_dx,sizeof(double),1,fp);
    fread(&c_space,sizeof(int),1,fp);
    fread(&c_space0,sizeof(int),1,fp);
    
    if(importparameters) {
      space = c_space;
      space0 = c_space0;
      dx = c_dx;
    }
        
    if(icount>0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    
    if(dcount>0) {
      dval = (double*)malloc(dcount*sizeof(int));
      fread(dval,sizeof(double),dcount,fp);
      free(dval);
    }
    
    if(space != c_space)print_error("lattice does not match! (read c)");
    nn = (double*)malloc(space*sizeof(double));
    fread(nn,sizeof(double),space,fp);
    for(i=0;i<space;i++) {
      nn[i] *= dx;
    }
    
    fclose(fp);
  }else{
    print_error("could not open c-infile");
  }
}



void write_popdens() {
  int icount=0,dcount=3;
  int i;
  FILE *fpc;
  
  fpc=fopen(c_outfile,"wb");
  if(fpc != NULL) {
    fwrite(&icount,sizeof(int),1,fpc);
    fwrite(&dcount,sizeof(int),1,fpc);
    fwrite(&dx,sizeof(double),1,fpc);
    fwrite(&space,sizeof(int),1,fpc);
    fwrite(&space0,sizeof(int),1,fpc);
    fwrite(&mutationrate,sizeof(double),1,fpc);
    fwrite(&wavespeed,sizeof(double),1,fpc);
    fwrite(&current_mean_fitness,sizeof(double),1,fpc);
    for(i=0;i<space;i++)nn[i] /= dx;
    fwrite(nn,sizeof(double),space,fpc);
    fclose(fpc);
  }else{
    print_error("could not open c-outfile");
  }
}





void read_constraint(int importparameters) {
  int icount,dcount;
  int *ival;
  double *dval;
  FILE *fp;
  
  fp = fopen(u_infile,"rb");
  if(fp!=NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    fread(&u_dx,sizeof(double),1,fp);
    fread(&u_space,sizeof(int),1,fp);
    fread(&u_space0,sizeof(int),1,fp);
    
    if(icount>0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    
    if(dcount>=2) {
      dval = (double*)malloc(dcount*sizeof(int));
      fread(dval,sizeof(double),dcount,fp);
      mutationrate = dval[0];
      wavespeed = dval[1];
      free(dval);
    }else{
      print_error("not enough values in constraint file! need at least 2 double parameters: mutationrate, wavespeed!");
    }
    
    if(importparameters) {
      space = u_space/dens_ustar_latticeratio;
      space0 = u_space0/dens_ustar_latticeratio;
      dx = u_dx*dens_ustar_latticeratio;
    }/*else{
      if((space*dens_ustar_latticeratio != u_space)||(space0*dens_ustar_latticeratio != u_space0))print_error("lattice does not match! u");
    }*/
    
    u_read = (double*)malloc(u_space*sizeof(double));
    fread(u_read,sizeof(double),u_space,fp);
    
    fclose(fp);
  }else{
    print_error("could not open ufile");
  }
    
  
  u = (double*)malloc(space*sizeof(double));

  
}



void flat_constraint(double size) {
  int i;
  u=(double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++)u[i] = 1./size;
}


void initialize_with_gaussian_popdens() {
  int i;
  double startvariance = wavespeed;
  if(noise==1) {
    // from Desai,Fisher (2007)
    startvariance = 2*dx*dx*log(populationsize*dx)/(log(dx/mutationrate)*log(dx/mutationrate));
  }
  for(i=0;i<space;i++) {
    nn[i] = exp(-(i-space0)*(i-space0)*dx*dx/(2.*startvariance));
  }
}

// ************************************************************
// **   timetrace
// ************************************************************

void read_timetrace() {
  int i;
  timetrace_length = (int)(maxSteps/timetrace_ps_interval)+1;
  FILE *fp;
  timetrace_popsize = (double*)malloc(timetrace_length*sizeof(double));
  fp = fopen(timetrace_filename,"r");
  if(fp != NULL) {
    for(i=0;i<timetrace_length;i++) {
      fscanf(fp,"%lf",&timetrace_popsize[i]);
//       printf("%.10e\n",timetrace_popsize[i]);
    }
  }else{
    print_error("could not open timetrace file");
  }
  fclose(fp);
}


double get_timetrace_popsize(int step) {
  int i;
  if(step%timetrace_ps_interval == 0) {
    return timetrace_popsize[step/timetrace_ps_interval];
  }else{
    int lasttimepoint = (int)floor(step/timetrace_ps_interval);
    int nexttimepoint = lasttimepoint + 1;
    double fractimepoint = (1.*(step%timetrace_ps_interval))/(1.*timetrace_ps_interval);
    return fractimepoint*timetrace_popsize[nexttimepoint] + (1-fractimepoint)*timetrace_popsize[lasttimepoint];  
  }
}




// ************************************************************
// **   screen output
// ************************************************************

 
void print_populationdensity(int timestep) {
  int i;
  double corr = allshifts*dx;
  if(correctformeanfitness)corr += current_mean_fitness;
  for(i=0;i<space;i++) {
    fprintf(stderr,"%lf %14.10lf %20.10e %20.10e\n",timestep*epsilon,(i-space0)*dx+corr,nn[i],u[i]);
    if(printhistotype)fprintf(stderr,"%lf %14.10lf %20.10e %20.10e\n",timestep*epsilon,(i-space0+1)*dx+corr,nn[i],u[i]);
  }
  fprintf(stderr,"\n");
}




// ************************************************************
// **   initialization
// ************************************************************

int initialize() {
  int i;
  
  if(read_from_file) {
    if(noise == 0) {
        read_popdens(1);
        flat_constraint(1.);
    }
    if(noise == 1) {
        read_popdens(1);
        flat_constraint(populationsize);
    }
    if(noise == 2) {
        read_constraint(1);
        read_popdens(0);
    }
    if(noise == 3) {
        read_popdens(1);
        flat_constraint(1.);
        read_timetrace();
    }
  }else{
    if(noise == 0)flat_constraint(1.);
    if(noise == 1)flat_constraint(populationsize);
    if(noise == 2)read_constraint(1);
    if(noise == 3) {
        flat_constraint(1.);
        read_timetrace();
    }
    nn = (double*)calloc(space,sizeof(double));
    initialize_with_gaussian_popdens();
  }
  tmp = (double*)malloc(space*sizeof(double));

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);

  twoepssqrt = sqrt(2.*epsilon);
  speedprefactor = wavespeed/dx*epsilon;
  
  x = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++)x[i] = (i-space0)*dx;
  mutation_jump = mutationrate;
  
  if(quiet<2) {
    printf("###############################################################################\n");
    printf("# stochastic simulation of adapting population with staircase mutation kernel #\n");
    printf("###############################################################################\n");
    printf("#    mutationrate     = %e\n",mutationrate);
    if(noise==0) {
      printf("#    constraint       = deterministic\n");
    }
    if(noise==1) {
      printf("#    constraint       = fixedN\n");
      printf("#    populationsize   = %e\n",populationsize);
    }
    if(noise==2) {
      printf("#    constraint       = ustar\n");
      printf("#    ufile            = %s\n",u_infile);
      printf("#    wavespeed        = %e\n",wavespeed);
    }
    if(noise==3) {
      printf("#    constraint       = fluctN\n");
      printf("#    N trajectory file= %s\n",timetrace_filename);
    }    
    printf("#    (lattice) space  = %d\n",space);
    printf("#    (lattice) space0 = %d\n",space0);
    printf("#    (lattice) dx     = %e\n",dx);
    printf("#    randseed         = %d\n",randseed);
  }

}


// ************************************************************
// **   average popdens
// ************************************************************

void init_averagepopdens() {
  int i;
  if(averagepopdens_havefile == 1) {
    averagepopdens_dens = (double*)calloc(space,sizeof(double));
    averagepopdens_densnose = (double*)calloc(space,sizeof(double));
  }
  
  if(averagepopdens_havecovarfile == 1) {
    averagepopdens_covar = (double**)malloc(space*sizeof(double*));
    averagepopdens_covarnose = (double**)malloc(space*sizeof(double*));
    for(i=0;i<space;i++) {
      averagepopdens_covar[i] = (double*)calloc(space,sizeof(double));
      averagepopdens_covarnose[i] = (double*)calloc(space,sizeof(double));
    }
  }
}

void update_averagepopdens() {
  int i,j,idx;
  int havenose = 0, nose_offset;

  for(j=space-1;j>=0;j--) {
    if(havenose == 0) {
      if(nn[j] > 1) {
	havenose = 1;
	nose_offset = j;
      }
    }else{
      averagepopdens_densnose[nose_offset-j] += nn[j];
    }
    averagepopdens_dens[j] += nn[j];
  }
  averagepopdens_count += 1.;
}


void update_averagecovar() {
  int i,j;
  int havenose = 0,nose_offset;
  for(j=space-1;j>=0;j--) {
    if((havenose == 0)  && (nn[j] > 1)) {
      havenose = 1;
      nose_offset = j;
      break;
    }
  }
    
  for(i=space-1;i>=0;i--) {
    for(j=space-1;j>=0;j--) {
      if((i<=nose_offset) && (j<=nose_offset))averagepopdens_covarnose[nose_offset-i][nose_offset-j] += nn[i]*nn[j];
      averagepopdens_covar[i][j] += nn[i]*nn[j];
    }
  }
  averagepopdens_covarcount += 1.;
}


void write_averagepopdens() {
  int i;
  FILE *fp;
  double norm = 1./(averagepopdens_count*dx);
  
  fp = fopen(averagepopdens_outputfile,"w");
  
  for(i=0;i<space;i++) {
    fprintf(fp,"%lf %e %e\n",(i-space0)*dx,averagepopdens_dens[i]*norm,averagepopdens_densnose[i]*norm);
  }
  
  fclose(fp);
  
  free(averagepopdens_dens);
  free(averagepopdens_densnose);
}

void write_averagecovar() {
  int i,j;
  FILE *fp;
  double norm = 1./(averagepopdens_covarcount*dx*dx);
  fp = fopen(averagepopdens_covarfile,"w");
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      fprintf(fp,"%e %e %.10e %.10e\n",(i-space0)*dx,(j-space0)*dx,averagepopdens_covar[i][j]*norm,averagepopdens_covarnose[i][j]*norm);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  for(i=0;i<space;i++) {
    free(averagepopdens_covar[i]);
    free(averagepopdens_covarnose[i]);
  }
  free(averagepopdens_covar);
  free(averagepopdens_covarnose[i]);
}

  
// ************************************************************
// **   main algorithm
// ************************************************************


void update_u(int timestep) {
   // bug if shiftthreshold -T >1, use shiftthreshold only with fixedN
  int i,j;
  int baseshift;
  double fracshift;
  
  baseshift = (int)floor(current_mean_fitness/u_dx);
  fracshift = current_mean_fitness/u_dx - 1.*baseshift;
  i = 0;
  for(j=0;j<u_space;j++) {
    if ( j + baseshift + fracshift > (i+.5)*dens_ustar_latticeratio) {
      u[i] += (1.-fracshift)*u_read[j];
      u[i] /= (1.*dens_ustar_latticeratio);
      i++;
      if(i>=space)break;
      u[i] = fracshift*u_read[j];
    }else{
      u[i] += u_read[j];
    }
  }
}


void update_mean_fit() {
  int i;
  current_mean_fitness = 0.;
  populationsize = 0.;
  for(i=0;i<space;i++) {
    current_mean_fitness += (i-space0)*nn[i];
    populationsize += nn[i];
  }
  current_mean_fitness *= dx/populationsize;
}
  


void shift_population_backward(int step) {
  int i;
  for(i=0;i<space-step;i++)     nn[i] = nn[i+step];
  for(i=space-step;i<space;i++) nn[i] = 0.;
}


void shift_population_forward(int step) {
  int i;
  for(i=space-1;i>-step;i--) nn[i] = nn[i+step];
  for(i=-step;i>=0;i--)      nn[i] = 0.;
}


void shift_population(int timestep) {
  int shift = (int)floor(current_mean_fitness/dx);
  if(shift >= shiftthreshold) {shift_population_backward(shift);}
  if(shift <= -shiftthreshold) {shift_population_forward(shift);}
  allshifts += shift;
  current_mean_fitness -= shift*dx;
}




void reproduce(int timestep) {
  int i,j;
  double tmpn;
  
  
  if(noise == 2) {
    current_mean_fitness += epsilon*(wavespeed + dv);
  }else{
    update_mean_fit();
  }
  
  if((current_mean_fitness/dx >= shiftthreshold) || (current_mean_fitness/dx <= -shiftthreshold)) {
    shift_population(timestep);
  }
      
  tmp[0] = 0;
  tmp[space-1] = 0;

  for(i=1;i<space-1;i++) {
    tmp[i] = nn[i]*(1.+(x[i]-current_mean_fitness)*epsilon);
    tmp[i] += epsilon*mutation_jump*(nn[i-1] - nn[i]);
    if(tmp[i]<0) {
      tmp[i] = 0;
    }else if(noise > 0) {
      if(tmp[i] < 1e9) { // Poisson-RNG breaks down for parameters > 1e9. see GSL doc.
			 // use smaller bins if this occurs too often or population size too large
			 // don't add noise if occupancy too large, it's anyway almost deterministic then
	nn[i] = tmp[i];
	tmp[i] += twoepssqrt*(gsl_ran_poisson(rg,nn[i])-nn[i]);
      }
    }
  }
  
  memcpy(nn,tmp,space*sizeof(double));
  
}


double populationconstraint(int timestep) {
  int i;
  double constraint = 0., inv = 1., sumxg = 0.0, sumg = 0.0;
  double tmp;
  
  if(noise == 2) {
    update_u(timestep);
  }
  
  popdens_0thmom = 0;
  popdens_1stmom = 0;
  popdens_2ndmom = 0;
  ucp = 0.;
  
  for(i=0;i<space;i++) {
    if(i>0)      ucp -= nn[i-1]*u[i];
    if(i<space-1)ucp += nn[i+1]*u[i];
    constraint += nn[i]*u[i];
    sumxg += (x[i] - current_mean_fitness)*nn[i]*u[i];
    popdens_0thmom += nn[i];
    popdens_1stmom += (i-space0)*nn[i];
    popdens_2ndmom += (i-space0)*(i-space0)*nn[i];
  }
  ucp /= 2.*dx;
  popdens_1stmom *= dx;
  popdens_2ndmom *= dx*dx;
  xg = sumxg/constraint;
  populationsize = popdens_0thmom;
  populationvariance = popdens_2ndmom/popdens_0thmom - popdens_1stmom*popdens_1stmom/(popdens_0thmom*popdens_0thmom);
  
  if(constraint>0)inv = 1./constraint;
  if(noise==3) {
      tmp = get_timetrace_popsize(timestep);
      inv *= tmp;
  }
  if(windcontrolled) {
    dv = -(1.-inv)/(ucp*epsilon);
  }else{
    for(i=0;i<space;i++) nn[i] *= inv;
  }
  populationsize *= inv;
  return inv;
}






// ************************************************************
// **   cleanup
// ************************************************************


  
void cleanup() {
  free(nn);
  free(tmp);
  free(x);
  free(mutation_inflow);
  free(u);
  if(noise == 2) {
    free(u_read);
  }
}



// ************************************************************
// **   main
// ************************************************************


int main(int argn, char *argv[]) {
  int i;
  double v;
  double lambda;
  
  parsecomamndline(argn,argv);
  initialize();
  populationconstraint(0);
  populationconstraint(0);
  if(quiet==0)print_populationdensity(0);
  if(averagepopdens) {
    init_averagepopdens();
  }
//   if (quiet<2) fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",0,0.,populationvariance,populationsize,1,xg,ucp,dv);
  if (quiet<2) fprintf(stdout,"%10.3lf %20.10e %.10e %.10e\n",0,0.,populationvariance,populationsize);
  for(i=1;i<=maxSteps;i++) {
    reproduce(i);
    lambda = populationconstraint(i);
    
    if(i%outputstep == 0) {
//       if (quiet<2)fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",i*epsilon,allshifts*dx + popdens_1stmom/popdens_0thmom,populationvariance,populationsize,lambda,xg,ucp,dv);
      if (quiet<2)fprintf(stdout,"%10.3lf %20.10e %.10e %.10e\n",i*epsilon,allshifts*dx + popdens_1stmom/popdens_0thmom,populationvariance,populationsize);
      if (quiet==0) print_populationdensity(i);
      if (averagepopdens_havefile) update_averagepopdens();
      if (averagepopdens_havecovarfile) update_averagecovar();
    }
  }
  
  if(write_to_file)write_popdens();
  if(averagepopdens_havefile)write_averagepopdens();
  if(averagepopdens_havecovarfile)write_averagecovar();
  cleanup();
  return 0;
}


/* #################################################
 * # stochastic simulation of adapting populations #
 * # with varying constraint type                  #
 * # and exponential mutation kernel               #
 * #################################################
 * 
 * 2012-2015, Lukas Geyrhofer, lukas.geyrhofer@ds.mpg.de
 * 
 * 
 * #################################################
 * 
 * simplest usage (population size 1e7):
 * ./travelingwavepeak_exp -N 1e7 1> out.txt 2> conf.txt
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

#include <mcheck.h>


int space = 300;
int space0 = 100;
double dx = 1e-2;
int maxSteps = 1000;

int outputstep = 100;
int quiet = 0;

int noise = 0;		// noise can assume 4 different values:
			// value in code	cmdline	description
			// 0 (default)			populationsize is constrained to 1, no noise is added
			//				constraint could be any value
			// 1			N	fixedN, set by parameter -N POPSIZE
			// 2			u U	ustar, set by parameter -u UFILENAME
			//				if ustar-profile is in finer grid than simulation, use parameter -U LATTICERATIO
			// 3			n g	timetrace, population size is constrained using a (single column) textfile,
			//				where the entries are the current population size
			//				use also parameter -g STEPS to indicate the number of timesteps (steps! not time in generations)
			//				between two values in this textfile
			// 4			y U	ustar, use file from -y UFILENAME
			//				no stochastic dynamics, noise is not added, only the constraint is used
			//				if ustar-profile is in finer grid than simulation, use parameter -U LATTICERATIO
double populationsize = 1.;
double populationvariance;
double xg,varg,xcmax,xnose;
int correctformeanfitness =0;
int printhistotype = 1;
int printfixationprob = 0;

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
double mutation_sigma = 1e-2;


gsl_rng* rg; /* gsl, global generator */
const gsl_rng_type* T;
unsigned long int randseed = 0;



int averagepopdens = 0;
int averagepopdens_center = 1;
int averagepopdens_resolution = 1;
int averagepopdens_havefile = 0;
double *averagepopdens_dens;
char averagepopdens_outputfile[128];
double averagepopdens_count = 0.;
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
  while((c = getopt(argn,argv,"s:z:d:S:e:D:O:qQi:o:R:N:CM:u:U:T:H:h:Pn:g:y:F")) != -1) {
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
		if(noise > 0)print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME or -n FILENAME_TIMETRACE)");
		noise = 2;
		break;
      case 'y':	strcpy(u_infile,optarg);
		if(noise > 0)print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME or -n FILENAME_TIMETRACE or -y FILENAME)");
		noise = 4;
		break;
      case 'U':	if((noise == 2)||(noise == 4)) {
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
		if(noise > 0)print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME or -n FILENAME_TIMETRACE)");
		noise = 1;
		break;
      case 'C':	correctformeanfitness = 1;
		break;
      case 'M':	mutation_sigma = atof(optarg);
		break;
      case 'T':	shiftthreshold = atoi(optarg);
		break;
      case 'h':	strcpy(averagepopdens_outputfile,optarg);
		averagepopdens_havefile = 1;
		averagepopdens = 1;
		break;
      case 'H':	averagepopdens = 1;
		averagepopdens_resolution = atoi(optarg);
		if(averagepopdens_resolution < 0) {
		  averagepopdens_resolution *= -1;
		  averagepopdens_center = 0;
		}
		break;
      case 'P':	printhistotype = 0;
		break;
      case 'F': printfixationprob = 1;
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
    
    if(space != c_space)print_error("lattice does not match!");
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
  int icount=0,dcount=0;
  int i;
  FILE *fpc;
  
  fpc=fopen(c_outfile,"wb");
  if(fpc != NULL) {
    fwrite(&icount,sizeof(int),1,fpc);
    fwrite(&dcount,sizeof(int),1,fpc);
    fwrite(&dx,sizeof(double),1,fpc);
    fwrite(&space,sizeof(int),1,fpc);
    fwrite(&space0,sizeof(int),1,fpc);
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
    
    if(dcount>=3) {
      dval = (double*)malloc(dcount*sizeof(int));
      fread(dval,sizeof(double),dcount,fp);
      mutationrate = dval[0];
      wavespeed = dval[1];
      mutation_sigma = dval[2];
      free(dval);
    }else{
      print_error("not enough values in constraint file! need at least 3 double parameters: mutationrate, wavespeed, mutationsigma!");
    }
    
    if(importparameters) {
      space = u_space/dens_ustar_latticeratio;
      space0 = u_space0/dens_ustar_latticeratio;
      dx = u_dx*dens_ustar_latticeratio;
    }else{
      if((space*dens_ustar_latticeratio != u_space)||(space0*dens_ustar_latticeratio != u_space0))print_error("lattice does not match! u");
    }
    
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
  double startvariance;
  if((noise == 2)||(noise == 4)) {
    startvariance = wavespeed;
  }else{
    // Good et al., PNAS (2012)
    startvariance = mutation_sigma*mutation_sigma*2*log(mutation_sigma*populationsize)/(log(mutation_sigma/mutationrate)*log(mutation_sigma/mutationrate));
    // behaves weird for extremal parameter values.
    // catch those cases and set on standard value: density takes half of the box
    if(startvariance < 0) {
      startvariance = (0.25*space*space*dx*dx);
    }
  }
  for(i=0;i<space;i++) {
    nn[i] = exp(-(i-space0)*(i-space0)*dx*dx/(2.*startvariance));
  }
}


void print_u(int time) {
  int i;
  for(i=0;i<space;i++) {
    fprintf(stderr,"%10.6lf %10.6lf %16.10e\n",time*epsilon,(i-space0+allshifts)*dx,u[i]);
  }
  fprintf(stderr,"\n");
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
    fprintf(stderr,"%lf %14.10lf %20.10e",timestep*epsilon,(i-space0)*dx+corr,nn[i]);
    if(printfixationprob)fprintf(stderr," %20.10e",u[i]);
    fprintf(stderr,"\n");
    
    if(printhistotype) {
      fprintf(stderr,"%lf %14.10lf %20.10e",timestep*epsilon,(i-space0+1)*dx+corr,nn[i]);
      if(printfixationprob)fprintf(stderr," %20.10e",u[i]);
      fprintf(stderr,"\n");
    }
  }
  fprintf(stderr,"\n");
}





// ************************************************************
// **   average popdens
// ************************************************************

void init_averagepopdens() {
  
  averagepopdens_dx = dx/(1.*averagepopdens_resolution);
  averagepopdens_space = space*averagepopdens_resolution;
  averagepopdens_space0 = space0*averagepopdens_resolution;
  
  averagepopdens_dens = (double*)calloc(averagepopdens_space,sizeof(double));
  
  averagepopdens_lower = (int)(-averagepopdens_resolution/2);
  averagepopdens_higher = (int)(averagepopdens_resolution/2);
  if(averagepopdens_higher - averagepopdens_lower < averagepopdens_resolution)averagepopdens_higher++;
  
}


void update_averagepopdens() {
  int i,j,idx;
  int offset=0;
  
  if(averagepopdens_center) {
    offset = (int)(current_mean_fitness/dx*averagepopdens_resolution);
  }
  
  for(i=0;i<space;i++) {
    for(j=averagepopdens_lower;j<averagepopdens_higher;j++) {
      idx = i*averagepopdens_resolution+j-offset;
      if((averagepopdens_space > idx) &&( idx >= 0))
	averagepopdens_dens[idx] += nn[i];
    }
  }
  
  averagepopdens_count += 1.;
  
}


void write_averagepopdens() {
  int i;
  FILE *fp;
  
  if(averagepopdens_havefile == 1) {
    fp = fopen(averagepopdens_outputfile,"w");
  }else{
    fp = stdout;
  }
  
  for(i=0;i<averagepopdens_space;i++) {
    fprintf(fp,"%lf %e\n",(i-averagepopdens_space0)*averagepopdens_dx,averagepopdens_dens[i]/averagepopdens_count);
  }
  
  if(averagepopdens_havefile == 1) {
    fclose(fp);
  }
  
  free(averagepopdens_dens);
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
  if(shift <= -shiftthreshold) {shift_population_forward(-shift);}
  allshifts += shift;
  current_mean_fitness -= shift*dx;
}




void reproduce(int timestep) {
  int i,j;
  double tmpn;
  long int poiss;
  
  if((noise == 2)||(noise==4)) {
    current_mean_fitness = timestep*wavespeed*epsilon - allshifts*dx;
  }	// else:	current_mean_fitness was already updated at end of last timestep
	//		with function "populationconstraint" below and is still at the correct value

  if((current_mean_fitness/dx >= shiftthreshold) || (current_mean_fitness/dx <= -shiftthreshold)) {
    shift_population(timestep);
  }
      
  tmp[0] = 0.;
  tmp[space-1] = 0.;
  for(i=1;i<space-1;i++) {
    tmpn = nn[i]*(1.+(x[i]-current_mean_fitness)*epsilon-mutation_outflow);
    for(j = 1;j<i;j++) tmpn += mutation_inflow[j]*nn[i-j];
    if(tmpn<0) {
      tmpn = 0.;
    }else if((noise > 0)&&(noise != 4)) {
      if(tmpn < 1.e9) {	//	Poisson-RNG breaks down for parameters > 1e9. see GSL doc.
			//	use smaller bins if this occurs too often or population size too large
	tmp[i] = tmpn + twoepssqrt*(gsl_ran_poisson(rg,tmpn)-tmpn);
      }else{		//	neglect noise term if occupancy is too large,
	tmp[i] = tmpn;	//	dynamics can be assumed to be almost deterministic then
      }
    }else{
      tmp[i] = tmpn;
    }
  }
  
  memcpy(&nn[1],&tmp[1],(space-2)*sizeof(double));
  
}


double populationconstraint(int timestep) {
  int i;
  double constraint = 0., inv, sumxg = 0., sumxxg = 0.;

  double popdens_0thmom = 0.;
  double popdens_1stmom = 0.;
  double popdens_2ndmom = 0.;
  double cmax = 0.;
  int havenose = 0;
  if((noise == 2)||(noise == 4)) {
    update_u(timestep);
  }
  
  for(i=space-1;i>=0;i--) {
    if (nn[i] > cmax) {
      cmax=nn[i];
      xcmax = x[i];
    }
    if((havenose==0)&&(nn[i]>1)) {
      havenose = 1;
      xnose = x[i];
    }
    constraint += nn[i]*u[i];
    sumxg += x[i] * nn[i] * u[i];
    sumxxg += x[i] * x[i] * nn[i] * u[i];
    popdens_0thmom += nn[i];
    popdens_1stmom += x[i]*nn[i];
    popdens_2ndmom += x[i]*x[i]*nn[i];
  }
  
  current_mean_fitness = popdens_1stmom/popdens_0thmom;
  populationvariance = popdens_2ndmom/popdens_0thmom - current_mean_fitness*current_mean_fitness;
  populationsize = popdens_0thmom;
  xg = sumxg/constraint;
  varg = sumxxg/constraint - xg*xg;
  
  inv = 1./constraint;
  if(noise==3) {
    inv *= get_timetrace_popsize(timestep);
  }
  for(i=0;i<space;i++) nn[i] *= inv;
  populationsize *= inv;
  
  return inv;
}




// ************************************************************
// **   initialization
// ************************************************************

int initialize() {
  int i;
  
  if(read_from_file) {
    read_popdens(1);
    if(noise == 0)flat_constraint(1.);
    if(noise == 1)flat_constraint(populationsize);
    if((noise == 2)||(noise == 4)) {
      read_constraint(0);
      update_u(0);
    }
    if(noise == 3) {
      flat_constraint(1.);
      read_timetrace();
    }
  }else{
    if(noise == 0)flat_constraint(1.);
    if(noise == 1)flat_constraint(populationsize);
    if((noise == 2) || (noise == 4)) {
      read_constraint(1);
      update_u(0);
    }
    if(noise == 3) {
      flat_constraint(1.);
      read_timetrace();
    }
    nn = (double*)malloc(space*sizeof(double));
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
  mutation_inflow = (double*)malloc(space*sizeof(double));
  mutation_outflow = 0.;
  mutation_inflow[0] = 0.;
  for(i=1;i<space;i++) {
    x[i] = (i-space0)*dx;
    mutation_inflow[i] = epsilon*mutationrate*exp(-i*dx/mutation_sigma)*(1.-exp(-dx/mutation_sigma));
    mutation_outflow += mutation_inflow[i];
  }
  
  if(quiet<2) {
    printf("#################################################################################\n");
    printf("# stochastic simulation of adapting population with exponential mutation kernel #\n");
    printf("#################################################################################\n");
    printf("#    mutationrate     = %e\n",mutationrate);
    printf("#    mutation_sigma   = %e\n",mutation_sigma);
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
      printf("#    constraint       = fixedN, timetrace\n");
      printf("#    timetracefile    = %s\n",timetrace_filename);
    }
    if(noise==4) {
      printf("#    constraint       = ustar, non-stochastic\n");
      printf("#    ufile            = %s\n",u_infile);
      printf("#    wavespeed        = %e\n",wavespeed);
    }      
    printf("#    (lattice) space  = %d\n",space);
    printf("#    (lattice) space0 = %d\n",space0);
    printf("#    (lattice) dx     = %e\n",dx);
    printf("#    randseed         = %d\n",randseed);
  }

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
  if((noise == 2)||(noise == 4)) {
    free(u_read);
  }
  if(noise == 3) {
    free(timetrace_popsize);
  }
  gsl_rng_free(rg);
}



// ************************************************************
// **   main
// ************************************************************


int main(int argn, char *argv[]) {
  int i,j;
  double v,lambda = 1.;
  
  
  parsecomamndline(argn,argv);
  initialize();
  lambda = populationconstraint(0);
  if (quiet<2) fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",0.0,current_mean_fitness,populationvariance,populationsize,lambda,xg,varg,allshifts*dx+xcmax,allshifts*dx+xnose);
  if(averagepopdens) {
    init_averagepopdens();
  }
  for(i=1;i<=maxSteps;i++) {
    reproduce(i);
    lambda = populationconstraint(i);

    if(i%outputstep == 0) {
      if (quiet<2)fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",i*epsilon,allshifts*dx + current_mean_fitness,populationvariance,populationsize,lambda,xg,varg,allshifts*dx+xcmax,allshifts*dx+xnose);
      if (quiet==0) print_populationdensity(i);
      if (averagepopdens) update_averagepopdens();
    }
  }
  
  if(write_to_file)write_popdens();
  if(averagepopdens)write_averagepopdens();
  cleanup();
  return 0;
}


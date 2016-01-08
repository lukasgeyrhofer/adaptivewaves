/* #################################################
 * # stochastic simulation of adapting populations #
 * # with variing constraint type                  #
 * # and exponential mutation kernel               #
 * #################################################
 * 
 * 2012-2014, Lukas Geyrhofer
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


int space = 300;
int space0 = 100;
double dx = 1e-2;
int maxSteps = 1000;

int outputstep = 100;
int quiet = 0;

int noise = 0;
double populationsize = 1.;
double populationvariance;
int correctformeanfitness =0;
int printhistotype = 0;

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

double mutationrate = 1e-9;
double mutation_2dx;


const gsl_rng* rg; /* gsl, global generator */
const gsl_rng_type* T;
unsigned long int randseed = 0;


double popdens_0thmom;
double popdens_1stmom;
double popdens_2ndmom;



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


int fixationextinction_events = 0;
int subpop_labeltype = 2;
double subpop_labelparameter = 1e-4;
double subpop_expected_fixationprobability = .5;
double subpop_predicted_fixationprobability;
double *ww,*vv;
double *subpop_start_ww, *subpop_start_vv;
double *tmpw,*tmpv;
double subpop_final_threshold = 1e-5; // if (current_fixation_prob > 1 - 1e-5) or (current_fixation_prob < 1e-5), then restart...
double subpop_popsize;
int subpop_count_extinctions = 0;
int subpop_count_fixations = 0;
double current_fixation_prob;
double subpop_starting_mean_fitness;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

// ************************************************************
// **   parameters
// ************************************************************


void parsecomamndline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn,argv,"s:z:d:S:e:D:O:qQi:o:R:N:Cu:U:T:H:h:PE:l:L:")) != -1) {
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
		  print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME)");
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
		  print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME)");
		noise = 1;
		break;
      case 'C':	correctformeanfitness = 1;
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
      case 'P':	printhistotype = 1;
		break;
      case 'E':	subpop_expected_fixationprobability = atof(optarg);
		break;
      case 'l': subpop_labeltype = atoi(optarg);
                break;
      case 'L': subpop_labelparameter = atof(optarg);
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
      if(dcount >= 2) {
          mutationrate = dval[0];
          wavespeed = dval[1];
          if(dcount >= 3)current_mean_fitness = dval[2];
      }
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
  if(noise == 2) {
    startvariance = wavespeed;
  }else{
    // from Neher, Hallatschek, PNAS (2013)
    startvariance = exp(0.33333*log(24.*mutationrate*mutationrate*log(populationsize*exp(0.3333*log(mutationrate))) ));
  }
  for(i=0;i<space;i++) {
    nn[i] = exp(-(i-space0)*(i-space0)*dx*dx/(2.*startvariance));
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
    fprintf(stderr,"%lf %14.10lf %20.10e %20.10e %20.10e\n",timestep*epsilon,(i-space0)*dx+corr,ww[i]+vv[i],ww[i],u[i]);
    if(printhistotype)fprintf(stderr,"%lf %14.10lf %20.10e %20.10e %20.10e\n",timestep*epsilon,(i-space0+1)*dx+corr,ww[i]+vv[i],ww[i],u[i]);
  }
  fprintf(stderr,"\n");
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
    if(noise == 2)read_constraint(0);
  }else{
    if(noise == 0)flat_constraint(1.);
    if(noise == 1)flat_constraint(populationsize);
    if(noise == 2)read_constraint(1);
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
  for(i=0;i<space;i++)x[i]=(i-space0)*dx;
  mutation_2dx = epsilon*mutationrate/(dx*dx);
  
  
  if(quiet<2) {
    printf("#################################################################################\n");
    printf("# stochastic simulation of adapting population with exponential mutation kernel #\n");
    printf("#################################################################################\n");
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
      printf("#    wavespeed        = %e\n",wavespeed);
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


void update_mean_fit() {
  int i;
  current_mean_fitness = 0.;
  populationsize = 0.;
  for(i=0;i<space;i++) {
    current_mean_fitness += (i-space0)*(ww[i]+vv[i]);
    populationsize += ww[i] + vv[i];
  }
  current_mean_fitness *= dx/populationsize;
}
  


void shift_population_backward(int step) {
  int i;
  for(i=0;i<space-step;i++){
    ww[i] = ww[i+step];
    vv[i] = vv[i+step];
  }
  for(i=space-step;i<space;i++) {
    ww[i] = 0.;
    vv[i] = 0.;
  }
}


void shift_population_forward(int step) {
  int i;
  for(i=space-1;i>step;i--) {
    ww[i] = ww[i-step];
    vv[i] = vv[i-step];
  }
  for(i=step;i>=0;i--) {
    ww[i] = 0.;
    vv[i] = 0.;
  }
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
  
  if(noise < 2) {
    update_mean_fit();
  }else if(noise == 2) {
    current_mean_fitness += epsilon*wavespeed;
  }
  
  if((current_mean_fitness > shiftthreshold*dx) || (current_mean_fitness < -shiftthreshold*dx)) {
    shift_population(timestep);
  }
      
  tmpw[0] = 0;
  tmpw[space-1] = 0;
  
  tmpv[0] = 0;
  tmpv[space-1] = 0;

  for(i=1;i<space-1;i++) {
    tmpw[i] = ww[i]*(1.+(x[i]-current_mean_fitness)*epsilon);
    tmpv[i] = vv[i]*(1.+(x[i]-current_mean_fitness)*epsilon);
    tmpw[i] += mutation_2dx*(ww[i-1]-2.*ww[i]+ww[i+1]);
    tmpv[i] += mutation_2dx*(vv[i-1]-2.*vv[i]+vv[i+1]);
    if(tmpw[i]<0) {
      tmpw[i] = 0;
    }else if(noise > 0) {
      if(tmpw[i] < 1e9) { // Poisson-RNG breaks down for parameters > 1e9. see GSL doc.
			 // use smaller bins if this occurs too often or population size too large
            ww[i] = tmpw[i];
            tmpw[i] += twoepssqrt*(gsl_ran_poisson(rg,ww[i])-ww[i]);
      }
    }
    
    if(tmpv[i]<0) {
      tmpv[i] = 0;
    }else if(noise > 0) {
      if(tmpv[i] < 1e9) { // Poisson-RNG breaks down for parameters > 1e9. see GSL doc.
			 // use smaller bins if this occurs too often or population size too large
            vv[i] = tmpv[i];
            tmpv[i] += twoepssqrt*(gsl_ran_poisson(rg,vv[i])-vv[i]);
      }
    }
    
  }
  
  memcpy(ww,tmpw,space*sizeof(double));
  memcpy(vv,tmpv,space*sizeof(double));
  
}


void populationconstraint(int timestep) {
  int i;
  double constraint = 0., inv;
  
  if(noise == 2) {
    update_u(timestep);
  }
  current_fixation_prob = 0;
  popdens_0thmom = 0;
  popdens_1stmom = 0;
  popdens_2ndmom = 0;
  subpop_popsize = 0;
  for(i=0;i<space;i++) {
    constraint += (ww[i] + vv[i])*u[i];
    current_fixation_prob += ww[i]*u[i];
    subpop_popsize += ww[i];
    popdens_0thmom += ww[i] + vv[i];
    popdens_1stmom += x[i]*(ww[i] + vv[i]);
    popdens_2ndmom += x[i]*x[i]*(ww[i] + vv[i]);
  }
  
  if(constraint>0)inv = 1./constraint;
  for(i=0;i<space;i++) {
    ww[i] *= inv;
    vv[i] *= inv;
  }
  current_fixation_prob *= inv;
  popdens_0thmom        *= inv;
  popdens_1stmom        *= inv;
  popdens_2ndmom        *= inv;
  
  populationsize = popdens_0thmom;
  populationvariance = popdens_2ndmom/popdens_0thmom - popdens_1stmom*popdens_1stmom/(popdens_0thmom*popdens_0thmom);
}



void initial_constraint() {
  int i;
  double sum=0;
  double inv;
  for(i=0;i<space;i++)sum += nn[i] * u[i];
  inv = 1./sum;
  for(i=0;i<space;i++)nn[i] *= inv;
}




// ************************************************************
// **   extinction and fixation events
// ************************************************************

double get_subpop_fraction() {
  int i;
  double nnsum = 0.,wwsum = 0.;
  for(i=0;i<space;i++) {
    nnsum += ww[i]+vv[i];
    wwsum += ww[i];
  }
  return wwsum/nnsum;
}

double get_current_fixprob() {
    int i;
    double nnsum = 0.,wwsum = 0.;
    for(i=0;i<space;i++) {
        nnsum += (ww[i] + vv[i])*u[i];
        wwsum += ww[i]*u[i];
    }
    return wwsum/nnsum;
}

double cutoff_function(double pos,double pos0,double scale) {
  return 1./(1.+exp(-(pos-pos0)/scale));
}

int init_subpopulation() {
    int i,j;
    int label_startindex = space;
    double fixprob = 0.,frac_lastbin;
    double subpop_crossover_step = dx;
    double subpop_crossoverpoint = -space0*dx;
    double subpop_last_fixationprob;
    double subpop_current_fixationprob;
    
    ww = (double*)calloc(space,sizeof(double));
    vv = (double*)calloc(space,sizeof(double));
    subpop_start_ww = (double*)malloc(space*sizeof(double));
    subpop_start_vv = (double*)malloc(space*sizeof(double));
    tmpw = (double*)malloc(space*sizeof(double));
    tmpv = (double*)malloc(space*sizeof(double));
    
    subpop_starting_mean_fitness = current_mean_fitness;
    update_u(0);
    initial_constraint();
    
    switch(subpop_labeltype) {
        case 0:	label_startindex = space;
                while (fixprob < subpop_expected_fixationprobability) {
                    label_startindex--;
                    fixprob += nn[label_startindex]*u[label_startindex];
                }

		if(label_startindex > space)print_error("subpop label threshold larger than simulationbox");
		memcpy(&ww[label_startindex],&nn[label_startindex],(space-label_startindex)*sizeof(double));
		memcpy(&vv[0],&nn[0],label_startindex*sizeof(double));
		
		frac_lastbin = (fixprob - subpop_expected_fixationprobability)/(nn[label_startindex]*u[label_startindex]);
		ww[label_startindex] = nn[label_startindex] * (1. - frac_lastbin);
		vv[label_startindex] = nn[label_startindex] * frac_lastbin;
		
		subpop_current_fixationprob = get_current_fixprob();
		break;
        case 1:	for(i=0;i<space;i++) {
		  ww[i] = subpop_expected_fixationprobability * nn[i];
		  vv[i] = (1.-subpop_expected_fixationprobability) * nn[i];
		}
		subpop_current_fixationprob = get_current_fixprob();
		break;
        case 2:	//fprintf(stderr,"# finding crossoverpoint\n");
		j = 0;
		if(subpop_labelparameter < 0)subpop_labelparameter = 10*dx; // default value
		for(i=0;i<space;i++) {
		  ww[i] = cutoff_function((i-space0)*dx,subpop_crossoverpoint,subpop_labelparameter)*nn[i];
		  vv[i] = nn[i] - ww[i];
		}
		subpop_current_fixationprob = get_current_fixprob();
		subpop_last_fixationprob = subpop_current_fixationprob;
		while(fabs(subpop_current_fixationprob - subpop_expected_fixationprobability)>1e-10) {
		  subpop_crossoverpoint += subpop_crossover_step;
		  subpop_current_fixationprob = 0;
		  for(i=0;i<space;i++) {
		    ww[i] = cutoff_function((i-space0)*dx,subpop_crossoverpoint,subpop_labelparameter)*nn[i];
		    vv[i] = nn[i] - ww[i];
		    subpop_current_fixationprob += ww[i]*u[i];
		  }
		  //fprintf(stderr,"# %d %.10e %.10e %.10e\n",j++,subpop_current_fixationprob,subpop_crossover_step,subpop_crossoverpoint);
		  if( (subpop_current_fixationprob - subpop_expected_fixationprobability)*(subpop_last_fixationprob - subpop_expected_fixationprobability) < 0)subpop_crossover_step *= -.5;
		  subpop_last_fixationprob = subpop_current_fixationprob;
		}
		break;
        default:print_error("label type not implemented");
		break;
    }
  
    memcpy(subpop_start_ww,ww,space*sizeof(double));
    memcpy(subpop_start_vv,vv,space*sizeof(double));
}


double reset_subpopulation() {
  allshifts = 0;
  current_mean_fitness = subpop_starting_mean_fitness;
  memcpy(ww,subpop_start_ww,space*sizeof(double));
  memcpy(vv,subpop_start_vv,space*sizeof(double));
}




// ************************************************************
// **   cleanup
// ************************************************************


  
void cleanup() {
  free(nn);
  free(tmp);
  free(x);
  free(u);
  if(noise == 2) {
    free(u_read);
  }
  free(ww);
  free(vv);
  free(tmpw);
  free(tmpv);
  free(subpop_start_ww);
  free(subpop_start_vv);
}



// ************************************************************
// **   main
// ************************************************************


int main(int argn, char *argv[]) {
  int i=1;
  double v;
  
  parsecomamndline(argn,argv);
  initialize();
  init_subpopulation();
  populationconstraint(0);
  
//   print_populationdensity(0);
//   exit(1);
  
  if(averagepopdens) {
    init_averagepopdens();
  }
  if (quiet<2) fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e %.10e\n",i*epsilon,0.,populationvariance,populationsize,subpop_popsize,current_fixation_prob);
  
  while(subpop_count_extinctions + subpop_count_fixations < maxSteps) {
    reproduce(i);
    populationconstraint(i);
    
    
    
    if(i%outputstep == 0) {
      if (quiet<2)fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e %.10e\n",i*epsilon,allshifts*dx + popdens_1stmom/popdens_0thmom,populationvariance,populationsize,subpop_popsize,current_fixation_prob);
      if (quiet==0) print_populationdensity(i);
      if (averagepopdens) update_averagepopdens();
    }
    if(current_fixation_prob < subpop_final_threshold) {
      reset_subpopulation();
      subpop_count_extinctions++;
      printf("# event ( %d of %d ): extinction %lf\n",subpop_count_extinctions+subpop_count_fixations,maxSteps,i*epsilon);
      i=0;
    }
    if(current_fixation_prob > 1-subpop_final_threshold) {
      reset_subpopulation();
      subpop_count_fixations++;
      printf("# event ( %d of %d ): fixation %lf\n",subpop_count_extinctions+subpop_count_fixations,maxSteps,i*epsilon);
      i=0;
    }
    i++;
  }
  if(write_to_file)write_popdens();
  if(averagepopdens)write_averagepopdens();
  cleanup();
  return 0;
}


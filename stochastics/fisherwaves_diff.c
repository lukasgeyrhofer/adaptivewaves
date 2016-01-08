/* #################################################
 * # stochastic simulation of fisherwaves          #
 * #################################################
 * 
 * 2012-2015, Lukas Geyrhofer
 * 
 * 
 * #################################################
 * 
 * simplest usage :
 * ./fisherwaves_diff -u ufile.conf 1> out.txt 2> conf.txt
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
double current_wave_center = 0.;

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

double diffusionconstant = 1e-5;
// double *mutation_inflow;
// double mutation_outflow;
double fisherwave_step = 1e-3;

double mutation_diff2dx = 0.;

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



int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

// ************************************************************
// **   parameters
// ************************************************************


void parsecomamndline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn,argv,"s:z:d:S:e:D:O:qQi:o:R:N:Cu:U:T:H:h:P")) != -1) {
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
      case 'D':	diffusionconstant = atof(optarg);
		break;
      case 'O':	outputstep = atoi(optarg);
		break;
      case 'q':	quiet = 2;
		break;
      case 'Q':	quiet = 1;
		break;
      case 'R':	randseed =atoi(optarg);
		break;
//       case 'N': populationsize = atof(optarg);
// 		if(noise > 0)
// 		  print_error("Only a single constraint-type can be used (either option -N POPSIZE or -u FILENAME)");
// 		noise = 1;
// 		break;
      case 'C':	correctformeanfitness = 1;
		break;
//       case 'M':	mutation_sigma = atof(optarg);
// 		break;
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
      diffusionconstant = dval[0];
      wavespeed = dval[1];
      fisherwave_step = dval[2];
      free(dval);
    }else{
      print_error("not enough values in constraint file! need at least 3 double parameters: diffusionconstant, wavespeed, stepsize!");
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
  double startvariance = wavespeed;
// **** holds for exp kernel: ****
//   if(noise == 2) {
//     startvariance = wavespeed;
//   }else{
//     // Good et al., PNAS (2012)
//     startvariance = mutation_sigma*mutation_sigma*2*log(mutation_sigma*populationsize)/(log(mutation_sigma/mutationrate)*log(mutation_sigma/mutationrate));
//   }
  for(i=0;i<space;i++) {
    nn[i] = exp(-(i-space0)*(i-space0)*dx*dx/(2.*startvariance));
  }
}

void initialize_popdens_with_stepfunction() {
  int i;
  for(i=0;i<space0;i++)nn[i] = 1;
  for(i=space0;i<space;i++)nn[i] = 0.;
}


// ************************************************************
// **   screen output
// ************************************************************

 
void print_populationdensity(int timestep) {
  int i;
  double corr = allshifts*dx;
  if(correctformeanfitness)corr += current_wave_center;
  for(i=0;i<space;i++) {
    fprintf(stderr,"%lf %14.10lf %20.10e\n",timestep*epsilon,(i-space0)*dx+corr,nn[i]);
    if(printhistotype)fprintf(stderr,"%lf %14.10lf %20.10e\n",timestep*epsilon,(i-space0+1)*dx+corr,nn[i]);
  }
  fprintf(stderr,"\n");
}




// ************************************************************
// **   initialization
// ************************************************************

int initialize() {
  int i;
  
  if(read_from_file) {
    if(noise == 0)flat_constraint(1.);
    if(noise == 1)flat_constraint(populationsize);
    if(noise == 2)read_constraint(1);
    read_popdens(0);
  }else{
    if(noise == 0)flat_constraint(1.);
    if(noise == 1)flat_constraint(populationsize);
    if(noise == 2)read_constraint(1);
    nn = (double*)calloc(space,sizeof(double));
    initialize_popdens_with_stepfunction();
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
  mutation_diff2dx = diffusionconstant/(dx*dx);
  
//   printf("D/(dx^2) = %e\n",mutation_diff2dx);
  
  if(quiet<2) {
    printf("#######################################################################\n");
    printf("# stochastic simulation of fisherwaves with diffusion mutation kernel #\n");
    printf("#######################################################################\n");
    printf("#    diffusionconstant = %e\n",diffusionconstant);
    if(noise==0) {
      printf("#    constraint        = deterministic\n");
    }
    if(noise==1) {
      printf("#    constraint        = fixedN\n");
      printf("#    populationsize    = %e\n",populationsize);
    }
    if(noise==2) {
      printf("#    constraint        = ustar\n");
      printf("#    ufile             = %s\n",u_infile);
      printf("#    wavespeed         = %e\n",wavespeed);
    }
    printf("#    (lattice) space   = %d\n",space);
    printf("#    (lattice) space0  = %d\n",space0);
    printf("#    (lattice) dx      = %e\n",dx);
    printf("#    randseed          = %d\n",randseed);
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
    offset = (int)(current_wave_center/dx*averagepopdens_resolution);
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
  
  if(dens_ustar_latticeratio > 1) {
    baseshift = (int)(timestep*epsilon*wavespeed/u_dx) - allshifts*dens_ustar_latticeratio;
    fracshift = timestep*epsilon*wavespeed/u_dx - 1.*(baseshift + allshifts*dens_ustar_latticeratio);
    i = (int)(baseshift/dens_ustar_latticeratio);
    j = 0;
    while( (i<space) && (j<u_space) ) {
      if((j+baseshift)%dens_ustar_latticeratio == 0) {
	u[i] = (fracshift) * u_read[j];
      }else if((j+baseshift)%dens_ustar_latticeratio == dens_ustar_latticeratio - 1) {
	u[i] += (1.-fracshift) * u_read[j];
	u[i] /= (1.*dens_ustar_latticeratio);
	i++;
      }else{
	u[i] += u_read[j];
      }
      j++;
    }
  }else{
    fracshift = timestep*epsilon*wavespeed/u_dx - 1.*allshifts;
    // ignore u[0], as it is unlikely that the population has a significant subpop there...
    for(i=1;i<space;i++)u[i] = fracshift * u_read[i-1] + (1.-fracshift)*u_read[i];
  }
}


void update_wave_center() {
  int i;
  double normalization_value = 0.;
  
  // based on the observation that nn(x)*nn(-x) is maximal at x \approx 0

  current_wave_center = 0.;
  for(i=0;i<space;i++) {
    current_wave_center += (i-space0)*nn[i]*nn[space-1-i];
    normalization_value += nn[i]*nn[space-1-i];
  }
  current_wave_center *= dx/normalization_value;
}
  


void shift_population_backward(int step) {
  int i;
  for(i=0;i<space-step;i++)     nn[i] = nn[i+step];
  for(i=space-step;i<space;i++) nn[i] = 0.;
}


void shift_population_forward(int step) {
  int i;
  for(i=space-1;i>step;i--) nn[i] = nn[i-step];
  for(i=step;i>=0;i--)      nn[i] = 0.;
}


void shift_population(int timestep) {
  int shift = (int)floor(current_wave_center/dx);
  if(shift >= shiftthreshold) {shift_population_backward(shift);}
  if(shift <= -shiftthreshold) {shift_population_forward(shift);}
  allshifts += shift;
  current_wave_center -= shift*dx;
}




void reproduce(int timestep) {
  int i,j;
  double tmpn;
  
  
  if(noise < 2) {
    update_wave_center();
  }else if(noise == 2) {
    current_wave_center = timestep*wavespeed*epsilon - allshifts*dx;
  }
  
  if((current_wave_center/dx >= shiftthreshold) || (current_wave_center/dx <= -shiftthreshold)) {
    shift_population(timestep);
  }
      
  tmp[0] = nn[0] + epsilon*mutation_diff2dx*(-nn[0]+nn[1]); //assume nn[-1] = nn[0]
  nn[0] = tmp[0];
  if(tmp[0]<0) {
    tmp[0] = 0;
  }else if(noise > 0) {
    if(tmp[0] < 1e9) { // Poisson-RNG breaks down for parameters > 1e9. see GSL doc.
			// use smaller bins if this occurs too often or population size too large
      tmp[0] += twoepssqrt*(gsl_ran_poisson(rg,nn[0])-nn[0]);
    }else{
      tmp[0] = 1e9;
// 	printf("# RESET BIN: time = %lf, n[%d] = 1e9\n",timestep*epsilon,i); 
    }
  }
  tmp[space-1] = 0;

  for(i=1;i<space-1;i++) {
    tmp[i] = nn[i]*(1.+((x[i]-current_wave_center>0)*fisherwave_step)*epsilon);
    tmp[i] += epsilon*mutation_diff2dx*(nn[i-1]-2.*nn[i]+nn[i+1]);
    nn[i] = tmp[i];
    if(tmp[i]<0) {
      tmp[i] = 0;
    }else if(noise > 0) {
      if(tmp[i] < 1e9) { // Poisson-RNG breaks down for parameters > 1e9. see GSL doc.
			 // use smaller bins if this occurs too often or population size too large
	tmp[i] += twoepssqrt*(gsl_ran_poisson(rg,nn[i])-nn[i]);
      }else{
	tmp[i] = 1e9;
// 	printf("# RESET BIN: time = %lf, n[%d] = 1e9\n",timestep*epsilon,i); 
      }
    }
  }
  
  memcpy(nn,tmp,space*sizeof(double));
  
}


double populationconstraint(int timestep) {
  int i;
  double sum = 0., inv;
  
  if(noise == 2) {
    update_u(timestep);
  }
  
  popdens_0thmom = 0;
  popdens_1stmom = 0;
  popdens_2ndmom = 0;
  
  for(i=0;i<space;i++) {
    sum += nn[i]*u[i];
    popdens_0thmom += nn[i];
    popdens_1stmom += (i-space0)*nn[i];
    popdens_2ndmom += (i-space0)*(i-space0)*nn[i];
  }
  popdens_1stmom *= dx;
  popdens_2ndmom *= dx*dx;
  
  populationsize = popdens_0thmom;
  populationvariance = popdens_2ndmom/popdens_0thmom - popdens_1stmom*popdens_1stmom/(popdens_0thmom*popdens_0thmom);
  
  inv = 1./(sum);
  for(i=0;i<space;i++) nn[i] *= inv;
  return inv;
}






// ************************************************************
// **   cleanup
// ************************************************************


  
void cleanup() {
  free(nn);
  free(tmp);
  free(x);
//   free(mutation_inflow);
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
  if(quiet==0)print_populationdensity(0);
  if(averagepopdens) {
    init_averagepopdens();
  }
  if (quiet<2) fprintf(stdout,"%10.3lf %20.10e %.10e %.10e\n",0,0.,populationvariance,populationsize);
  for(i=1;i<=maxSteps;i++) {
    reproduce(i);
    lambda = populationconstraint(i);
    
    if(i%outputstep == 0) {
      if (quiet<2)fprintf(stdout,"%10.3lf %20.10e %.10e %.10e %.10e\n",i*epsilon,allshifts*dx + popdens_1stmom/popdens_0thmom,populationvariance,populationsize,lambda);
      if (quiet==0) print_populationdensity(i);
      if (averagepopdens) update_averagepopdens();
    }
  }
  
  if(write_to_file)write_popdens();
  if(averagepopdens)write_averagepopdens();
  cleanup();
  return 0;
}


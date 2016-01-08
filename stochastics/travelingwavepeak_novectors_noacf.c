#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>


int space = 100;
int space0 = 50;
double dx = 1e-2;
int maxSteps = 1000;

int outputstep = 100;
int quiet = 0;

int noise = 0;
double populationsize = 1.;

int correctformeanfitness =0;

int allshifts = 0;

double current_mean_fitness = 0.;

int read_from_file = 0;
int write_to_file = 0;
char c_infile[128],c_outfile[128];


double epsilon = 1e-2;
double mutationrate = 1e-5;
double twoepssqrt;
double *nn;
double *tmp;

const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
unsigned long int randseed = 0;


int acf = 0;
double *acf_sum, *acf_count;
double *acf_fv_pastvalues;
double acf_maxtime;
double acf_fv_sum = 0., acf_fv_count = 0.;
int acf_maxindex;
char acf_outfile[128];


int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


void parsecomamndline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn,argv,"s:z:d:S:e:D:O:qQi:o:R:N:Ca:A:")) != -1) {
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
		noise = 1;
		break;
      case 'C':	correctformeanfitness = 1;
		break;
      case 'a':	acf = 1;
		strcpy(acf_outfile,optarg);
		break;
      case 'A':	if( acf == 0 ){
		  print_error("option -a OUTFILE necessary before -A");
		}else{
		  acf_maxtime = atof(optarg);
		}
		break;
    }
  }
  if(randseed==0)randseed=time(NULL);
}



void read_popdens() {
  int icount,dcount;
  int *ival;
  double *dval;
  FILE *fp;
  
  fp = fopen(c_infile,"rb");
  fclose(fp);
  
}



void write_popdens() {
}









int initialize() {
  int i;
  
  nn = (double*)calloc(space,sizeof(double));
  nn[space0] = populationsize;
  tmp = (double*)malloc(space*sizeof(double));

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);

  twoepssqrt = 2.*sqrt(epsilon);
  
  
  if(acf) {
    acf_maxindex      = (int)(acf_maxtime/(outputstep*epsilon));
    acf_count         = (double*)calloc(acf_maxindex,sizeof(double));
    acf_sum           = (double*)calloc(acf_maxindex,sizeof(double));
    acf_fv_pastvalues = (double*)malloc(acf_maxindex*sizeof(double));
  }
  
}


void update_acf(int timestep,double fv) {
  int acf_curindex = (timestep/outputstep)%acf_maxindex;
  int i;
  acf_fv_pastvalues[acf_curindex] = fv;
  
//   printf("update acf: step %d\n",timestep);
  
//   printf("maxindex = %d, curindex = %d\n",acf_maxindex,acf_curindex);
  
  acf_fv_sum += fv;
  acf_fv_count += 1.;

  if(timestep/outputstep > acf_maxindex) {
    
    
    for(i=acf_curindex;i>=0;i--) {
      acf_sum[acf_curindex-i] += fv*acf_fv_pastvalues[i];
      acf_count[acf_curindex-i] += 1.;
    }
    for(i=acf_maxindex;i>acf_curindex;i--) {
      acf_sum[acf_curindex-i+acf_maxindex] += fv*acf_fv_pastvalues[i];
      acf_count[acf_curindex-i+acf_maxindex] += 1.;
    }
  }
}


void write_acf() {
  int i;
  FILE *fp_acf;
  double fv2 = acf_fv_sum*acf_fv_sum/(acf_fv_count*acf_fv_count);
  double var0 = acf_sum[0]/acf_count[0] - fv2;
  
  fp_acf = fopen(acf_outfile,"w");
  if((fp_acf != NULL)&&(acf_count[0]>0)) {
    for(i=0;i<acf_maxindex;i++) {
      fprintf(fp_acf,"%lf %20.10e %d\n",i*outputstep*epsilon,(acf_sum[i]/acf_count[i] - fv2)/var0,(int)acf_count[i]);
    }
  }
  fclose(fp_acf);
}

  
void update_mean_fit() {
  int i;
  current_mean_fitness = 0.;
  for(i=0;i<space;i++)current_mean_fitness += (i-space0)*nn[i];
  current_mean_fitness *= dx/populationsize;
}
  
  

 
void print_populationdensity(int timestep) {
  int i;
  double corr = current_mean_fitness;
  if(correctformeanfitness)corr += allshifts*dx;
  for(i=0;i<space;i++) {
    fprintf(stderr,"%lf %14.10lf %20.10e\n",timestep*epsilon,(i-space0)*dx-corr,nn[i]);
  }
  fprintf(stderr,"\n");
}


double get_population_variance() {
  int i;
  double xn=0., xxn=0.;
  
  for(i=0;i<space;i++) {
    xn += i*dx*nn[i];
    xxn += i*i*dx*dx*nn[i];
  }
  
  return xxn/populationsize - xn*xn/(populationsize*populationsize);
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
  int shift = (int)floor(current_mean_fitness/dx);
  if(shift > 0) {shift_population_backward(shift);}
  if(shift < 0) {shift_population_forward(shift);}
  allshifts += shift;
  current_mean_fitness -= shift*dx;
}

void reproduce(int timestep) {
  int i;
  
  update_mean_fit();
  
  if((current_mean_fitness > dx) || (current_mean_fitness < -dx))shift_population(timestep);
    
  tmp[0]  = nn[0]*(1.+epsilon*(-space0*dx-current_mean_fitness-mutationrate));
  if(noise)tmp[0] += twoepssqrt*(gsl_ran_poisson(rg,tmp[0])-tmp[0]);
  
  for(i=1;i<space;i++) {
    tmp[i]  = nn[i]*(1.+epsilon*((i-space0)*dx-current_mean_fitness-mutationrate)) + nn[i-1]*epsilon*mutationrate;
    if(noise)tmp[i] += twoepssqrt*(gsl_ran_poisson(rg,tmp[i])-tmp[i]);
  }
  
  memcpy(&nn[0],&tmp[0],space*sizeof(double));
}


void populationconstraint() {
  double currentpopsize = 0., inv_ps;
  int i;
  for(i=0;i<space;i++) currentpopsize += nn[i];
  inv_ps = populationsize/currentpopsize;
  for(i=0;i<space;i++) nn[i] *= inv_ps;
}






  
void cleanup() {
  free(nn);
  free(tmp);
}


int main(int argn, char *argv[]) {
  int i;
  double v;
  int last_allshifts;
  
  parsecomamndline(argn,argv);
  initialize();

  for(i=0;i<maxSteps;i++) {
    
    reproduce(i);
    populationconstraint();
    
    if(i%outputstep == 0) {
//       if ((acf) || (quiet<2))v = get_population_variance();
//       if (acf) update_acf(i,v);
      if (quiet<2) {
	v = get_population_variance();
	fprintf(stdout,"%lf %14.10lf %14.10lf\n",i*epsilon,current_mean_fitness+allshifts*dx,v);
      }
      if (quiet==0) print_populationdensity(i);
    }
    last_allshifts = allshifts;
  }
  
//   if(acf)write_acf();
  
  cleanup();
  return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_erf.h>

int space,space0;
double dx;

double *c,*u;

double *uoriginal;
double *uexp;

int haveufile = 0;
char ufilename[128];
int havecfile = 0;
char cfilename[128];

int adjust_u_profile = 0;

double *fixedClones;
double *fixedMutations;
double *Background;
double *winningTicket;
double *haldanePi;

double fitnessvariance = 1.e-5;	// note: std-parameters do not match!!
double adaptationspeed = 1.e-5;
double mutationrate    = 1.e-5;
double sigma           = 1.e-2;
double populationsize  = 1.e7;


int normalizedistributions = 0;
double max_u = 0.,max_c = 0.;
double max_fixedClones = 0.;
double max_fixedMutations = 0.;
double max_Background = 0.;
double max_winningTicket = 0.;
double max_haldanePi = 0.;

int meanvalues = 0;
double mean_fixedClones;
double mean_fixedMutations;


double adjust_threshold = 1e-3;
int adjust_step = 2;

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"u:c:AV:v:N:D:s:MmT:t:")) != -1){
    switch(c) {
      case 'u':	strcpy(ufilename,optarg);
		haveufile = 1;
		break;
      case 'c':	strcpy(cfilename,optarg);
		havecfile = 1;
		break;
      case 'A':	adjust_u_profile = 1;
		break;
      case 'V':	fitnessvariance = atof(optarg);
		break;
      case 'v':	adaptationspeed = atof(optarg);
		break;
      case 'N':	populationsize = atof(optarg);
		break;
      case 'D':	mutationrate = atof(optarg);
		break;
      case 's':	sigma = atof(optarg);
		break;
      case 'M':	normalizedistributions = 1;
		break;
      case 'm':	meanvalues = 1;
		break;
      case 'T':	adjust_threshold = atof(optarg);
		break;
      case 't':	adjust_step = atoi(optarg);
		break;
    }
  }
  if(haveufile == 0)print_error("need at least u*-profile, option -u FILENAME");
}


int read_u() {
  int icount,dcount;
  int *ival;
  double *dval;
  FILE *fp;
  
  fp = fopen(ufilename,"rb");
  if(fp != NULL) {
    fread(&icount,sizeof(int),1,fp);
    fread(&dcount,sizeof(int),1,fp);
    fread(&dx,sizeof(double),1,fp);
    fread(&space,sizeof(int),1,fp);
    fread(&space0,sizeof(int),1,fp);
    
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fp);
      free(ival);
    }
    
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(double),dcount,fp);
      free(dval);
    }
    
    u = (double*)malloc(space*sizeof(double));
    fread(u,sizeof(double),space,fp);
    
    fclose(fp);
  }else{
    print_error("could not open file with u*-profile");
  }
}

  
int read_c() {
  // dummy
  print_error("not yet implemented");
}


int approximate_c() {
  int i;
  double prefactor = populationsize*dx/sqrt(2.*3.141*fitnessvariance);
  double x;
  
  c = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++) {
    x = (i-space0)*dx;
    c[i] = prefactor*exp(-0.5*x*x/fitnessvariance);
  }
}

int adjust() {
  int i;
  int step = adjust_step;
  int startindex_approximation = -1;
  double threshold2 = adjust_threshold*adjust_threshold;
  double diff = 1.;
  double decay;
  double *logu = (double*)malloc(space*sizeof(double));
  
  uoriginal = (double*)malloc(space*sizeof(double));
  memcpy(uoriginal,u,space*sizeof(double));
  
  for(i=space-1;i>=0;i--) {
    if(u[i]>0) {
      logu[i] = log(u[i]);
    }else{
      if(startindex_approximation == -1)startindex_approximation = i;
    }
  }

  if(startindex_approximation == -1)return -1;
  
  while(diff*diff > threshold2) {
    startindex_approximation += step;
    diff = logu[startindex_approximation] - 2.*logu[startindex_approximation+step] + logu[startindex_approximation+2*step];
//     fprintf(stderr,"%lf %lf\n",(startindex_approximation-space0)*dx,diff);
  }
  startindex_approximation += step;
  decay = (logu[startindex_approximation+step] - logu[startindex_approximation+2*step])/(step*dx);
  for(i=startindex_approximation-1;i>=0;i--) {
    u[i] = u[startindex_approximation]*exp(decay*(startindex_approximation-i)*dx);
  }
  
  free(logu);
  
  return startindex_approximation;
}


double get_g() {
  int i;
  double norm = 0.;
  fixedClones = (double*)malloc(space*sizeof(double));
  for(i=0;i<space;i++) {
    fixedClones[i] = u[i]*c[i];
    norm += fixedClones[i];
  }
  for(i=0;i<space;i++)fixedClones[i] /= norm;
  return norm;
}



double get_fixedMutations(int startindex_approximation) {
  int i,j;
  double xms;
  double c_prefactor = populationsize/sqrt(2.*3.1315*fitnessvariance);
  
  fixedMutations = (double*)malloc(space*sizeof(double));
  
  // for x << 0, u*(x) ~ A * Exp(B*x)
  // A ... u_expprefactor
  // B ... u_expdecay
  
  double u_expdecay = (log(u[startindex_approximation]) - log(u[startindex_approximation-adjust_step]))/(adjust_step*dx);
  double u_expprefactor = u[startindex_approximation] * exp(-u_expdecay*(startindex_approximation-space0)*dx);
  
  double lowerbound = -space0*dx;
  double upperbound = (space-space0)*dx;
  
  double s;
  
  double correction_lowerbound;
  double correction_upperbound;

  uexp = (double*)malloc(space*sizeof(double));
  haldanePi = (double*)malloc(space*sizeof(double));

  
  for(i=0;i<space;i++) {
    
    uexp[i] = u_expprefactor*exp(s*u_expdecay);
    s = (i-space0)*dx;
    correction_lowerbound = 0.5*populationsize*u_expprefactor*exp(u_expdecay*s)*exp(u_expdecay*u_expdecay*fitnessvariance*0.5)*gsl_sf_erfc( (s+u_expdecay*fitnessvariance - lowerbound)/sqrt(2.*fitnessvariance));
//     fixedMutations[i] = 0.5*populationsize*u_expprefactor*exp(u_expdecay*s)*exp(u_expdecay*u_expdecay*fitnessvariance*0.5)*gsl_sf_erfc( (s+u_expdecay*fitnessvariance - lowerbound)/sqrt(2.*fitnessvariance));
    fixedMutations[i] = 0.;
    for(j=0;j<space;j++) {
      xms = (j-space0)*dx - s;
      fixedMutations[i] += u[j]* (c_prefactor*exp(-.5*xms*xms/fitnessvariance)) * dx;
    }
    correction_upperbound = 0.25*populationsize*(sqrt(2.*fitnessvariance/3.1415)*exp(- (s-upperbound)*(s-upperbound)/(2.*fitnessvariance)) + s*(1.+gsl_sf_erf((s-upperbound)/sqrt(2.*fitnessvariance))));
    
    
    fixedMutations[i] += correction_lowerbound + correction_upperbound;
    
    haldanePi[i] = fixedMutations[i];
    
    if(i>=space0) {
      fixedMutations[i] *= 1/sigma * exp(-(i-space0)*dx/sigma);
    }else{
      fixedMutations[i] *= 0.;
    }
    
//     fprintf(stderr,"%d %e %e\n",i,correction_lowerbound,correction_upperbound);
  }
}


void cleanup() {
  free(u);
  free(c);
  free(fixedClones);
  free(fixedMutations);
  free(uexp);
  free(uoriginal);
}


int main(int argn, char *argv[]) {
  int i;
  
  double sum_xg = 0.,sum_g = 0.;
  double sum_xf = 0.,sum_f = 0.;
  double x;
  
  int tmp;
  
  parsecommandline(argn,argv);
  
  read_u();
  
  if(havecfile == 1) {
    read_c();
  }else{
    approximate_c();
  }
  
  if(adjust_u_profile)tmp = adjust();
  
  double a = get_g();
  
  get_fixedMutations(tmp);
  
    
  if(normalizedistributions) {
    for(i=0;i<space;i++) {
      if(c[i] > max_c)max_c = c[i];
      if(u[i] > max_u)max_u = u[i];
      if(fixedClones[i] > max_fixedClones) max_fixedClones = fixedClones[i];
      if(fixedMutations[i] > max_fixedMutations) max_fixedMutations = fixedMutations[i];
      if(haldanePi[i] > max_haldanePi) max_haldanePi = haldanePi[i];
    }
    for(i=0;i<space;i++) {
      c[i] /= max_c;
      u[i] /= max_u;
      fixedClones[i] /= max_fixedClones;
      fixedMutations[i] /= max_fixedMutations;
      haldanePi[i] /= max_haldanePi;
    }
  }
  
  
  if(meanvalues) {
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      sum_f  +=     fixedMutations[i];
      sum_xf += x * fixedMutations[i];
      sum_g  +=     fixedClones[i];
      sum_xg += x * fixedClones[i];
    }
    fprintf(stderr,"%e %e\n",sum_xg/sum_g,sum_xf/sum_f);
  }
  
  
  for(i=0;i<space;i++) {
    printf("%lf\t%e\t%e\t%e\t%e\t%e",(i-space0)*dx,u[i],c[i],fixedClones[i],fixedMutations[i],haldanePi[i]);
    if(adjust_u_profile)printf("\t%e\t%e",uoriginal[i],uexp[i]);
    printf("\n");
  }
  
  
  
  return 0;
}

  
  


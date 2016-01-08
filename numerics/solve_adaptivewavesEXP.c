/* **********************************************************************
 Solve for numerical profiles of fixation probabilities u and
 mean stationary population density c in tuned models
 
 Mutations are drawn from exponential kernel with only
 beneficial mutations

 Algorithm is based on transformed equation:
  (Differential operator (+- 1 + \sigma\partial_X) is applied to
  governing equation for c and u)

 

 compiling:
    gcc -o solve_adaptivewavesEXP solve_adaptivewavesEXP.c -lm -lgsl
    
 usage:
    ./solve_adaptivewavesEXP -U -s 2000 -z 1000 -d 1e-5 -M 1e-3 -D 1e-9 -V 2e-6 -o u.conf -S 1000000
    ./solve_adaptivewavesEXP -C -u u.conf -o c.conf -S 1000000
  
  
 2014, Lukas Geyrhofer, lukas.geyrhofer@ds.mpg.de
  
********************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_exp.h>

// generic algorithm parameters
int maxsteps = 10000;
double alpha = 1.;
int usefivepointstencil = 0;
int adjust_step = 50;
int derivativestencil = 3;

// quiet = 0:	print everything
// quiet = 1:	print only options        (cmdline flag -Q)
// quiet = 2:	print final configuration (cmdline flag -q)
int quiet = 0;

// switchc == 1 -> calculate density c(x)
// switchu == 1 -> calculate u*(x)
int switchc = 0;
int switchu = 0;

// names and flags for binary input- and output-file
char c_infile[128], c_outfile[128];
char u_infile[128], u_outfile[128];
int have_c_infile = 0, have_c_outfile = 0;
int have_u_infile = 0, have_u_outfile = 0;

// lattice
double dx = 1e-4;
int space = 2000,space0 = 1000;

// model parameters
double variance = 1e-6;

// mutations
double diffusionconstant = 1e-10;
double sigma = 1e-3;
double mutationalbias = 0.;

// eigenvalue equation,  Lc = Gc (where L ... Louivillean, G ... eigenvalue Gamma)
// default: evgamma = 0 ... stationary solution
double evgamma = 0.;
int usegamma = 0;

// initial conditions
int initialcond = 1; // 1 - known approximations, 2 - random, 3 - flat, 4 - gaussian (only c, unavailable for u)
double maxval;

double *x,*xg;

double *u;
double *fu,*f;

double *c;
double *fc;

// ===========================================================================================
// print_error
// ===========================================================================================

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}

int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"UCi:o:u:I:S:a:s:d:z:M:V:D:QqA:fT:G:")) != -1){
    switch(c) {
      case 'C':	if(switchu >= 1) {
		  print_error("can only use one of the options -U -C");
		}else{
		  switchc = 1;
		}
		break;
      case 'U':	if(switchc >= 1) {
		  print_error("can only use one of the options -U -C");
		}else{
		  switchu = 1;
		}
		break;
      case 'i':	if(switchc == 1) {
		  strcpy(c_infile,optarg);
		  have_c_infile = 1;
		}else if(switchu == 1) {
		  strcpy(u_infile,optarg);
		  have_u_infile = 1;
		}else{
		  print_error("option -C or -U has to be specified first");
		}
		break;
      case 'o':	if(switchc == 1) {
		  strcpy(c_outfile,optarg);
		  have_c_outfile = 1;
		}else if(switchu == 1) {
		  strcpy(u_outfile,optarg);
		  have_u_outfile = 1;
		}else{
		  print_error("option -C or -U has to be specified first");
		}
		break;
      case 'u':	if(switchc == 1) {
		  strcpy(u_infile,optarg);
		  have_u_infile = 1;
		}else{
		  print_error("option -u requires option -C");
		}
		break;
		
      case 'M':	sigma = atof(optarg);
		break;
      case 'I':	if(strcmp(optarg,"approx") == 0) {
		  initialcond = 1;
		}else if(strcmp(optarg,"random") == 0) {
		  initialcond = 2;
		}else if(strcmp(optarg,"flat") == 0) {
		  initialcond = 3;
		}else if(strcmp(optarg,"zero") == 0) {
		  if(switchu==1) {
		    initialcond = 4;
		  }else{
		    print_error("'zero' initial conditions only available for u(x), use option -U");
		  }
		}else{
		  print_error("argument provided to option -I is not valid: choose from 'approx', 'random' or 'flat'");
		}
		break;
		
      case 'S':	maxsteps = atoi(optarg);
		break;
      case 'a':	alpha = atof(optarg);
		break;

      case 'd':	dx = atof(optarg);
		break;
      case 's':	space = atoi(optarg);
		break;
      case 'z':	space0 = atoi(optarg);
		break;
		
      case 'V':	variance = atof(optarg);
		break;
      case 'D':	diffusionconstant = atof(optarg);
		break;
      case 'Q':	quiet = 1;
		break;
      case 'q':	quiet = 2;
		break;
      case 'A':	if(switchc == 1) {
		  adjust_step = atoi(optarg);
		}
		break;
      case 'T':	derivativestencil = atoi(optarg);
		if( (derivativestencil != 3) && (derivativestencil != 5)){
		  print_error("only stencils of order 3 and 5 implemented yet");
		}
		break;
      case 'G': if(switchc == 1) {
		  usegamma = 1;
		  evgamma = atof(optarg);
		}else{
		  print_error("Eigenvalue Gamma can only used with population density (mode -C)");
		}
		break;
    }
  }
  if(switchc + switchu == 0) {
    print_error("have to use at least one of the options, -C -U");
  }
  if(switchc== 1) {
    if(have_u_infile == 0)print_error("need u* file to calculate density profile c(x), option -u FILE");
  }
  return 0;
}


void initial_c_approx() {
  int i;
  for(i=0;i<space;i++) {
    if(usegamma == 1) {
      c[i] = exp(-xg[i]*xg[i]/(2.*variance));
    }else{
      c[i] = exp(-x[i]*x[i]/(2.*variance));
    }
  }
}

void initial_c_flat() {
  int i;
  for(i=0;i<space;i++)c[i] = 1.;
}

void initial_c_random() {
  int i;
  for(i=0;i<space;i++)c[i] = rand();
}

void initial_u_approx() {
  int i;
  for(i=0;i<=space0;i++)u[i] = 0.5*dx*exp(x[i]/sigma);
  for(i=space0+1;i<space;i++)u[i] = 0.5*x[i];
}

void initial_u_flat() {
  int i;
  for(i=0;i<space;i++)u[i] = dx;
}

void initial_u_random() {
  int i;
  for(i=0;i<space;i++)u[i] = (dx*rand())/(1.*RAND_MAX);
}


void initial_u_zero() {
  int i;
  for(i=0;i<space0;i++)u[i] = 0.;
  for(i=space0;i<space;i++)u[i] = 0.5*x[i];
}


// ===========================================================================================
// read external files
// ===========================================================================================

void read_u_file(int extractparameters) {
  FILE *fpu;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(u_infile,"rb");
  if(fpu != NULL) {
    fread(&icount,1,sizeof(int),fpu);
    fread(&dcount,1,sizeof(int),fpu);
    fread(&u_dx,1,sizeof(double),fpu);
    fread(&u_space,1,sizeof(int),fpu);
    fread(&u_space0,1,sizeof(int),fpu);
    if(icount > 0) {
      tmpi = (int*)malloc(icount*sizeof(int));
      fread(&tmpi[0],icount,sizeof(int),fpu);
      free(tmpi);
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpu);
      
      if((dcount >= 3)&&(extractparameters==1)) {
	diffusionconstant = tmpd[0]*tmpd[2]*tmpd[2];
	variance = tmpd[1] - tmpd[0]*tmpd[2];
	sigma = tmpd[2];
	
	// for reference: other programs use the following paramters, not the ones used above
	// 	mutation_rate     = tmpd[0];
	// 	wavespeed         = tmpd[1];
	// 	mutation_expdecay = tmpd[2];
      }
      free(tmpd);
    }
    if(extractparameters == 1) {
      dx = u_dx;
      space = u_space;
      space0 = u_space0;
    }else if (( space != u_space) || (space0 != u_space0)) {
      print_error("lattice does not match");
    }
    u=(double*)malloc(space*sizeof(double));;
    fread(u,sizeof(double),space,fpu);
    fclose(fpu);
  }else{
    print_error("could not open u-file");
  }
}
  
void read_c_file(int extractparameters) {
  FILE *fpc;
  int icount,dcount;
  int *ival;
  double *dval;
  int c_space,c_space0;
  double c_dx;
  
  fpc = fopen(c_infile,"rb");
  if(fpc != NULL) {
    fread(&icount,sizeof(int),1,fpc);
    fread(&dcount,sizeof(int),1,fpc);
    fread(&c_dx,sizeof(double),1,fpc);
    fread(&c_space,sizeof(int),1,fpc);
    fread(&c_space0,sizeof(int),1,fpc);
    
    if(icount > 0) {
      ival = (int*)malloc(icount*sizeof(int));
      fread(ival,sizeof(int),icount,fpc);
      free(ival);
    }
    
    if(dcount > 0) {
      dval = (double*)malloc(dcount*sizeof(double));
      fread(dval,sizeof(double),dcount,fpc);
      free(dval);
    }
    
    if(extractparameters) {
      space  = c_space;
      space0 = c_space0;
      dx     = c_dx;
    }
    
    c = (double*)malloc(space*sizeof(double));
    fread(c,sizeof(double),space,fpc);
    
    fclose(fpc);
  }else{
    print_error("could not open c-file");
  }
}

// ===========================================================================================
// write output files
// ===========================================================================================

void write_u_file() {
  FILE *fpu;
  int icount = 0,dcount = 3;
  double *tmpd;
  
  tmpd = (double*)malloc(3*sizeof(double));
  tmpd[0] = diffusionconstant/(sigma*sigma);
  tmpd[1] = variance + diffusionconstant/sigma;
  tmpd[2] = sigma;
  
  fpu = fopen(u_outfile,"wb");
  if(fpu != NULL) {
    fwrite(&icount,1,sizeof(int),fpu);
    fwrite(&dcount,1,sizeof(int),fpu);
    fwrite(&dx,1,sizeof(double),fpu);
    fwrite(&space,1,sizeof(int),fpu);
    fwrite(&space0,1,sizeof(int),fpu);
    fwrite(&tmpd[0],sizeof(double),3,fpu);
    fwrite(u,sizeof(double),space,fpu);
    fclose(fpu);
  }
  free(tmpd);
}

void write_c_file() {
  FILE *fpc;
  int icount = 0,dcount = 3 + usegamma;
  double *tmpd;
  
  
  
  tmpd = (double*)malloc((3+usegamma)*sizeof(double));
  tmpd[0] = diffusionconstant/(sigma*sigma);
  tmpd[1] = variance + diffusionconstant/sigma;
  tmpd[2] = sigma;
  if(usegamma) {
    tmpd[3] = evgamma;
  }
    

  fpc = fopen(c_outfile,"wb");
  if(fpc != NULL) {
    fwrite(&icount,1,sizeof(int),fpc);
    fwrite(&dcount,1,sizeof(int),fpc);
    fwrite(&dx,1,sizeof(double),fpc);
    fwrite(&space,1,sizeof(int),fpc);
    fwrite(&space0,1,sizeof(int),fpc);
    fwrite(&tmpd[0],sizeof(double),3+usegamma,fpc);
    fwrite(c,sizeof(double),space,fpc);
    fclose(fpc);
  }
  free(tmpd);
}



// ===========================================================================================
// derivatives
// ===========================================================================================

double derivativeu(int i,int order,int stencil) {
  switch(order){
    case 1:	switch(stencil) {
		  case 3:	return (u[i+1] - u[i-1])/(dx); break;
		  case 5:	return (-u[i+2] + 8.*u[i+1] -8.*u[i-1] + u[i-2])/(12.*dx); break;
		}
		break;
    case 2:	switch(stencil) {
		  case 3:	return (u[i+1] - 2.*u[i] + u[i-1])/(dx*dx); break;
		  case 5:	return (-u[i+2] + 16.*u[i+1] -30.*u[i] + 16.*u[i-1] - u[i-2])/(12.*dx*dx); break;
		}
		break;
  }
}


// ===========================================================================================
// iteration procedures
// ===========================================================================================

int iterate_u (int timestep) {
				/* nonlinear red-black gauss-seidel iteration */
				/* taken from oskar's code ... */
  int i;
  double du1,du2,duu1;
  
//   if(u[3]>0) {
//     u[1]=u[2] * u[2]/u[3];
//   }else{
//     u[1] = 0;
//   }
  if(u[2]>0) {
    u[0]=u[1] * u[1]/u[2];
  }else{
    u[0] = 0;
  }
//   u[0] = 0;
  u[space-1] = 2*u[space-2]-u[space-3];

  
  // use known asymptotic solution
//   u[space-2]=0.5*(x[space-2]-variance/x[space-2]);
//   u[space-1]=0.5*(x[space-1]-variance/x[space-1]);
  
//   u[space-2] = 2*u[space-3]-u[space-2

  for (i=1; i<space-1; i+=2) {	/* then odd-numbered items */

    du2   = derivativeu(i,2,derivativestencil);
    du1   = derivativeu(i,1,derivativestencil);
    duu1  = 2.*u[i]*derivativeu(i,1,derivativestencil);

    f[i]  = (sigma*variance + diffusionconstant)*du2 - (sigma*x[i]+variance)*du1 + (x[i]-sigma)*u[i] - 2.*u[i]*u[i] + 2.*sigma*duu1;
    fu[i] = -2.*(sigma*variance + diffusionconstant)/(dx*dx) + (x[i] - sigma) - 4.*u[i] + 4*sigma*du1;

    u[i]=u[i]-alpha*f[i]/fu[i];
    if (u[i]>maxval)u[i]=maxval;
    if (u[i]<1e-300)u[i]=0;

  }

  
  for (i=2; i<space-1; i+=2) {	/* first even-numbered items */

    du2   = derivativeu(i,2,derivativestencil);
    du1   = derivativeu(i,1,derivativestencil);
    duu1  = 2.*u[i]*derivativeu(i,1,derivativestencil);

    f[i]  = (sigma*variance + diffusionconstant)*du2 - (sigma*x[i]+variance)*du1 + (x[i]-sigma)*u[i] - 2.*u[i]*u[i] + 2.*sigma*duu1;
    fu[i] = -2.*(sigma*variance + diffusionconstant)/(dx*dx) + (x[i] - sigma) - 4.*u[i] + 4*sigma*du1;

    u[i]=u[i]-alpha*f[i]/fu[i];
    if (u[i]>maxval)u[i]=maxval;
    if (u[i]<1e-300)u[i]=0;

  }
  return 0;
}


int iterate_c (int timestep) {
				/* nonlinear red-black gauss-seidel iteration */
  int i;
  double dc1,dc2,duc1;

  c[0]       = 0;
  c[1]       = 0;
  c[space-2] = 0;
  c[space-1] = 0;
  
  for(i=2;i<space-2;i+=2) {
    dc1  = 0.5*(c[i+1] - c[i-1])/dx;
    dc2  = (c[i+1] - 2.*c[i] + c[i-1])/(dx*dx);
    duc1 = 0.5*(u[i+1]*c[i+1] - u[i-1]*c[i-1])/dx;
    
    f[i]  = (sigma*variance + diffusionconstant)*dc2 + (variance + sigma*xg[i])*dc1 + (xg[i] + sigma)*c[i] - 2.*u[i]*c[i] - 2.*sigma*duc1;
    fc[i] = -2.*(sigma*variance + diffusionconstant)/(dx*dx) + (xg[i] + sigma - 2.*u[i]);
    
    c[i] = c[i] - alpha*f[i]/fc[i];
    if((c[i]<0)&&(usegamma==0))c[i]=0.; // restriction to positive values only for stationary solution
  }

  for(i=3;i<space-2;i+=2) {
    dc1  = 0.5*(c[i+1] - c[i-1])/dx;
    dc2  = (c[i+1] - 2.*c[i] + c[i-1])/(dx*dx);
    duc1 = 0.5*(u[i+1]*c[i+1] - u[i-1]*c[i-1])/dx;
    
    f[i]  = (sigma*variance + diffusionconstant)*dc2 + (variance + sigma*xg[i])*dc1 + (xg[i] + sigma)*c[i] - 2.*u[i]*c[i] - 2.*sigma*duc1;
    fc[i] = -2.*(sigma*variance + diffusionconstant)/(dx*dx) + (xg[i] + sigma - 2.*u[i]);
    
    c[i] = c[i] - alpha*f[i]/fc[i];
    if((c[i]<0)&&(usegamma==0))c[i]=0.; // restriction to positive values only for stationary solution
  }
  return 0;
}




void adjust_constraint(int timestep) {
  int i;
  double sum = 0.,invsum;
  for(i=0;i<space;i++)sum += u[i]*c[i];
  invsum = 1./(sum*dx);
  for(i=0;i<space;i++)c[i] *= invsum;
}


// ===========================================================================================
// initialize all variables
// ===========================================================================================

void initialize() {
  int i,j;
  double y,tmp;

  
  // ==================================================
  // initilization for U mode
  // ==================================================
  
  if (switchu == 1) {
    if (have_u_infile == 1) {
      read_u_file(1);
      x = (double*)malloc(space*sizeof(double));
      for(i=0;i<space;i++) x[i] = (i-space0)*dx;
    }else{
      u = (double*)malloc(space*sizeof(double));
      x = (double*)malloc(space*sizeof(double));
      for(i=0;i<space;i++) x[i] = (i-space0)*dx;
      switch(initialcond) {
	case 1:	initial_u_approx();
		break;
	case 2:	initial_u_random();
		break;
	case 3:	initial_u_flat();
		break;
	case 4:	initial_u_zero();
		break;
	default:print_error("something went wrong, option for initial condition not available");
      }
    }
    
    f = (double*)malloc(space*sizeof(double));
    fu = (double*)malloc(space*sizeof(double));
  }

  
  
  // ==================================================
  // initialization for C mode
  // ==================================================
  
  if(switchc == 1) {
    read_u_file(1);
    x=(double*)malloc(space*sizeof(double));
    xg = (double*)malloc(space*sizeof(double));
    for(i=0;i<space;i++) x[i] = (i-space0)*dx;
    memcpy(xg,x,space*sizeof(double));
    if(usegamma == 1) {
      for(i=0;i<space;i++) xg[i] -= evgamma;
    }
      
    if( have_c_infile == 1) {
      read_c_file(0);
    }else{
      c  = (double*)malloc(space*sizeof(double));
      
      switch(initialcond) {
	case 1:	initial_c_approx();
		break;
	case 2:	initial_c_random();
		break;
	case 3:	initial_c_flat();
		break;
	default:print_error("something went wrong, option for initial condition not available");
      }
    }
    
    adjust_constraint(0);
    
    f  = (double*)malloc(space*sizeof(double));
    fc = (double*)malloc(space*sizeof(double));
  }
  
  // set maximal value the algorithm admits for a value of u[i]
  // 10x larger values than expected allowed, then culled
  maxval = 5.*(space-space0)*dx;
  
}



void print_options() {
  printf("###################################################\n");
  printf("# solve profiles and weight-functions numerically #\n");
  printf("###################################################\n");  
  printf("#\n");
  printf("# LATTICE\n");
  printf("#   space           = %d\n",space);
  printf("#   space0          = %d\n",space0);
  printf("#   dx              = %lf\n",dx);
  printf("# mutationmodel     = exp, transformed equation\n");
  printf("# variance          = %e\n",variance);
  printf("# diffusionconst    = %e\n",diffusionconstant);
  printf("# sigma             = %e\n",sigma);
  if(usegamma == 1) {
    printf("# Gamma             = %e\n",evgamma);
  }
    
}

void print_u() {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %e\n",x[i],u[i]);
  }
}

void print_c() {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %e\n",x[i],c[i]);
  }
}



// ===========================================================================================
// cleanup
// ===========================================================================================

void cleanup() {
  int i;
  if(switchu == 1) {
    free(u);
    free(fu);
    free(f);
    free(x);
  }else if(switchc == 1){
    free(c);
    free(fc);
    free(f);
    free(x);
    free(xg);
  }
}



// ===========================================================================================
// main
// ===========================================================================================

int main(int argn, char *argv[]) {
  int i;
  
  
  parsecommandline(argn,argv);
  initialize();
  
  
  if(quiet<1)print_options();
  if(switchu == 1) {
    for(i=1;i<=maxsteps;i++) {
      iterate_u(i);
    }
    if(have_u_outfile == 1)write_u_file();
    if(quiet<2)print_u();
  }else if(switchc == 1) {
    for(i=1;i<=maxsteps;i++) {
      iterate_c(i);
      if(adjust_step>0)if(i%adjust_step==0)adjust_constraint(i);
    }
    adjust_constraint(i);
    if(have_c_outfile == 1)write_c_file(i);
    if(quiet<2)print_c();
  }
  cleanup();
  return 0;
}



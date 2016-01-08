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

int quiet = 0;


// switchc == 1 -> calculate density c(x)
// switchu == 1 -> calculate u*(x)
// switchw == 1 -> calculate w(x), u*(x) = w(x)^2, hence positive
// switchx == 1 -> calculate spectral properties of linear operator for c(x)
int switchc = 0;
int switchu = 0;
int switchx = 0;
int switchw = 0;

// names and flags for binary input- and output-file
char c_infile[128], c_outfile[128];
char u_infile[128], u_outfile[128];
int have_c_infile = 0, have_c_outfile = 0;
int have_u_infile = 0, have_u_outfile = 0;

// lattice
double dx = 1e-4;
int space = 10000,space0 = 5000;


// densities, etc...
// double *c,*cnew;
// double *u,*unew;

// double *c_coeff_p,*c_coeff_0,*c_coeff_m,*c_coeff_j;
// double *u_coeff_p,*u_coeff_0,*u_coeff_m,*u_coeff_j;


// model parameters
double wavespeed = 0.0001;


// mutations
int    mutation_type = 3; // 1 - diffusion, 2 - jumps, 3 - expdecay, 4 - external file
double mutation_rate = 0.00001;
int    mutation_jumpwidth = 20;
char   mutation_kernelfile[128];
double *mutation_kernel;
double mutation_expdecay = 0.01;
int    allowmutationoptions = 0;
double mutationalbias = 0.;
// initial conditions
int initialcond = 1; // 1 - known approximations, 2 - random, 3 - flat, 4 - gaussian (only c, unavailable for u)


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
  while((c = getopt(argn, argv,"UCi:o:u:I:S:a:s:d:z:m:M:v:D:QqA:f")) != -1){
    switch(c) {
      case 'C':	if(switchu+switchx+switchw >= 1) {
		  print_error("can only use one of the options -U -C -X -W");
		}else{
		  switchc = 1;
		}
		break;
      case 'U':	if(switchc+switchx+switchw >= 1) {
		  print_error("can only use one of the options -U -C -X -W");
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
		  print_error("option -C or -U or -W has to be specified first");
		}
		break;
      case 'o':	if(switchc == 1) {
		  strcpy(c_outfile,optarg);
		  have_c_outfile = 1;
		}else if(switchu + switchw == 1) {
		  strcpy(u_outfile,optarg);
		  have_u_outfile = 1;
		}else{
		  print_error("option -C or -U or -W has to be specified first");
		}
		break;
      case 'u':	if(switchc+switchx == 1) {
		  strcpy(u_infile,optarg);
		  have_u_infile = 1;
		}else{
		  print_error("option -u requires option -C/-X");
		}
		break;
		
      case 'm':	if(strcmp(optarg,"diffusion") == 0) {
		  mutation_type = 1;
		  allowmutationoptions = 1;
		}else if(strcmp(optarg,"jump") == 0) {
		  mutation_type = 2;
		  allowmutationoptions = 1;
		}else if(strcmp(optarg,"expdecay") == 0) {
		  mutation_type = 3;
		  allowmutationoptions = 1;
		}else if(strcmp(optarg,"file") == 0) {
		  mutation_type = 4;
		  allowmutationoptions = 1;
		}else{
		  print_error("argument provided to -m option is not valid: choose from 'diffusion', 'jump', 'expdecay' or 'file'");
		}
		break;
      case 'M':	if(allowmutationoptions == 1) {
		  switch(mutation_type) {
		    case 1:	print_error("'diffusion' mutation type does not have additional parameters");
				break;
		    case 2:	mutation_jumpwidth = atoi(optarg);
				break;
		    case 3:	mutation_expdecay = atof(optarg);
				break;
		    case 4:	strcpy(mutation_kernelfile,optarg);
				break;
		  }
		}else{
		  print_error("mutation type is not specified yet, use option -m");
		}
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
		
      case 'v':	wavespeed = atof(optarg);
		break;
      case 'D':	mutation_rate = atof(optarg);
		break;
      case 'Q':	quiet = 1;
		break;
      case 'q':	quiet = 2;
		break;
      case 'A':	if(switchc == 1) {
		  adjust_step = atoi(optarg);
		}
		break;
      case 'f':	usefivepointstencil = 1;
		break;
    }
  }
  if(switchc + switchu + switchx + switchw== 0) {
    print_error("have to use at least one of the options, -C -U -X -W");
  }
  if(switchc+switchx == 1) {
    if(have_u_infile == 0)print_error("need u* file to calculate density profile c(x)");
  }
  return 0;
}


void initial_c_approx() {
  int i;
  double x;
  for(i=0;i<space;i++) {
    x = (i-space0)*dx;
    c[i] = exp(-x*x/(2.*wavespeed));
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
  double x;
  if(mutation_type == 3) {
    for(i=0;i<=space0;i++)u[i] = 0.5*dx*gsl_sf_exp((i-space0)*dx/mutation_expdecay);
  }else{
    for(i=0;i<=space0;i++)u[i] = dx;
  }
  for(i=space0+1;i<space;i++)u[i] = 0.5*(i-space0)*dx;
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
  for(i=space0;i<space;i++)u[i] = 0.5*(i-space0)*dx;
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
      
      if((dcount > 4)&&(extractparameters==1)) {
	mutation_rate = tmpd[0];
	wavespeed = tmpd[1];
	mutation_expdecay = tmpd[2];
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
  fpu = fopen(u_outfile,"wb");
  fwrite(&icount,1,sizeof(int),fpu);
  fwrite(&dcount,1,sizeof(int),fpu);
  fwrite(&dx,1,sizeof(double),fpu);
  fwrite(&space,1,sizeof(int),fpu);
  fwrite(&space0,1,sizeof(int),fpu);
  fwrite(&mutation_rate,1,sizeof(double),fpu);
  fwrite(&wavespeed,1,sizeof(double),fpu);
  fwrite(&mutation_expdecay,1,sizeof(double),fpu);
  fwrite(u,sizeof(double),space,fpu);
  fclose(fpu);
}

void write_c_file() {
  FILE *fpc;
  int icount = 0,dcount = 0;
  fpc = fopen(c_outfile,"wb");
  if(fpc != NULL) {
    fwrite(&icount,1,sizeof(int),fpc);
    fwrite(&dcount,1,sizeof(int),fpc);
    fwrite(&dx,1,sizeof(double),fpc);
    fwrite(&space,1,sizeof(int),fpc);
    fwrite(&space0,1,sizeof(int),fpc);
    fwrite(c,sizeof(double),space,fpc);
    fclose(fpc);
  }
}

// ===========================================================================================
// iteration procedures
// ===========================================================================================

int iterate_u (int timestep) {
				/* nonlinear red-black gauss-seidel iteration */
				/* from oskar.... */
  int i;
  double p2,p1,p3,p4,p5;
  double x;


  /* u[1]=2*u[2]-u[3]; */
  /* u[0]=2*u[1]-u[2]; */

  /* for (i=1; i<k+1; i+=1) { */
  /*     u[N+1-i]=2*u[N-i]-u[N-i-1]; */
  /* } */

//   for (i=2; i>-1; i-=1) {
//     u[i]=0;
//   }
     u[1]=u[2] * u[2]/u[3];
     u[0]=u[1] * u[1]/u[2];

  u[space-2]=2*u[space-3]-u[space-4]; 
  u[space-1]=2*u[space-2]-u[space-3]; 
  				/* if f is to be zero, we iterate p <- p-f/f' */
				
				
//   x = (-space0)*dx;
//   p2 = (u[1] - 2.*u[0] + gsl_sf_exp(-dx/mutation_expdecay)*u[0])/(dx*dx);
//   p3 = (u[1] - gsl_sf_exp(-dx/mutation_expdecay)*u[0])/(2.*dx);
//   p4 = (u[1]*u[1] - gsl_sf_exp(-2.*dx/mutation_expdecay)*u[0]*u[0])/(2.*dx);
//   
//   f[0]=wavespeed*p2 - (x+wavespeed/mutation_expdecay-mutation_rate)*p3 + (x/mutation_expdecay-1)*u[0] + 2.*p4 - 2.*u[0]*u[0]/mutation_expdecay;
//   fu[0]=-wavespeed*2./(dx*dx)+(x/mutation_expdecay-1.)-4.*u[0]/mutation_expdecay;
// 
//   u[0]=u[0]-alpha*f[0]/fu[0];
//   if (u[0]>2)u[0]=2;
//   if (u[0]<0)u[0]=0;
// 

  for (i=2; i<space-2; i+=2) {	/* first even-numbered items */

    x=(i-space0)*dx;
    p2=(u[i+1]-2.*u[i]+u[i-1])/(dx*dx); /* 2nd derivative*/
    p3=(u[i+1]-u[i-1])/(2.*dx); /* 1st derivative */
    p4=(u[i+1]*u[i+1]-u[i-1]*u[i-1])/(2.*dx); /*  derivative of u^2*/

    f[i]=wavespeed*p2 - (x+wavespeed/mutation_expdecay-mutation_rate)*p3 + (x/mutation_expdecay-1)*u[i] + 2.*p4 - 2.*u[i]*u[i]/mutation_expdecay;
    fu[i]=-wavespeed*2./(dx*dx)+(x/mutation_expdecay-1.)-4.*u[i]/mutation_expdecay;

    u[i]=u[i]-alpha*f[i]/fu[i];
    if (u[i]>2)u[i]=2;
    if (u[i]<0)u[i]=0;

  }
  for (i=3; i<space-2; i+=2) {	/* then odd-numbered items */

    x=(i-space0)*dx;
    p2=(u[i+1]-2.*u[i]+u[i-1])/(dx*dx); /* 2nd derivative*/
    p3=(u[i+1]-u[i-1])/(2.*dx); /* 1st derivative */
    p4=(u[i+1]*u[i+1]-u[i-1]*u[i-1])/(2.*dx); /*  derivative of u^2*/

    f[i]=wavespeed*p2 - (x+wavespeed/mutation_expdecay-mutation_rate)*p3 + (x/mutation_expdecay-1)*u[i] + 2.*p4 - 2.*u[i]*u[i]/mutation_expdecay;
    fu[i]=-wavespeed*2./(dx*dx)+(x/mutation_expdecay-1.)-4.*u[i]/mutation_expdecay;

    u[i]=u[i]-alpha*f[i]/fu[i];
    if (u[i]>2)u[i]=2;
    if (u[i]<0)u[i]=0;

  }
  return 0;
}


int iterate_c (int timestep) {
				/* nonlinear red-black gauss-seidel iteration */
  int i;
  double p2,p1,p3,p4;
  double x;

  c[0]       = 0;
  c[1]       = 0;
  c[space-2] = 0;
  c[space-1] = 0;

  /* c[N-1]=2*c[N-2]-c[N-3];  */
  /* c[N]=2*c[N-1]-c[N-2];  */

  for (i=2; i<space-2; i+=2) {	/* first even-numbered items */
    x=(i-space0)*dx; 
    p2=(c[i+1]-2.*c[i]+c[i-1])/(dx*dx); /* 2nd derivative */
    p3=(c[i+1]-c[i-1])/(2*dx); /*1st drift derivative */
    p4=(c[i+1]*u[i+1]-c[i-1]*u[i-1])/(2*dx); /* */
    f[i]=wavespeed*p2+(x+wavespeed/mutation_expdecay-mutation_rate)*p3+(x/mutation_expdecay+1.)*c[i] - 2*p4-2*u[i]*c[i]/mutation_expdecay;
    fc[i]=-wavespeed*2./(dx*dx)+(x/mutation_expdecay+1.)-2*u[i]/mutation_expdecay;
    c[i]=c[i]-alpha*f[i]/fc[i];
  }
  for (i=1; i<space-2; i+=2) {	/* then odd-numbered items */
    x=(i-space0)*dx; 
    p2=(c[i+1]-2.*c[i]+c[i-1])/(dx*dx); /* 2nd derivative */
    p3=(c[i+1]-c[i-1])/(2*dx); /*1st drift derivative */
    p4=(c[i+1]*u[i+1]-c[i-1]*u[i-1])/(2*dx); /* */
    f[i]=wavespeed*p2+(x+wavespeed/mutation_expdecay-mutation_rate)*p3+(x/mutation_expdecay+1.)*c[i] - 2*p4-2*u[i]*c[i]/mutation_expdecay;
    fc[i]=-wavespeed*2./(dx*dx)+(x/mutation_expdecay+1.)-2*u[i]/mutation_expdecay;
    c[i]=c[i]-alpha*f[i]/fc[i];
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
  double x,y,tmp;
  
  // ==================================================
  // initilization for U mode
  // ==================================================
  
  if (switchu == 1) {
    if (have_u_infile == 1) {
      read_u_file(1);
    }else{
      u = (double*)malloc(space*sizeof(double));
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
  
  
  
  if(switchc == 1) {
    if( have_c_infile == 1) {
      read_c_file(1);
      read_u_file(0);
    }else{
      read_u_file(1);
      c = (double*)malloc(space*sizeof(double));
      switch(initialcond) {
	case 1:	initial_c_approx();
		break;
	case 2:	initial_c_random();
		break;
	case 3:	initial_c_flat();
		break;
	default:print_error("something went wrong, option for initial condition not available");
      }
      adjust_constraint(0);
    }
    
    f  = (double*)malloc(space*sizeof(double));
    fc = (double*)malloc(space*sizeof(double));
  }
}



void print_options() {
  printf("###################################################\n");
  printf("# solve profiles and weight-functions numerically #\n");
  printf("###################################################\n");  
  printf("#\n");
  printf("# mutationmodel     = %d\n",mutation_type);
  printf("# wavespeed         = %lf\n",wavespeed);
  printf("# mutationrate      = %lf\n",mutation_rate);
  if(mutation_type == 2) {
    printf("# option            = %d\n",mutation_jumpwidth);
  }
}

void print_u() {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %e\n",(i-space0)*dx,u[i]);
  }
}

void print_c() {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %e\n",(i-space0)*dx,c[i]);
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
  }else if(switchc == 1){
    free(c);
    free(fc);
    free(f);
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
      if(i%adjust_step==0)adjust_constraint(i);
    }
    adjust_constraint(i);
    if(have_c_outfile == 1)write_c_file(i);
    if(quiet<2)print_c();
  }
  cleanup();
  return 0;
}



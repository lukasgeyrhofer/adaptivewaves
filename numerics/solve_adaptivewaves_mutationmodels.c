// solve_adaptivewaves
// numerical solution of the differential equations of u* and c in CRBWs
// lukas.geyrhofer@ds.mpg.de

// using newton's method at each lattice point to find solution

// v0.1		120604	first version
// v0.2		120619  minor updates
// v0.3		120704	flag for periodic/reflecting bc


// v0.4 	120820	modified code from oasis-problem (S = a \delta(x) -r) to adaptive waves (S= x -r)
//			using reduced variables

// v0.5		120906	new fileformat for binary files
// 			include output for densities

// v0.6		120910	can now read previous configurations,
//			-u -c output
//			-U -C input
//			option -m MODE: (MODE&1>0 - calc u), (MODE&2>0 - calc c)


// v0.6.01	120917	option -N, take no parameters, read all of them from previous c-file


// v0.7		130131	different mutation models
//			

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



double r = 0.0;
double v = 0.0;
double diffusionconstant = 1.0;
double dx = 0.1;
double xmax = 500.;
int itersteps = 10000;
int verbose = 0;
int outputstep = 1000;
int startwithgaussian = 0;

double clparameter_v = -1.;
double clparameter_mutationrate = -1.;

int set_c_to_zero = 1;

int ccheck_interval = 1000;

int periodic_bc = 0;

int write_u_to_file = 0;
char ufilename[128];
int u_outputstep = 1;

int write_c_to_file = 0;
char cfilename[128];

int calc_u = 1;
int calc_c = 1;

double r_sqrt;

int space, space0;
double zeropos = 0.5;

int sizecontrainedby_space = 0;


double *c,*u,*unew,*cnew;

double *change_u,*change_u_ratio;

double *c_coeff_m1, *c_coeff_0, *c_coeff_p1, *c_coeff_j;
double *u_coeff_m1, *u_coeff_0, *u_coeff_p1, *u_coeff_j;
double *uder_coeff_0,*cder_const;

int randominitialcond = 0;
const gsl_rng_type* T;
const gsl_rng* rg; /* gsl, global generator */
unsigned long randseed = 0;

int uselatticeparameters = 0;
double clparameter_dx = 1e-4;
int clparameter_space = 1e4;
int clparameter_space0 = 4e3;



int initialcondition = 2;

int read_u_from_file = 0;
char uinfilename[128];
int u_space,u_space0;
double u_dx;

int read_c_from_file = 0;
char cinfilename[128];
int c_space,c_space0;
double c_dx;

int mutationtype = 1;	// 1 - diffusion
			//      -m diffusion
			// 2 - single (beneficial) jumps
			//      -m jumps
			// 3 - exponential kernel, only beneficial
			//      -m exp
			// 4 - kernel from external file
			//      -m file
double diffusionconstant;
int mutation_allowoptions = 0;
double mutationrate = 0.001;
double mutation_exprate = 0.001;
int mutation_jumpwidth = 1;
char mutation_kernelfile[128];
double *mutationkernel;
int space_m,space0_m;

// function pointers

typedef int (*iteration_fp)(int timestep);
iteration_fp iterate_u;
iteration_fp iterate_c;

// ===========================================================================================
// print_error
// ===========================================================================================

int print_error(char *msg) {
  // print an error-msg, then quit program
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


// ===========================================================================================
// parsecommandline
// ===========================================================================================

void parsecommandline(int argn, char *argv[]) {
  char c;
  int haveinfile = 0;
  while((c = getopt(argn, argv,"c:C:u:U:d:s:z:v:D:S:O:o:m:M:A:PI:i:")) != -1){
    switch(c) {

      // input & output of c and u* files, lowercase=output, uppercase=input
      case 'c': write_c_to_file = 1;
		strcpy(cfilename,optarg);
		break;
      case 'C':	read_c_from_file = 1;
		strcpy(cinfilename,optarg);
		break;
      case 'u': write_u_to_file =1;		// write the solution u(x) to a file (name as parameter to "-u")
		strcpy(ufilename,optarg);
		break;
      case 'U':	read_u_from_file = 1;
		strcpy(uinfilename,optarg);
		break;

      // define lattice
      case 'P':	uselatticeparameters = 1;
		break;
      case 'd':	if(uselatticeparameters == 1) {
		  clparameter_dx = atof(optarg);		// lattice constant (distance between two points in the discretized lattice)
		}else{
		  print_error("lattice parameters (options -d DELTAX -s SPACE -S SPACE0) can only be set when using the parameter flag (-P)");
		}
		break;
      case 's':	if(uselatticeparameters == 1) {
		  sizecontrainedby_space = 1;
		  clparameter_space = atoi(optarg);
		}else{
		  print_error("lattice parameters (options -d DELTAX -s SPACE -S SPACE0) can only be set when using the parameter flag (-P)");
		}
		break;
      case 'z':	if(uselatticeparameters == 1) {
		  clparameter_space0 = atoi(optarg);
		}else{
		  print_error("lattice parameters (options -d DELTAX -s SPACE -S SPACE0) can only be set when using the parameter flag (-P)");
		}
		break;
      
      case 'I':	if(read_c_from_file == 0) {
		  if(strcmp(optarg,"random") == 0) {
		    initialcondition = 1;
		  }else if(strcmp(optarg,"uexp") == 0) {
		    initialcondition = 2;
		  }else{
		    print_error("valid initialcondition values '-I random' or '-I uexp'");
		  }
		}else{
		  print_error("cannot use option -I for initialcondition together with already reading profile from c-file");
		}
		break;

      // other external parameters
      case 'v':	clparameter_v = atof(optarg);		// convection velocity
		break;
      case 'D': clparameter_mutationrate = atof(optarg);	// formerly diffusionconstant, hence option -D
		break;

      // option for the algorithm
      case 'S':	itersteps = atoi(optarg);	// number of iterations to solve c and u*
		break;
      case 'O':	outputstep = atoi(optarg);	// output of convergence criteria every "outputstep" steps (criteria usually intergral over u/c)
		break;
      case 'o': u_outputstep = atoi(optarg);	// only used with option "-u FILENAME"
		break;				// output only every U_OUTPUTSTEP lattice point for the u-file
						// useful when the numerical solution is on a finer grid than the stochastic simulation

      // options for the mutation model, first provide option -m to set the mutation type,
      // then use option -M to specify further details
      case 'm': if(strcmp(optarg,"diffusion") == 0) {
		  mutationtype = 1;
		}else if(strcmp(optarg,"jumps") == 0) {
		  mutationtype = 2;
		}else if(strcmp(optarg,"exp") == 0) {
		  mutationtype = 3;
		}else if(strcmp(optarg,"file") == 0) {
		  mutationtype = 4;
		}else{
		  print_error("could not find mutationtype option! valid strings for option -m: 'diffusion', 'jumps', 'exp', 'file'");
		}
		mutation_allowoptions = 1;
		break;
      case 'M':	if(mutation_allowoptions==1) {
		  switch(mutationtype) {
		    case 2:	mutation_jumpwidth = atoi(optarg);
				break;
		    case 3:	mutation_exprate = atof(optarg);
				break;
		    case 4:	strcpy(mutation_kernelfile,optarg);
				break;
		  }
		}
		break;

      case 'A': if(strcmp(optarg,"u") == 0) {
		  calc_u = 1;
		  calc_c = 0;
		}else if(strcmp(optarg,"c") == 0) {
		  calc_u = 0;
		  calc_c = 1;
		}else{
		  calc_u = 1;
		  calc_c = 1;
		}
		break;

      case 'i': ccheck_interval = atoi(optarg);
		break;

//       case 'Z': set_c_to_zero = 1;
// 		break;
		
      
    }
  }
  if((uselatticeparameters==0)&&(read_c_from_file==0)&&(read_u_from_file==0)) 
    print_error("simulation lattice has to be specified either by commandline parameters (flag -P, options -d DELTAX -s SPACE -S SPACE0) or by either continuing with calculations on u* or c solutions (options -C FILENAME or -U FILENAME to load intermediate solutions");
  randseed = time(NULL);
}


// ===========================================================================================
// get_gaussian
// ===========================================================================================

double get_gaussian(double x, double sigma2) {
  return exp(-x*x/sigma2);
}



// ===========================================================================================
// read_u
// ===========================================================================================

void read_u() {
  FILE *fpu;
  double ux;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(uinfilename,"rb");
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
      free(tmpd);
    }
    
    dx = u_dx;
    space = u_space;
    space0 = u_space0;
    
    u=(double*)malloc(space*sizeof(double));
    fread(&u[0],space,sizeof(double),fpu);
    
    fclose(fpu);
  }
}



// ===========================================================================================
// read_c
// ===========================================================================================

void read_c() {
  FILE *fpc;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int c_space,c_space0;
  double c_dx;
  
  fpc = fopen(cinfilename,"rb");
  if(fpc != NULL) {
    fread(&icount,1,sizeof(int),fpc);
    fread(&dcount,1,sizeof(int),fpc);
    
    fread(&c_dx,1,sizeof(double),fpc);
    fread(&c_space,1,sizeof(int),fpc);
    fread(&c_space0,1,sizeof(int),fpc);
    
    if(icount > 0) {
      tmpi = (int*)malloc(icount*sizeof(int));
      fread(&tmpi[0],icount,sizeof(int),fpc);
      free(tmpi);
    }
    if(dcount >= 4) {
      tmpd = (double*)malloc((dcount-2)*sizeof(double));
      fread(&tmpd[0],2,sizeof(double),fpc);
      fread(&mutationrate,1,sizeof(double),fpc);
      fread(&v,1,sizeof(double),fpc);
      if(dcount>4)fread(&tmpd[2],dcount-4,sizeof(double),fpc);
      free(tmpd);
    }else{
      print_error("reading c: configuration file does not contain necessary parameters");
    }

    space = c_space;
    space0 = c_space0;
    dx = c_dx;

    c=(double*)malloc(space*sizeof(double));
    fread(&c[0],space,sizeof(double),fpc);
    
    fclose(fpc);
  }
}



// ===========================================================================================
// read_mutationkernel
// ===========================================================================================
void read_mutationkernel() {
  int i;
  int icount,dcount;
  int *tmpi;
  double *tmpd;
  double dx_m;
  FILE *fpm;
  
  fpm = fopen(mutation_kernelfile,"rb");
  if(fpm!=NULL) {
    fread(&icount,1,sizeof(int),fpm);
    fread(&dcount,1,sizeof(int),fpm);
    
    fread(&dx_m,1,sizeof(double),fpm);
    fread(&space_m,1,sizeof(int),fpm);
    fread(&space0_m,1,sizeof(int),fpm);
    
    if(icount > 0) {
      tmpi = (int*)malloc(icount*sizeof(int));
      fread(&tmpi[0],icount,sizeof(int),fpm);
      free(tmpi);
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpm);
      free(tmpd);
    }
    if((space_m!=2*space+1)||(space0_m!=space))print_error("mutationkernel lattice from file does not match!");
    
    mutationkernel = (double*)malloc(space_m*sizeof(double));
    fread(&mutationkernel[0],space_m,sizeof(double),fpm);
    
    fclose(fpm);
  }
}
    


// ===========================================================================================
// iterate_u
// ===========================================================================================

int iterate_u_diffusion(int timestep) {
  int i,j;
  double f,fu;
  double u_int = 0.0;

  // periodic boundary condition, hence u[0] and u[space-1] (first and last points) have to be treated separately (to use correct indices)
  // by using the predefined coefficients "u_coeff_?" and not doing calculations every step, the speedup is something between 15-20%
  f = u_coeff_m1[0] * u[space-1] + u_coeff_0[0]*u[0] + u_coeff_p1[0]*u[1] - 2.*u[0]*u[0];
  fu = uder_coeff_0[0]-4.*u[0];
  unew[0] = u[0] - f/fu;
#pragma omp parallel for 
  for(i=1;i<space-1;i++) {
    f = u_coeff_m1[i]*u[i-1] + u_coeff_0[i]*u[i] + u_coeff_p1[i]*u[i+1]-2.*u[i]*u[i];
    fu = uder_coeff_0[i]-4.*u[i];
    unew[i] = u[i] - f/fu;
  }
#pragma omp barrier
  
  f = u_coeff_m1[space-1] * u[space-2] + u_coeff_0[space-1]*u[space-1] + u_coeff_p1[space-1]*u[0] - 2.*u[space-1]*u[space-1];
  fu = uder_coeff_0[space-1]-4.*u[space-1];
  unew[space-1] = u[space-1] - f/fu;
  
  // after updating unew, everything copied to u
  memcpy(&u[0],&unew[0],space*sizeof(double));
  
  
  // output integral over u to see convergence of algorithm
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)u_int += u[i];
    u_int *= dx;
    printf("%8d %.10e nan\n",timestep,u_int);
  }
  return 0;
}

int iterate_u_jumps(int timestep) {
  int i,j;
  double f,fu;
  double u_int = 0.;
  double alpha = 1.;
  double ratiomax = 0.;

  // reflecting boundaries
  // assuming space >> mutation_jumpwidth
  f = u_coeff_0[0] * u[0] + u_coeff_j[0]*u[mutation_jumpwidth] - 2.*u[0]*u[0];
  fu = uder_coeff_0[0] - 4.*u[0];
//   change_u[0] = f/(fu*u[0]);
  unew[0] = u[0] - alpha*f/fu;
  
  for(i=1;i<space-mutation_jumpwidth;i++) {
    f = 0.5*v/dx* u[i-1] + (-mutationrate + (i-space0)*dx- 2.*u[i])*u[i] - 0.5*v/dx*u[i+1] + mutationrate*u[i+mutation_jumpwidth];
    fu = uder_coeff_0[i] - 4.*u[i];
//     printf("timestep = %8d i =%8d fu =%14.5e f = %14.5e f/fu = %14.5e\n",timestep,i,fu,f,f/fu);
    unew[i] = u[i] - alpha*f/fu;
//     change_u[i] = f/(fu*u[i]);
  }
  
  for(i=space-mutation_jumpwidth;i<space-1;i++) {
    f = 0.5*v/dx* u[i-1] + (-mutationrate + (i-space0)*dx-2.*u[i])*u[i] - 0.5*v/dx*u[i+1];
    f += mutationrate*0.5*(i+mutation_jumpwidth)*dx; //cheating: assuming u(x) = x/2 far out...
    fu = uder_coeff_0[i] - 4.*u[i];
//     change_u[i] = f/(fu*u[i]);
    unew[i] = u[i] - alpha*f/fu;
  }
  
  f = u_coeff_0[space-1]*u[space-1] + u_coeff_j[space-1]*0.5*(space-1+mutation_jumpwidth)*dx - 2.*u[space-1]*u[space-1];
  fu = uder_coeff_0[space-1]-4.*u[space-1];
//   change_u[space-1] = f/(fu*u[space-1]);
  unew[space-1] = u[space-1] - alpha*f/fu;
  
//   for(i=1;i<space-1;i++) {
//     if (change_u[i] > ratiomax)ratiomax = change_u[i];
//   }
//   
//   if(ratiomax > 0.8)alpha = 0.8/ratiomax;
// //   printf("alpha = %14.6e\n",alpha);
// //   alpha = 0.01;
//   for(i=0;i<space;i++) {
//     unew[i] = (1.-alpha*change_u[i])*u[i];
//   }
  
  
  memcpy(&u[0],&unew[0],space*sizeof(double));
  
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)u_int += u[i];
    u_int *= dx;
    printf("%8d %.10e\n",timestep,u_int);
  }
  
//   for(i=0;i<space;i++) {
//     printf("%lf %lf\n",(i-space0)*dx,u[i]);
//   }
//   exit(1);
//   
  
  return 0;
}

// ===========================================================================================
// iterate_c
// ===========================================================================================

int iterate_c_diffusion(int timestep) {
  int i,j;
  double f,fc;
  double c_int = 0.0;

  
  // same as for u, first and last treated separately
  
  if(set_c_to_zero == 1) {
    cnew[0] = 0.;
    cnew[space-1] = 0.;
  }else{
    f = c_coeff_m1[0] * c[space-1] + c_coeff_0[0]*c[0] + c_coeff_p1[0]*c[1];
    cnew[0] = c[0] - f/cder_const[0];
    
    f = c_coeff_m1[space-1] * c[space-2] + c_coeff_0[space-1]*c[space-1] + c_coeff_p1[space-1]*c[0];
    cnew[space-1] = c[space-1] - f/cder_const[space-1];
  }

  for(i=1;i<space-1;i++) {
    f = c_coeff_m1[i]*c[i-1] + c_coeff_0[i]*c[i] + c_coeff_p1[i]*c[i+1];
    cnew[i] = c[i] - f/cder_const[i];
  }
  
  memcpy(&c[0],&cnew[0],space*sizeof(double));

  // output population size to see convergence
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)c_int += c[i];
    c_int *= dx;
    printf("%8d nan %.10e\n",timestep,c_int);
  }
  
  return 0;
}



int iterate_c_jumps(int timestep) {
  int i,j;
  double f,fc;
  double c_int = 0.0;
  
  if(set_c_to_zero == 1) {
    cnew[0] = 0.;
    cnew[space-1] = 0.;
  }else{
    f = c_coeff_0[0]*c[0] + 27.388; // blabla, c should be zero anyway...
    
  }
  for(i=1;i<mutation_jumpwidth;i++) {
    f = c_coeff_m1[i]*c[i-1] + c_coeff_0[i]*c[i] + c_coeff_p1[i]*c[i+1];
    cnew[i] = c[i] - f/cder_const[i];
  }
  for(i=mutation_jumpwidth;i<space;i++) {
    f = c_coeff_m1[i]*c[i-1] + c_coeff_0[i]*c[i] + c_coeff_p1[i]*c[i+1] + c_coeff_j[i] * c[i-mutation_jumpwidth];
    cnew[i] = c[i] - f/cder_const[i];
  }
  
  memcpy(&c[0],&cnew[0],space*sizeof(double));
  
  if(timestep%outputstep==0) {
    for(i=0;i<space;i++)c_int += c[i];
    c_int *= dx;
    printf("%8d %.10e\n",timestep,c_int);
  }
  return 0;
}
  


// ===========================================================================================
// print_output
// ===========================================================================================

void print_output() {
  int i,j;
  double x;
  
  double popsize = 0.0;
  if(calc_c == 1) {
    for(i=0;i<space;i++) {
      popsize += c[i];
    }
    popsize *= dx;
    
    fprintf(stderr,"# N = %lf\n",popsize);
  }
  
  
  if((calc_u == 1)&&(calc_c == 0)) {
    fprintf(stderr,"# x\t\t--\t\tu(x)\n");
  }else{
    fprintf(stderr,"# x\t\tc(x)\t\t\tu(x)\t\t\tc(x)*u(x)\n");
  }
  for(i=0;i<space;i++) {
    x = (i-space0)*dx;
    if((calc_u == 1)&&(calc_c == 0)) {	// no solution for c present with option "-U", so second column is just "nan"
					// (second column is still written, so that u(x) is ALWAYS in third column)
      fprintf(stderr,"%lf\t%.10e\t%lf\t%lf\t%lf\t%lf\n",x,u[i],u_coeff_m1[i],u_coeff_0[i],u_coeff_p1[i],u_coeff_j[i]);
    }else{						// second column now c, fourth u*c
      fprintf(stderr,"%lf\t%.10e\t%.10e\t%.10e\n",x,c[i],u[i],c[i]*u[i]);
    }
  }
}



// ===========================================================================================
// write_u
// ===========================================================================================

void write_u() {
  // write u(x) solution to an additional file
  // which can then be read into the stochastic simulation
  // to have an arbitrary constraint
  
  FILE* fp;
  int i,j;
  
  int p_int    = 0; // number of integer parameters
  int p_double = 0; // number of double parameters
  
  double *u_write;
  int space_write, space0_write;
  double dx_write;
  
  int offset = space0 % u_outputstep;
  
  space_write = space/u_outputstep;
  dx_write = dx*u_outputstep;
  u_write = (double*)malloc(space_write*sizeof(double));
  
  j=0;
  for(i=0;i<space;i++) {
    if(i%u_outputstep == offset) {
      u_write[j] = u[i];
      if(i==space0)space0_write = j;
      j++;
    }
  }
  
  fp = fopen(ufilename,"wb");
  if(fp!=NULL) {
    fwrite(&p_int,1,sizeof(int),fp);
    fwrite(&p_double,1,sizeof(int),fp);
    fwrite(&dx_write,1,sizeof(double),fp);
    fwrite(&space_write,1,sizeof(int),fp);
    fwrite(&space0_write,1,sizeof(int),fp);
    fwrite(&u_write[0],space,sizeof(double),fp);
  }
  fclose(fp);
  free(u_write);
}



// ===========================================================================================
// write_c
// ===========================================================================================
  
void write_c() {
  FILE* fp;
  int i,j;
  
  int p_int    = 0; // number of integer parameters
  int p_double = 4; // number of double parameters
  
  double tmp_rho = 0.0;
  double tmp_starttime = 0.0;
  
  double *c_write;
  int space_write, space0_write;
  double dx_write;
  
  int offset = space0 % u_outputstep;
  
  space_write = space/u_outputstep;
  dx_write = dx*u_outputstep;
  c_write = (double*)malloc(space_write*sizeof(double));
  
  j=0;
  for(i=0;i<space;i++) {
    if(i%u_outputstep == offset) {
      c_write[j] = c[i];
      if(i==space0)space0_write = j;
      j++;
    }
  }
  
  fp = fopen(cfilename,"wb");
  if(fp!=NULL) {
    fwrite(&p_int,1,sizeof(int),fp);
    fwrite(&p_double,1,sizeof(int),fp);
    
    fwrite(&dx_write,1,sizeof(double),fp);
    fwrite(&space_write,1,sizeof(int),fp);
    fwrite(&space0_write,1,sizeof(int),fp);
    
    fwrite(&tmp_rho,1,sizeof(double),fp);
    fwrite(&tmp_starttime,1,sizeof(double),fp);
    fwrite(&diffusionconstant,1,sizeof(double),fp);
    fwrite(&v,1,sizeof(double),fp);
    
    fwrite(&c_write[0],space,sizeof(double),fp);
  }
  fclose(fp);
  free(c_write);
}
  
// ===========================================================================================
// consistency_check
// ===========================================================================================
  
void consistency_check(int timestep) {
  
  // integral (uc)dx should be 1
  
  int i,j;
  double popsize = 0.;
  double sum = 0.;
  
  for(i=0;i<space;i++) {
    popsize += c[i];
    sum += u[i]*c[i];
  }
  popsize *= dx;
  sum *= dx;
  
  // ODE for c(x) is only linear, so doesnt give correct prefactor
  // hence rescale c(x) to comply constraint:
  for(i=0;i<space;i++) {
    c[i] /= sum;
  }
  popsize /= sum; // population size has to be rescaled, too
}
 


// ===========================================================================================
// set initial conditions when starting computation from scratch
// ===========================================================================================


void initialcond_u_random() {
  int i;
  double xmax = (space-space0)*dx;
  for(i=0;i<space0;i++)u[i] = 2e-4*(4.+gsl_rng_uniform(rg));
  for(i=space0;i<space;i++)u[i] = 0.2*(i-space0)*dx*(4.+gsl_rng_uniform(rg));
}

void initialcond_u_approx() {
  int i;
  for(i=0;i<=space0;i++)u[i] = 0.0001;
  for(i=space0+1;i<space;i++)u[i]=0.5*(i-space0)*dx;
}  

void initialcond_c_approx() {
  int i;
  for(i=0;i<space;i++)c[i] = u[i]*exp(-(i-space0)*dx*v/mutationrate);
  consistency_check(0);
}

void initialcond_c_random() {
  int i;
  for(i=0;i<space;i++)c[i] = 0.001*gsl_rng_uniform(rg);
  consistency_check(0);
}

void initialcond_c_gaussian() {
  int i;
  double s2 = space*space*dx*dx*.001953125; // FWHM (full-width-half-maximum) of gaussian is roughly 1/6 of simulation box, 3FWHM is roughly the whole shape ...
  for(i=0;i<space;i++)c[i] = get_gaussian((i-space0)*dx,s2);
  consistency_check(0);
}

// ===========================================================================================
// initialize
// ===========================================================================================

void initialize() {
  int i,j;
  double x;
  double u2int_inv;
  double sigma2;
  double N_expectednovelocity,rho;

  double Dl2i;
  double v2li;

  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  
  if(uselatticeparameters == 1) {
    dx     = clparameter_dx;
    space  = clparameter_space;
    space0 = clparameter_space0;
//     printf("dx = %lf, space=%d, space0=%d\n",dx,space,space0);
    mutationrate = clparameter_mutationrate;
    v = clparameter_v;
  }else{
    if((calc_c == 1)&&(read_c_from_file==1)) {
      read_c();
      dx = c_dx;
      space = c_space;
      space0 = c_space0;
    }else if((calc_u == 1)&&(read_u_from_file == 1)) {
      read_u();
      dx = u_dx;
      space = u_space;
      space0 = u_space0;
      mutationrate = clparameter_mutationrate;
      v = clparameter_v;
    }else{
      print_error("could not determine calculation lattice.");
    }
  }

  
  if((calc_u == 1)&&(read_u_from_file == 0)) {
    u = (double*)malloc(space*sizeof(double));
    switch(initialcondition) {
      case 1:	initialcond_u_random();
		break;
      default:	initialcond_u_approx();
		break;
    }
  }
  
  
  if((calc_c == 1)&&(read_c_from_file == 0)) {
    c = (double*)malloc(space*sizeof(double));
    switch(initialcondition) {
      case 1:	initialcond_c_random();
		break;
      case 2:	initialcond_c_approx();
		break;
      case 3:	initialcond_c_gaussian();
		break;
    }
  }

  cnew = (double*)malloc(space*sizeof(double));
  unew = (double*)malloc(space*sizeof(double));
  
  c_coeff_0  = (double*)malloc(space*sizeof(double));	// coefficients in the discretized ODE for c[i], c[i-1] and c[i+1]
  c_coeff_m1 = (double*)malloc(space*sizeof(double));
  c_coeff_p1 = (double*)malloc(space*sizeof(double));
  c_coeff_j  = (double*)malloc(space*sizeof(double));

  u_coeff_0  = (double*)malloc(space*sizeof(double));	// coefficients in the discretized ODE for u[i], u[i-1] and u[i+1]
  u_coeff_m1 = (double*)malloc(space*sizeof(double));
  u_coeff_p1 = (double*)malloc(space*sizeof(double));
  u_coeff_j  = (double*)malloc(space*sizeof(double));
  
  uder_coeff_0 = (double*)malloc(space*sizeof(double));	// ODE differentiated with respect to u, for newton iteration
  cder_const = (double*)malloc(space*sizeof(double));	//				      c

  switch(mutationtype) {
    case 1:	iterate_u = iterate_u_diffusion;
		iterate_c = iterate_c_diffusion;
		Dl2i = mutationrate/(dx*dx);
		v2li = 0.5*v/dx;
  
		for(i=0;i<space;i++) {
		  x = (i-space0)*dx;

		  u_coeff_m1[i] = Dl2i + v2li;
		  u_coeff_0[i] = -2.*Dl2i - r + x;
		  u_coeff_p1[i] = Dl2i - v2li;
		  
		  c_coeff_m1[i] = Dl2i - v2li;
		  c_coeff_0[i] = -2.*Dl2i - r + x;
		  c_coeff_p1[i] = Dl2i + v2li;
		  
		  uder_coeff_0[i] = -2.*Dl2i-r + x;
		  
		  cder_const[i] = uder_coeff_0[i];
    
		  // periodic_bc == 0 means reflecting boundary condition,
		  // no dependence on convection anymore, coefficient at boundaries is doubled or zero.
		  if(i==0) {
		    u_coeff_m1[i] = 0.;
		    u_coeff_p1[i] = 2.*Dl2i;
		    c_coeff_m1[i] = 0.;
		    c_coeff_p1[i] = 2.*Dl2i;
		  }
		  if(i==space-1) {
		    u_coeff_m1[i] = 2.*Dl2i;
		    u_coeff_p1[i] = 0.;
		    c_coeff_m1[i] = 2.*Dl2i;
		    c_coeff_p1[i] = 0.;
		  }
		}
		break;
    case 2:	iterate_u = iterate_u_jumps;
		iterate_c = iterate_c_jumps;
		change_u = (double*)malloc(space*sizeof(double));
		change_u_ratio = (double*)malloc(space*sizeof(double));
		for(i=0;i<space;i++) {
		  change_u[i] = 0.;
		  change_u_ratio[i] = 0.;
		}
		v2li = 0.5*v/dx;
		for(i=0;i<space;i++) {
		  x = (i-space0)*dx;
		  
		  u_coeff_m1[i] = v2li;
		  u_coeff_0[i]  = -mutationrate + x - r;
		  u_coeff_p1[i] = -v2li;
		  u_coeff_j[i]  = mutationrate;
		  
		  c_coeff_m1[i] = -v2li;
		  c_coeff_0[i]  = -mutationrate + x -r;
		  c_coeff_p1[i] = v2li;
		  c_coeff_j[i]  = mutationrate;
		  
		  uder_coeff_0[i] = -mutationrate + x -r;
		  cder_const[i]   = uder_coeff_0[i];
		  // boundaries will be written directly into the iteration
		}
		fprintf(stderr,"#options\n");
		fprintf(stderr,"#mutationrate = %lf\n",mutationrate);
		fprintf(stderr,"#wavespeed = %lf\n",v);
		fprintf(stderr,"#jumpwidth = %d\n",mutation_jumpwidth);
		break;
    case 3:	// dummy
		break;
    case 4:	// dummy
		break;
  }
}



// ===========================================================================================
// cleanup
// ===========================================================================================

void cleanup() {
  // free all the allocated memory again
  free(c);
  free(u);
  free(unew);
  free(cnew);
  free(c_coeff_0);
  free(c_coeff_m1);
  free(c_coeff_p1);
  free(u_coeff_0);
  free(u_coeff_m1);
  free(u_coeff_p1);
  free(uder_coeff_0);
  free(cder_const);
}






// ===========================================================================================
// MAIN
// ===========================================================================================

int main(int argn, char *argv[]) {
  int i,j;
  parsecommandline(argn,argv);
  initialize();
  
  if(calc_u == 1) {
    for(i=0;i<itersteps;i++) iterate_u(i);
  }
  if(write_u_to_file == 1) {
    write_u();
  }
  if(calc_c == 1) {
    printf("calc_c\n");
    for(i=0;i<space;i++) {
      c_coeff_0[i] -= 2.*u[i];
      cder_const[i] -= 2.*u[i];
    }
    printf("\n");
    for(i=0;i<itersteps;i++) {
      printf("i = %d\n",i);
      iterate_c(i);
      if(i%ccheck_interval == 0)consistency_check(i);
    }
    consistency_check(itersteps);
  }
  print_output();
  
  if(write_c_to_file == 1) {
    write_c();
  }
  cleanup();
  return 0;
}


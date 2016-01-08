#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// generic algorithm parameters
int maxsteps = 10000;
double alpha = 1.;

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
int space = 10000,space0 = 5000;


// densities, etc...
double *c,*cnew;
double *u,*unew;

double *c_coeff_p,*c_coeff_0,*c_coeff_m,*c_coeff_j;
double *u_coeff_p,*u_coeff_0,*u_coeff_m,*u_coeff_j;


// model parameters
double wavespeed = 0.0001;
int useselectiongradient = 1;
double fisherwavestep = 1e-3;


// mutations
int    mutation_type = 1; // 1 - diffusion, 2 - jumps, 3 - expdecay, 4 - external file
double mutation_rate = 0.00001;
double mutation_asymmetry = 1e-10;
int    mutation_jumpwidth = 100;
char   mutation_kernelfile[128];
double *mutation_kernel;
double mutation_expdecay = 0.001;
int    allowmutationoptions = 0;

// initial conditions
int initialcond = 1; // 1 - known approximations, 2 - random, 3 - flat, 4 - gaussian (only c, unavailable for u)


// eigenvalue
double evgamma = 0.0;
int usegamma = 0;


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



int parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"UCi:o:u:I:S:a:s:d:z:m:M:v:D:QqG:F:")) != -1){
    switch(c) {
      case 'C':	if(switchu == 1) {
		  print_error("cannot use both options, -C and -U");
		}else{
		  switchc = 1;
		}
		break;
      case 'U':	if(switchc == 1) {
		  print_error("cannot use both options, -C and -U");
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
		}else if(strcmp(optarg,"gaussian") == 0) {
		  if(switchc==1) {
		    initialcond = 4;
		  }else{
		    print_error("gaussian initial conditions only available for density c(x), use option -C");
		  }
		}else{
		  print_error("argument provided to option -I is not valid: choose from 'approx', 'random', 'flat' or 'gaussian' (only mode -C)");
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
      case 'G':	if (switchc == 1) {
		  usegamma = 1;
		  evgamma = atof(optarg);
		}else{
		  print_error("eigenvalue only available in -C mode");
		}
		break;
      case 'F':	useselectiongradient = 0;
		fisherwavestep = atof(optarg);
		break;
    }
  }
  if(switchc + switchu == 0) {
    print_error("have to use at least one of the options, -C and -U");
  }
  if(switchc == 1) {
    if(have_u_infile == 0)print_error("need u* file to calculate density profile c(x)");
  }
  return 0;
}



// ===========================================================================================
// initial conditions
// ===========================================================================================

double get_gaussian(double pos, double sigma2) {
  return exp(-pos*pos/(2.*sigma2));
}

void initial_c_gaussian() {
  int i;
  double sigma2 = wavespeed;
  for(i=0;i<space;i++)c[i] = get_gaussian((i-space0)*dx,sigma2);
}

void initial_c_approx() {
  int i;
  if(mutation_type == 1) { // this approximation is only valid for diffusion type mutations ...
    for(i=0;i<space;i++)c[i] = exp(-wavespeed * (i-space0) *dx / mutation_rate)*u[i];
  }else{
    initial_c_gaussian();
  }
}

void initial_c_random() {
  int i;
  for(i=0;i<space;i++)c[i] = (dx*rand())/(1.*RAND_MAX);
}

void initial_c_flat() {
  int i;
  for(i=0;i<space;i++)c[i] = dx;
}

void initial_u_approx() {
  int i;
  double x;
  if(useselectiongradient == 1) {
    for(i=0;i<=space0;i++)u[i] = dx*exp(wavespeed*(i-space0)*dx/mutation_rate);
    for(i=space0+1;i<space;i++)u[i] = 0.5*(i-space0)*dx;
  }else{
    for(i=0;i<=space0;i++)u[i] = fisherwavestep*exp(wavespeed*(i-space0)*dx/mutation_rate);
    for(i=space0+1;i<space;i++)u[i] = fisherwavestep;
  }
}

void initial_u_flat() {
  int i;
  for(i=0;i<space;i++)u[i] = dx;
}

void initial_u_random() {
  int i;
  for(i=0;i<space;i++)u[i] = (dx*rand())/(1.*RAND_MAX);
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
      free(tmpd);
    }
    if(extractparameters == 1) {
      dx = u_dx;
      space = u_space;
      space0 = u_space0;
    }else if (( space != u_space) || (space0 != u_space0)) {
      print_error("lattice does not match");
    }
    u=(double*)malloc(space*sizeof(double));
    fread(&u[0],space,sizeof(double),fpu);
    fclose(fpu);
  }
}
  
  
void read_c_file(int extractparameters) {
  FILE *fpc;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int c_space,c_space0;
  double c_dx;
  
  fpc = fopen(c_infile,"rb");
  if(fpc != NULL) {
    fread(&icount,1,sizeof(int),fpc);
    fread(&dcount,1,sizeof(int),fpc);
    fread(&c_dx,1,sizeof(double),fpc);
    fread(&c_space,1,sizeof(int),fpc);
    fread(&c_space0,1,sizeof(int),fpc);
    if(icount > 0) {
      if(extractparameters==1) {
	fread(&mutation_type,1,sizeof(int),fpc);
	if(icount>=2) {
	  switch(mutation_type) {
	    case 2:	fread(&mutation_jumpwidth,1,sizeof(int),fpc);
			break;
	    default:	if(icount>=2) {
			  tmpi = (int*)malloc((icount-1)*sizeof(int));
			  fread(&tmpi[0],icount-1,sizeof(int),fpc);
			  free(tmpi);
			}
			break;
	  }
	}
      }else{
	tmpi = (int*)malloc(icount*sizeof(int));
	fread(&tmpi[0],icount,sizeof(int),fpc);
	free(tmpi);
      }
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpc);
      free(tmpd);
    }
    if(extractparameters == 1) {
      dx = c_dx;
      space = c_space;
      space0 = c_space0;
    }else if (( space != c_space) || (space0 != c_space0)) {
      print_error("lattice does not match");
    }
    c=(double*)malloc(space*sizeof(double));
    fread(&c[0],space,sizeof(double),fpc);
    fclose(fpc);
  }
}




// ===========================================================================================
// write output files
// ===========================================================================================

void write_c_file() {
  FILE *fpc;
  int icount = 1,dcount = 4;
  double tmp_starttime = 0.,tmp_rho = 0.;
  fpc = fopen(c_outfile,"wb");
  if(mutation_type==2)icount=2;
  fwrite(&icount,1,sizeof(int),fpc);
  fwrite(&dcount,1,sizeof(int),fpc);
  fwrite(&dx,1,sizeof(double),fpc);
  fwrite(&space,1,sizeof(int),fpc);
  fwrite(&space0,1,sizeof(int),fpc);
  fwrite(&mutation_type,1,sizeof(int),fpc);
  if(mutation_type==2)fwrite(&mutation_jumpwidth,1,sizeof(int),fpc);
  fwrite(&tmp_rho,1,sizeof(double),fpc);
  fwrite(&tmp_starttime,1,sizeof(double),fpc);
  fwrite(&mutation_rate,1,sizeof(double),fpc);
  fwrite(&wavespeed,1,sizeof(double),fpc);
  fwrite(&c[0],space,sizeof(double),fpc);
  fclose(fpc);
}

void write_u_file() {
  FILE *fpu;
  int icount = 0,dcount = 3-useselectiongradient;
  fpu = fopen(u_outfile,"wb");
  fwrite(&icount,1,sizeof(int),fpu);
  fwrite(&dcount,1,sizeof(int),fpu);
  fwrite(&dx,1,sizeof(double),fpu);
  fwrite(&space,1,sizeof(int),fpu);
  fwrite(&space0,1,sizeof(int),fpu);
  fwrite(&mutation_rate,sizeof(double),1,fpu);
  fwrite(&wavespeed,sizeof(double),1,fpu);
  if(useselectiongradient == 0) {
    fwrite(&fisherwavestep,sizeof(double),1,fpu);
  }
  fwrite(&u[0],space,sizeof(double),fpu);
  fclose(fpu);
}



// ===========================================================================================
// iteration procedures
// ===========================================================================================

int iterate_c_diffusion(int timestep) {
  int i;
  double f,fc;
  
  //assume exponential decay on both sides...
  f = c_coeff_m[0]*(c[0]*c[0]/c[1]) + c_coeff_0[0]*c[0] + c_coeff_p[0]*c[1];
  fc = c_coeff_0[0];
  cnew[0] = c[0] - alpha * f/fc;

  for(i=1;i<space-1;i++) {
    f = c_coeff_m[i] * c[i-1]  + c_coeff_0[i]* c[i] + c_coeff_p[i]*c[i+1];
    fc = c_coeff_0[i];
    cnew[i] = c[i] - alpha * f/fc;
  }

  f = c_coeff_m[space-1]*c[space-2] + c_coeff_0[space-1]*c[space-1] + c_coeff_p[space-1] * (c[space-1]*c[space-1]/c[space-2]);
  fc = c_coeff_0[space-1];
  cnew[space-1] = c[space-1] - alpha * f/fc;
  
  memcpy(&c[0],&cnew[0],space*sizeof(double));
}


int iterate_c_jumps(int timestep) {
  int i;
  double f,fc;
  
  // ignore indices 0 and (space-1), such that c[0] = 0 and c[space-1] = 0
  for(i=1;i<mutation_jumpwidth;i++) {
    f = c_coeff_m[i] * c[i-1] + c_coeff_0[i] * c[i] + c_coeff_p[i] * c[i+1];
    fc = c_coeff_0[i];
    fprintf(stderr,"%d %lf %lf %lf %lf\n",timestep,(i-space0)*dx,f,fc,f/fc);
    cnew[i] = c[i] - alpha*f/fc;
  }
  for(i=mutation_jumpwidth;i<space-1;i++) {
    f = c_coeff_m[i] * c[i-1] + c_coeff_0[i] * c[i] + c_coeff_p[i] * c[i+1] + c_coeff_j[i] * c[i-mutation_jumpwidth];
    fc = c_coeff_0[i];
    fprintf(stderr,"%d %lf %lf %lf %lf\n",timestep,(i-space0)*dx,f,fc,f/fc);
    cnew[i] = c[i] - alpha*f/fc;
  }
  fprintf(stderr,"\n");
  memcpy(&c[0],&cnew[0],space*sizeof(double));
}

int iterate_u_diffusion(int timestep) {
  int i;
  
  double f,fu;
  
  
  // u[-1] = exp decay from u[0]
  f =  u_coeff_m[0]*u[0]*u[0]/u[1] + u_coeff_0[0] * u[0] + u_coeff_p[0] * u[1] - 2.*u[0]*u[0];
  fu = u_coeff_0[0] - 4.*u[0];
  unew[0] = u[0] - alpha * f/fu;
  for(i=1;i<space-1;i++) {
    f = u_coeff_m[i]*u[i-1] + u_coeff_0[i]*u[i] + u_coeff_p[i]*u[i+1] - 2.*u[i]*u[i];
    fu = u_coeff_0[i] - 4.*u[i];
    unew[i] = u[i] - alpha * f/fu;
  }
  // u[space] = linear extrapolation from u[space-1]
  f = u_coeff_m[space-1] * u[space-2] + u_coeff_0[space-1]*u[space-1] + u_coeff_p[space-1] * (2*u[space-1]-u[space-2]) - 2.*u[space-1]*u[space-1];
  fu = u_coeff_0[space-1] - 4.*u[space-1];
  unew[space-1] = u[space-1] - alpha * f/fu;
  
  memcpy(&u[0],&unew[0],space*sizeof(double));
}

int iterate_u_jumps(int timestep) {
  int i;
  double f,fu;
  
//   u[0] = 0.;
  for(i=1;i<space-mutation_jumpwidth;i++) {
    // u ~ 0 for small x
    f = u_coeff_m[i]*u[i-1] + u_coeff_0[i]*u[i] + u_coeff_p[i]*u[i+1] + u_coeff_j[i]*u[i+mutation_jumpwidth] - 2.*u[i]*u[i];
    fu = u_coeff_0[i] - 4.*u[i];
    unew[i] = u[i] - alpha*f/fu;
  }
  for(i=space-mutation_jumpwidth;i<space-1;i++) {
    f = u_coeff_m[i]*u[i-1] + u_coeff_0[i]*u[i] + u_coeff_p[i]*u[i+1] + u_coeff_j[i]*(i-space0)*dx*0.5 - 2.*u[i]*u[i]; 
    fu = u_coeff_0[i] - 4.*u[i];
    unew[i] = u[i] - alpha*f/fu;
  }
//   u[space-1] = (space-1-space0)*dx*.5;
  memcpy(&u[0],&unew[0],space*sizeof(double));  
}




// ===========================================================================================
// initialize all variables
// ===========================================================================================

void set_c_bc() {
  c[0] = 0.;
  c[space-1] = 0.;
}

void set_u_bc() {
  u[0] = 0.;
  unew[0] = 0.;
  u[space-1] = (space-1-space0)*dx*0.5;
  unew[space-1] = u[space-1];
}

void initialize() {
  int i;
  double x;
  
  // ==================================================
  // initialization for C mode
  // ==================================================
  
  if (switchc == 1) {
    if(have_c_infile == 1) {
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
	case 4:	initial_c_gaussian();
		break;
	default:print_error("something went wrong, option for initial condition not available");
      }
    }
    cnew = (double*)malloc(space*sizeof(double));
    
    c_coeff_m = (double*)malloc(space*sizeof(double));
    c_coeff_0 = (double*)malloc(space*sizeof(double));
    c_coeff_p = (double*)malloc(space*sizeof(double));
    
    if(mutation_type == 2)c_coeff_j = (double*)malloc(space*sizeof(double));
    
    
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      
      c_coeff_m[i] = -0.5*wavespeed/dx;
      if(useselectiongradient == 1) {
	c_coeff_0[i] = x - 2.*u[i];
      }else{
	c_coeff_0[i] = (i>=space0)*fisherwavestep - 2*u[i];
      }
      c_coeff_p[i] = 0.5*wavespeed/dx;
      
      if(usegamma == 1) {
	c_coeff_0[i] -= evgamma;
      }
      
      switch(mutation_type) {
	case 1:	c_coeff_m[i] += mutation_rate/(dx*dx);
		c_coeff_0[i] += -2.*mutation_rate/(dx*dx);
		c_coeff_p[i] += mutation_rate/(dx*dx);
		iterate_c = iterate_c_diffusion;
		break;
	case 2:	c_coeff_0[i] += -mutation_rate;
		c_coeff_j[i]  = mutation_rate;
		iterate_c = iterate_c_jumps;
		break;
	default:print_error("other mutation models not yet implemented");
      }
    }
    
    set_c_bc();
    
  }

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
	default:print_error("something went wrong, option for initial condition not available");
      }
    }
    unew = (double*)malloc(space*sizeof(double));
    
    u_coeff_m = (double*)malloc(space*sizeof(double));
    u_coeff_0 = (double*)malloc(space*sizeof(double));
    u_coeff_p = (double*)malloc(space*sizeof(double));
    u_coeff_j = (double*)malloc(space*sizeof(double));
    
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      u_coeff_m[i] = 0.5*wavespeed/dx;
      if(useselectiongradient == 1) {
	u_coeff_0[i] = x;
      }else{
	u_coeff_0[i] = (i>=space0)*fisherwavestep;
      }
      u_coeff_p[i] = -0.5*wavespeed/dx;
      
      switch(mutation_type) {
	case 1:	u_coeff_m[i] += mutation_rate/(dx*dx);
		u_coeff_0[i] += -2.*mutation_rate/(dx*dx);
		u_coeff_p[i] += mutation_rate/(dx*dx);
		iterate_u = iterate_u_diffusion;
		break;
	case 2:	u_coeff_0[i] += -mutation_rate;
		u_coeff_j[i]  =  mutation_rate;
		iterate_u = iterate_u_jumps;
		break;
	default:print_error("other mutation models not yet implemented");
      }
    }
    set_u_bc();
  }
}



void print_options() {
  printf("###################################################\n");
  printf("# solve profiles and weight-functions numerically #\n");
  printf("###################################################\n");  
  printf("#\n");
  printf("# mutationmodel     = %d\n",mutation_type);
  if(mutation_type != 1) {
    printf("# option            = %d\n",mutation_jumpwidth);
  }
  printf("# mutation_rate     = %.6e\n",mutation_rate);
  printf("# wavespeed         = %.6e\n",wavespeed);
  if(useselectiongradient == 1) {
    printf("# selection         = gradient\n");
  }else{
    printf("# selection         = step (Fisher waves)\n");
    printf("# wavestep          = %.6e\n",fisherwavestep);
  }
}

void print_u() {
  int i;
  for(i=0;i<space;i++) {
    printf("%e %e\n",(i-space0)*dx,u[i]);
  }
}

void print_c(int timestep) {
  int i;
  for(i=0;i<space;i++) {
    printf("%d %lf %e\n",timestep,(i-space0)*dx,c[i]);
  }
  printf("\n");
}

void print_coeffs() {
  int i;
  
  for(i=0;i<space;i++) {
    printf("%lf %lf %lf %lf",(i-space0)*dx,c_coeff_m[i],c_coeff_0[i],c_coeff_p[i]);
    if(mutation_type == 2) printf(" %lf",c_coeff_j[i]);
    printf("\n");
  }
}


void print_coeff_matrix() {
  int i,j;
  double x;
  for(i=0;i<space;i++) {
    for(j=0;j<space;j++) {
      if(i==j) {
	printf(" %lf",c_coeff_0[i]);
      }else if(j==i-1) {
	printf(" %lf",c_coeff_m[i]);
      }else if(j==i+1) {
	printf(" %lf",c_coeff_p[i]);
      }else if(j==i-mutation_jumpwidth) {
	printf(" %lf",c_coeff_j[i]);
      }else{
	printf(" 0.");
      }
    }
    printf("\n");
  }
}
      

// ===========================================================================================
// adjust the density to fulfil the constraint \int dx u(x) c(x) = 1
// ===========================================================================================

void adjust_constraint() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++)sum += u[i]*c[i];
  sum *= dx;
  for(i=0;i<space;i++)c[i] /= sum;
}

// ===========================================================================================
// cleanup
// ===========================================================================================

void cleanup() {
  int i;
  if(switchc == 1) {
    free(c);
    free(cnew);
    free(c_coeff_m);
    free(c_coeff_0);
    free(c_coeff_p);
    if(mutation_type == 2)free(c_coeff_j);
    free(u);
  }else{
    free(u);
    free(unew);
    free(u_coeff_m);
    free(u_coeff_0);
    free(u_coeff_p);
    if(mutation_type == 2)free(u_coeff_j);
  }
  if(mutation_type == 4) {
    free(mutation_kernel);
  }
}



// ===========================================================================================
// main
// ===========================================================================================

int main(int argn, char *argv[]) {
  int i;
  
  parsecommandline(argn,argv);
  initialize();
  
  if(quiet<2)print_options();
  if(switchu == 1) {
    for(i=1;i<=maxsteps;i++) {
      iterate_u(i);
    }
    if(quiet==0)print_u();
    if(have_u_outfile == 1)write_u_file();
  }else if(switchc == 1) {
    adjust_constraint();
    for(i=1;i<=maxsteps;i++) {
      iterate_c(i);
      adjust_constraint();
    }
    adjust_constraint();
    if(quiet==0)print_c(i);
    if(have_c_outfile == 1)write_c_file();
  }
  cleanup();
  return 0;
}



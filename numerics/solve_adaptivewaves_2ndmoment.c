/*
solve for 2point-function of adaptivewaves problems

v0.01	somewhen mar/apr13
	crude implementation of 2nd moment, basically just copy&paste from 1d problem

v0.02	130424
	fixing bugs, improve performance:
	  * calculate only upper half of xy-plane, problem is symmetric!

v0.03	130430
	implemented numerical solution for left eigenvector PHI^(2)
v0.03.1	130502
	rescale PHI^(2) with C^(2), < PHI^(2) | C^(2) > == 1
	(optional)

*/
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



// switchc == 1 -> calculate density c(x)
// switchu == 1 -> calculate u*(x)
int switchc = 0;
int switchu = 0;

// names and flags for binary input- and output-file
char c_infile[128], c_outfile[128];
char u_infile[128], u_outfile[128];
char phi2_infile[128], phi2_outfile[128];
char c2_infile[128];
int have_c_infile = 0, have_c_outfile = 0;
int have_u_infile = 0, have_u_outfile = 0;
int have_phi_infile = 0, have_phi_outfile = 0;
int have_c2_infile = 0;

// lattice
double dx = 1e-4;
int space = 10000,space0 = 5000;


// densities, etc...
double *c,*cnew;
double *u,*unew;
double *u4;

double **c2,**c2new;

double **phi2,**phi2new;


double *c_coeff_p,*c_coeff_0,*c_coeff_m,*c_coeff_j;
double *u_coeff_p,*u_coeff_0,*u_coeff_m,*u_coeff_j;


// model parameters
double wavespeed = 0.0001;

int adjustcontraintsteps = 10000;


// mutations
int    mutation_type = 1; // 1 - diffusion, 2 - jumps, 3 - expdecay, 4 - external file
double mutation_rate = 0.00001;
int    mutation_jumpwidth = 100;
char   mutation_kernelfile[128];
double *mutation_kernel;
double mutation_expdecay = 0.001;
int    allowmutationoptions = 0;

// initial conditions
int initialcond = 1; // 1 - known approximations, 2 - random, 3 - flat, 4 - gaussian (only c, unavailable for u)


// function pointers
typedef int (*iteration_fp)(int timestep);
iteration_fp iterate_u;
iteration_fp iterate_c;


// output
int verbosity = 2;

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
  while((c = getopt(argn, argv,"UCi:o:u:c:I:S:a:s:d:z:m:M:v:D:A:Qq")) != -1){
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
		  strcpy(phi2_infile,optarg);
		  have_phi_infile = 1;
		}else{
		  print_error("option -C or -U has to be specified first");
		}
		break;
      case 'o':	if(switchc == 1) {
		  strcpy(c_outfile,optarg);
		  have_c_outfile = 1;
		}else if(switchu == 1) {
		  strcpy(phi2_outfile,optarg);
		  have_phi_outfile = 1;
		}else{
		  print_error("option -C or -U has to be specified first");
		}
		break;
      case 'u':	strcpy(u_infile,optarg);
		have_u_infile = 1;
		break;
      case 'c':	if(switchu == 1) {
		  strcpy(c_infile,optarg);
		  have_c2_infile = 1;
		}else{
		  print_error("option -c only specifies constraint < PHI^(2) | C^(2) > == 1 in mode -U");
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
      case 'A': adjustcontraintsteps = atoi(optarg);
		break;
      case 'Q':	verbosity = 1;
		break;
      case 'q':	verbosity = 0;
		break;
    }
  }
  if(switchc + switchu == 0) {
    print_error("have to use at least one of the options, -C and -U");
  }
  if(switchc == 1) {
    if(have_u_infile == 0)print_error("need u* file to calculate density profile C^(2)");
  }
  if(switchu == 1) {
    if(have_u_infile == 0)print_error("need u* file to calculate PHI^(2)");
//     if(have_c2_infile == 0)print_error("need C^(2) file to normalize PHI^(2), option -c");
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
  int i,j;
  double x,y;
  if(mutation_type == 1) { // this approximation is only valid for diffusion type mutations ...
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      for(j=0;j<i;j++) {
	y = (j-space0)*dx;
	c2[i][j] = exp(-wavespeed * x / mutation_rate - wavespeed*y/mutation_rate)*u[i]*u[j];
      }
    }
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
  for(i=0;i<=space0;i++)u[i] = 0.;
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
  
  
void read_c2_file(int extractparameters) {
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
    if(icount >= 3) {
      if(extractparameters==1) {
	tmpi = (int*)malloc(icount*sizeof(int));
	fread(&tmpi[0],icount,sizeof(int),fpc);
	mutation_type = tmpi[0];
	if(mutation_type == 2)mutation_jumpwidth = tmpi[1];
	if(tmpi[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
	free(tmpi);
      }else{
	tmpi = (int*)malloc(icount*sizeof(int));
	fread(&tmpi[0],icount,sizeof(int),fpc);
	free(tmpi);
      }
    }else{
      print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
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
    c2=(double**)malloc(space*sizeof(double*));
    tmpd = (double*)malloc(space*sizeof(double));
    for(i=0;i<space;i++) {
      c2[i] = (double*)malloc((i+1)*sizeof(double));
      fread(&c2[i][0],i+1,sizeof(double),fpc);
      if(space-i-1>0)fread(&tmpd[0],space-i-1,sizeof(double),fpc);
    }
    free(tmpd);
    fclose(fpc);
  }
}



void read_phi2_file(int extractparameters) {
  FILE *fpp;
  int i,icount,dcount;
  int *tmpi;
  double *tmpd;
  int p_space,p_space0;
  double p_dx;
  
  fpp = fopen(phi2_infile,"rb");
  if(fpp!=NULL) {
    fread(&icount,1,sizeof(int),fpp);
    fread(&dcount,1,sizeof(int),fpp);
    fread(&p_dx,1,sizeof(double),fpp);
    fread(&p_space,1,sizeof(int),fpp);
    fread(&p_space0,1,sizeof(int),fpp);
    if(icount>=3) {
      if(extractparameters==1) {
	tmpi=(int*)malloc(icount*sizeof(int));
	fread(&tmpi[0],icount,sizeof(int),fpp);
	mutation_type = tmpi[0];
	if(mutation_type == 2)mutation_jumpwidth = tmpi[1];
	if(tmpi[2] != 1)print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
	free(tmpi);
      }else{
	tmpi=(int*)malloc(icount*sizeof(int));
	fread(&tmpi[0],icount,sizeof(int),fpp);
	free(tmpi);
      }
    }else{
      print_error("wrong parameters in c-file. 3rd int-parameter has to be set to '1', in order for the file to be 2d.");
    }
    if(dcount > 0) {
      tmpd = (double*)malloc(dcount*sizeof(double));
      fread(&tmpd[0],dcount,sizeof(double),fpp);
      free(tmpd);
    }
    if(extractparameters == 1) {
      dx = p_dx;
      space = p_space;
      space0 = p_space0;
    }else if (( space != p_space) || (space0 != p_space0)) {
      print_error("lattice does not match");
    }
    phi2 = (double**)malloc(space*sizeof(double*));
    phi2new = (double**)malloc(space*sizeof(double*));
    tmpd = (double*)malloc(space*sizeof(double));
    for(i=0;i<space;i++) {
      phi2[i] = (double*)malloc((i+1)*sizeof(double));
      phi2new[i] = (double*)malloc((i+1)*sizeof(double));
      fread(&phi2[i][0],i+1,sizeof(double),fpp);
      if(space-i-1>0)fread(&tmpd[0],space-i-1,sizeof(double),fpp);
    }
    free(tmpd);
    fclose(fpp);
  }
}
    
  


// ===========================================================================================
// write output files
// ===========================================================================================

void write_c2_file() {
  FILE *fpc;
  int i,j;
  int icount = 3,dcount = 4;
  int tmpi_zeros = 0,twodim=1;
  double tmp_starttime = 0.,tmp_rho = 0.;
  fpc = fopen(c_outfile,"wb");
  //if(mutation_type==2)icount=2;
  fwrite(&icount,1,sizeof(int),fpc);
  fwrite(&dcount,1,sizeof(int),fpc);
  fwrite(&dx,1,sizeof(double),fpc);
  fwrite(&space,1,sizeof(int),fpc);
  fwrite(&space0,1,sizeof(int),fpc);
  fwrite(&mutation_type,1,sizeof(int),fpc);
  if(mutation_type==2) {
    fwrite(&mutation_jumpwidth,1,sizeof(int),fpc);
  }else{
    fwrite(&tmpi_zeros,1,sizeof(int),fpc);
  }
  fwrite(&twodim,1,sizeof(int),fpc);
  fwrite(&tmp_rho,1,sizeof(double),fpc);
  fwrite(&tmp_starttime,1,sizeof(double),fpc);
  fwrite(&mutation_rate,1,sizeof(double),fpc);
  fwrite(&wavespeed,1,sizeof(double),fpc);
  for(i=0;i<space;i++) {
    fwrite(&c2[i][0],i+1,sizeof(double),fpc);
    for(j=i+1;j<space;j++) {
      fwrite(&c2[j][i],1,sizeof(double),fpc);
    }
  }
  fclose(fpc);
}


void write_phi2_file() {
  FILE *fpp;
  int i,j;
  int icount = 3,dcount = 0;
  int tmpi_zeros = 0,twodim = 1;
  double tmp_starttime = 0., tmp_rho = 0.;
  fpp = fopen(phi2_outfile,"wb");
  fwrite(&icount,1,sizeof(int),fpp);
  fwrite(&dcount,1,sizeof(int),fpp);
  fwrite(&dx,1,sizeof(double),fpp);
  fwrite(&space,1,sizeof(int),fpp);
  fwrite(&space0,1,sizeof(int),fpp);
  fwrite(&mutation_type,1,sizeof(int),fpp);
  if(mutation_type==2) {
    fwrite(&mutation_jumpwidth,1,sizeof(int),fpp);
  }else{
    fwrite(&tmpi_zeros,1,sizeof(int),fpp);
  }
  fwrite(&twodim,1,sizeof(int),fpp);
  for(i=0;i<space;i++) {
    fwrite(&phi2[i][0],i+1,sizeof(double),fpp);
    for(j=i+1;j<space;j++) {
      fwrite(&phi2[j][i],1,sizeof(double),fpp);
    }
  }
  fclose(fpp);
}


// ===========================================================================================
// iteration procedures
// ===========================================================================================


int iterate_c_diffusion(int timestep) {
  int i,j;
  double f,fc;
  
  for(i=1;i<space-1;i++) {
    c[i] = 0.;
    
    // offdiagonal elements
    for(j=1;j<i;j++) {
      f  = c_coeff_m[i]*c2[i-1][j] + c_coeff_0[i]*c2[i][j] + c_coeff_p[i]*c2[i+1][j];	// first dimension (x)
      f += c_coeff_m[j]*c2[i][j-1] + c_coeff_0[j]*c2[i][j] + c_coeff_p[j]*c2[i][j+1];	// second dimension (y)
      fc = c_coeff_0[i] + c_coeff_0[j];
      
      c2new[i][j] = c2[i][j] - alpha*f/fc;
      
      c[i] += c2[i][j]*u[j];	// integrate c2 to c1 (at least first part)
    }
    
    for(j=i;j<space;j++)c[i] += c2[j][i]*u[j];
				// finalize integration from c2 to c1, save i iterations with this split
//     c[i] *= dx;
    
    // i == j, diagonal elements
    f  = 2.*(c_coeff_m[i]*c2[i][i-1] + c_coeff_0[i]*c2[i][i] + c_coeff_p[i]*c2[i+1][i] + c[i]);
    fc = 2.*c_coeff_0[i] - 2*u[i];
    c2new[i][i] = c2[i][i] - alpha*f/fc;

  }
  
  for(i=1;i<space-1;i++)memcpy(&c2[i][1],&c2new[i][1],i*sizeof(double));
}

int iterate_c_jumps(int timestep) {
  int i;
}

int iterate_u_diffusion(int timestep) {
  int i,j;
  
  double f,fu;
  
  for(i=1;i<space-1;i++) {
    for(j=1;j<i;j++) {
      f  = u_coeff_m[i]*phi2[i-1][j] + u_coeff_0[i]*phi2[i][j] + u_coeff_p[i]*phi2[i+1][j];
      f += u_coeff_m[j]*phi2[i][j-1] + u_coeff_0[j]*phi2[i][j] + u_coeff_p[j]*phi2[i][j+1];
      f += u[i]*phi2[j][j] + u[j]*phi2[i][i];
      fu = u_coeff_0[i]+u_coeff_0[j];
      phi2new[i][j] = phi2[i][j] - alpha*f/fu;
    }
    f  = 2.*(u_coeff_m[i]*phi2[i][i-1] + u_coeff_0[i]*phi2[i][i] + u_coeff_p[i]*phi2[i+1][i] + u[i]*phi2[i][i]);
    fu = 2.*u_coeff_0[i] + 2.*u[i];
    phi2new[i][i] = phi2[i][i] - alpha*f/fu;
  }
  for(i=1;i<space-1;i++)memcpy(&phi2[i][1],&phi2new[i][1],i*sizeof(double));
}

int iterate_u_jumps(int timestep) {
  int i;
  double f,fu;
  
  f = (-mutation_rate + (0-space0)*dx)*u[0] + - 0.5*wavespeed/dx*u[1] + mutation_rate*u[mutation_jumpwidth]- 2.*u[0]*u[0];
  fu = (-mutation_rate + (0-space0)*dx) - 4.*u[0];
  unew[0] = u[0] - alpha * f/fu;
  for(i=1;i<=space-1-mutation_jumpwidth;i++) {
    f = 0.5*wavespeed/dx*u[i-1] + (-mutation_rate + (i-space0)*dx)*u[i] - 0.5*wavespeed/dx*u[i+1] + mutation_rate*u[i+mutation_jumpwidth] - 2.*u[i]*u[i];
    fu = (-mutation_rate + (i-space0)*dx) - 4.*u[i];
    unew[i] = u[i] - alpha * f/fu;
  }
  for(i=space-mutation_jumpwidth;i<space-1;i++) {
    f = 0.5*wavespeed/dx*u[i-1] + (-mutation_rate + (i-space0)*dx)*u[i] - 0.5*wavespeed/dx*u[i+1] + mutation_rate*0.5*(i-space0+mutation_jumpwidth)*dx - 2.*u[i]*u[i];
    fu = (-mutation_rate + (i-space0)*dx) - 4.*u[i];
    unew[i] = u[i] - alpha * f/fu;
  }
  f = 0.5*wavespeed/dx*u[space-2] + (-mutation_rate + (space-1-space0)*dx)*u[space-1] - 0.5*wavespeed/dx*(space-space0)*0.5*dx + mutation_rate*0.5*(space-1-space0+mutation_jumpwidth)- 2.*u[space-1]*u[space-1];
  fu = (-mutation_rate + (space-1-space0)*dx) - 4.*u[space-1];
  unew[space-1] = u[space-1] - alpha * f/fu;

  memcpy(&u[0],&unew[0],space*sizeof(double));  
}






void set_c2_bc() {
  int i;
  for(i=0;i<space;i++) {
    c2[i][0] = 0.;
    c2[space-1][i] = 0.;
  }
}


void set_phi2_bc() {
//   int i,j;
//   double xmax,y;
//   xmax = 2.*u[space-3] - u[space-2];
//   for(i=0;i<space0;i++) {
//     phi2[i][0] = 0.;
//     phi2[space-1][i] = 0.;
//   }
//   for(i=space0;i<space;i++) {
//     phi2[i][0] = 0.;
//     phi2[space-1][i] = u[i]*xmax;
//   }
}
   
    


// ===========================================================================================
// initialize all variables
// ===========================================================================================

void initialize() {
  int i;
  double x;
  
  // ==================================================
  // initilization for C mode
  // ==================================================
  
  
  if (switchc == 1) {
    if(have_c_infile == 1) {
      read_c2_file(1);
      read_u_file(0);
    }else{
      read_u_file(1);
      c = (double*)malloc(space*sizeof(double));
      c2=(double**)malloc(space*sizeof(double*));
      for(i=0;i<space;i++)c2[i] = (double*)malloc((i+1)*sizeof(double));
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
    c2new = (double**)malloc(space*sizeof(double*));
    for(i=0;i<space;i++)c2new[i] = (double*)malloc((i+1)*sizeof(double));
    
    c_coeff_m = (double*)malloc(space*sizeof(double));
    c_coeff_0 = (double*)malloc(space*sizeof(double));
    c_coeff_p = (double*)malloc(space*sizeof(double));
    
    
    if(mutation_type == 2)c_coeff_j = (double*)malloc(space*sizeof(double));
    
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      
      c_coeff_m[i] = -0.5*wavespeed/dx;
      c_coeff_0[i] = x - 4.*u[i];
      c_coeff_p[i] = 0.5*wavespeed/dx;
      
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
    
    set_c2_bc();
    
  }

  // ==================================================
  // initilization for U mode
  // ==================================================
  
  if (switchu == 1) {
    if (have_u_infile == 1) {
      read_phi2_file(1);
      read_u_file(0);
      if(have_c2_infile==1)read_c2_file(0);
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
    
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      u_coeff_m[i] = 0.5*wavespeed/dx;
      u_coeff_0[i] = x - 4.*u[i];
      u_coeff_p[i] = -0.5*wavespeed/dx;
      
      switch(mutation_type) {
	case 1:	u_coeff_m[i] += mutation_rate/(dx*dx);
		u_coeff_0[i] += -2.*mutation_rate/(dx*dx);
		u_coeff_p[i] += mutation_rate/(dx*dx);
		iterate_u = iterate_u_diffusion;
		break;
	case 2:	u_coeff_0[i] += -mutation_rate;
		u_coeff_p[i] +=  mutation_rate;
		iterate_u = iterate_u_jumps;
		break;
	default:print_error("other mutation models not yet implemented");
      }
    }
    
    set_phi2_bc();
    
    
  }
}



void print_options() {
  printf("###################################################\n");
  printf("# solve profiles and weight-functions numerically #\n");
  if(switchc == 1) {
    printf("# calculating C^(2) profile                       #\n");
  }else if(switchu==1) {
    printf("# calculating PHI^(2) profile                     #\n");
  }
  printf("###################################################\n");  
  printf("#\n");
  printf("# mutationmodel     = %d\n",mutation_type);
  if(mutation_type > 1) {
    printf("# option            = %d\n",mutation_jumpwidth);
  }
  if(switchc == 1) {
    printf("#\n");
    printf("###################################################\n");
    if(verbosity>0) {
      printf("#     step\t xc              \t xg              \t popsize \t Tc\n");
    }
  }
}

void print_u() {
  int i;
  for(i=0;i<space;i++) {
    printf("%lf %e\n",(i-space0)*dx,u[i]);
  }
}


// ===========================================================================================
// adjust the density to fulfil the constraint \int dx dy u(x) u(y) c2(x,y) = 1
// ===========================================================================================

void adjust_constraint() {
  int i,j;
  double sum_offd = 0.,sum_d = 0.;
  double sum = 0.;
  double inv;
  double xc = 0.,xg = 0.,popsize = 0.;
  for(i=0;i<space;i++) {
    for(j=0;j<i;j++) {
      sum_offd += u[i]*u[j]*c2[i][j];	// calculate offdiagonal elements
    }
    sum_d += u[i]*u[i]*c2[i][i];	// diagonal elements
  }
  sum = (2.*sum_offd + sum_d)*dx*dx;	// offdiagonal elements appear twice (!) in the sum, diagonal ones not
  inv = 1./sum;
  
  for(i=0;i<space;i++) {
    for(j=0;j<=i;j++) {
      c2[i][j] *= inv;
    }
  }
}


void adjust_constraint_phi() {
  int i,j;
  double sum_offd = 0.,sum_d = 0.;
  double sum = 0.;
  double inv;
  
  for(i=0;i<space;i++) {
    for(j=0;j<i;j++) {
      sum_offd += phi2[i][j]*c2[i][j];
    }
    sum_d += phi2[i][i]*c2[i][i];
  }
  inv = 1./(dx*dx*(2.*sum_offd+sum_d));
  for(i=0;i<space;i++){
    for(j=0;j<=i;j++)phi2[i][j] *= inv;
  }
}
  

void get_g_peak(int step) {
  int i;
  double xg = 0.,xc = 0.,popsize = 0.;
  double tc = 0;
  for(i=0;i<space;i++) {
    xg += (i-space0)*u[i]*c[i];
    xc += (i-space0)*c[i];
    popsize += c[i];
    tc += u[i]*u[i]*c[i];
  }
  xc *= dx/popsize;
  xg *= dx*dx;
  tc *= dx;
  popsize *= dx;
  printf("%10d\t%17.10e\t%17.10e\t%17.10e\t%17.10e\n",step,xc,xg,popsize,tc);
}
  


// ===========================================================================================
// cleanup
// ===========================================================================================

void cleanup() {
  int i;
  if(switchc == 1) {
    free(c);
    free(cnew);
    for(i=0;i<space;i++) {
      free(c2[i]);
      free(c2new[i]);
    }
    free(c2);
    free(c2new);
    free(c_coeff_m);
    free(c_coeff_0);
    free(c_coeff_p);
    if(mutation_type == 2)free(c_coeff_j);
    free(u);
  }else{
    free(u);
    for(i=0;i<space;i++) {
      free(phi2[i]);
      free(phi2new[i]);
      if(have_c2_infile==1)free(c2[i]);
    }
    free(phi2);
    free(phi2new);
    if(have_c2_infile==1)free(c);
    if(have_c2_infile==1)free(c2);
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
  
  
  print_options();
  if(switchu == 1) {
    for(i=1;i<=maxsteps;i++) {
      iterate_u(i);
      if(have_c2_infile==1)if(i%adjustcontraintsteps == 0)adjust_constraint_phi();
    }
    if(have_phi_outfile == 1)write_phi2_file();
  }else if(switchc == 1) {
    adjust_constraint();
    for(i=1;i<=maxsteps;i++) {
      iterate_c(i);
      if(i%adjustcontraintsteps == 0) {
	adjust_constraint();
	if(verbosity>0)get_g_peak(i);
      }
    }
    if(have_c_outfile == 1)write_c2_file();
  }
  cleanup();
  return 0;
}



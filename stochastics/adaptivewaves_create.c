// adaptivewaves_create.c
// v0.01	120827	modify for adaptive waves
//			include versioning of configuration files
//			can provide -u FILENAME, so that density c[i] is already scaled correctly

// v0.02	120906	modified the binary file structure:
// 			first two (int) entries are number of integer and double parameters
//			then dx, space, space0
//			all integer parameters
//			all double parameters
//			density[space]


// v0.03,	120907	format = (0,4) (rho,starttime,diffusionconstant,driftvelocity)
//			new options -D, -V

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int p_int = 0,p_double=4;
double rho = 0.0;
double dx = -1;
int space = 10000;
int space0;
double starttime = 0.;

double diffusionconstant = 1.0;
double convectionvelocity = 0.0;

double zeropos = 0.5;

double sigma2 = -1.;

int mode = 3;
char filename[128];

int quiet = 0;


int usefile = 0;

char ufilename[128];
int read_u_from_file = 0;

double *c,*u;


double populationsize=1.;
double populationvariance = 1e-4;

int setvariance = 0;

int print_error(char *msg) {
  fprintf(stderr,"ERROR: %s\n",msg);
  exit(1);
}


void parsecommandline(int argn, char *argv[]) {
  char c;
  while((c = getopt(argn, argv,"r:d:s:z:m:o:qN:u:S:v:D:N:V:")) != -1){
    switch(c) {
      case 'r': rho = atof(optarg);
		break;
      case 'd': dx = atof(optarg);
		break;
      case 's': space = atoi(optarg);
		break;
      case 'z': space0 = atoi(optarg);
		break;
      case 'm': if(strcmp(optarg,"flat")==0) {
		  mode = 2;
		}else if(strcmp(optarg,"gauss")==0) {
		  mode = 3;
		}else if(strcmp(optarg,"delta")==0) {
		  mode = 4;
		}else print_error("mode not known!");
		break;
      case 'o':	strcpy(filename,optarg);
		usefile = 1;
		break;
      case 'q':	quiet = 1;
		break;
      case 'u': strcpy(ufilename,optarg);
		read_u_from_file = 1;
		break;
      case 'S': sigma2 = atof(optarg);
		break;
      case 'D': diffusionconstant = atof(optarg);
		break;
      case 'v': convectionvelocity = atof(optarg);
		break;
      case 'N':	populationsize = atof(optarg);
		break;
      case 'V':	populationvariance = atof(optarg);
		setvariance = 1;
		break;
    }
  }
  if(dx < 0.)print_error("dx not set! use option -d");
  space0 = zeropos * space;
  if(usefile == 0)strcpy(filename,"-");
}

int get_u_from_file() {
  FILE *fpu;
  double ux;
  int i,tmp_count_params[2];
  int *tmpi;
  double *tmpd;
  int u_space,u_space0;
  double u_dx;
  
  fpu = fopen(ufilename,"rb");
  if(fpu != NULL) {
    fread(&tmp_count_params[0],2,sizeof(int),fpu);
    
    fread(&u_dx,1,sizeof(double),fpu);
    fread(&u_space,1,sizeof(int),fpu);
    fread(&u_space0,1,sizeof(int),fpu);
    
    if(tmp_count_params[0] > 0) {
      tmpi = (int*)malloc(tmp_count_params[0]*sizeof(int));
      fread(&tmpi[0],tmp_count_params[0],sizeof(int),fpu);
      free(tmpi);
    }
    if(tmp_count_params[1] > 0) {
      tmpd = (double*)malloc(tmp_count_params[1]*sizeof(double));
      fread(&tmpd[0],tmp_count_params[1],sizeof(double),fpu);
      free(tmpd);
    }
    
    if((space!= u_space) || (space0 != u_space0)) {
      print_error("reading u: lattice does not match");
    }
    u=(double*)malloc(space*sizeof(double));
    fread(&u[0],space,sizeof(double),fpu);
    
    fclose(fpu);
  }
  return 0;
}



double get_gaussian(double x, double sigma2) {
  int i;
  return 0.398942/sqrt(sigma2)*exp(-x*x/(2.*sigma2));
}

void initialize() {
  int i,j;
  double u2int_inv,uint_inv;
  double x;
  
  if(sigma2<0.)sigma2 = 0.000025*space*space*dx*dx;
  if(setvariance == 1)sigma2 = populationvariance;
  printf("sigma2 = %lf\n",sigma2);
  
  c=(double*)malloc(space*sizeof(double));
  if(mode == 3) {
    for(i=0;i<space;i++) {
      x = (i-space0)*dx;
      c[i] = populationsize*get_gaussian(x,sigma2);
    }
  }else if(mode == 2) {
    for(i=0;i<space;i++) {
      c[i] = 1.;
    }
  }else if(mode == 4) {
    for(i=0;i<space;i++) {
      c[i]=0.;
    }
    c[space0] = 1.;
  }else{
    print_error("could not initialize density c(x)");
  }
  
}


void rescale_c_with_u() {
  int i;
  double int_uc = 0.;
  
  for(i=0;i<space;i++) {
    int_uc += u[i] * c[i];
  }
  int_uc *= dx;
  
  for(i=0;i<space;i++) {
    c[i] /= int_uc;
  }
}



double get_N() {
  int i;
  double sum = 0.;
  for(i=0;i<space;i++) {
    sum += c[i];
  }
  return sum*dx;
}

void write_data() {
  FILE *fp;
  if(usefile == 0) {
    fp = stdout;
  }else{
    fp = fopen(filename,"wb");
  }
  fwrite(&p_int,1,sizeof(int),fp);
  fwrite(&p_double,1,sizeof(int),fp);

  fwrite(&dx,1,sizeof(double),fp);
  fwrite(&space,1,sizeof(int),fp);
  fwrite(&space0,1,sizeof(int),fp);
  
  fwrite(&rho,1,sizeof(double),fp);
  fwrite(&starttime,1,sizeof(double),fp);
  fwrite(&diffusionconstant,1,sizeof(double),fp);
  fwrite(&convectionvelocity,1,sizeof(double),fp);
  
  fwrite(&c[0],space,sizeof(double),fp);
  if(usefile == 1)fclose(fp);
}

void print_options() {
  double tmp_n;
  fprintf(stderr,"creating binary data file for stochastic simulation of an oasis...\n");
  fprintf(stderr,"options:\n");
  fprintf(stderr,"  r           = %6.4e\n",rho);
  fprintf(stderr,"  dx          = %lf\n",dx);
  fprintf(stderr,"  space       = %d\n",space);
  fprintf(stderr,"  space0      = %d\n",space0);
  fprintf(stderr,"other parameters:\n");
  if(read_u_from_file == 1) {
    fprintf(stderr,"  u*-FILENAME = %s\n",ufilename);
  }
  tmp_n = get_N();
  fprintf(stderr,"  N           = %lf\n",tmp_n);
}


void cleanup() {
  free(c);
  if(read_u_from_file)free(u);
}

int main(int argn, char* argv[]) {
  parsecommandline(argn,argv);
  initialize();
  if(read_u_from_file == 1) {
    get_u_from_file();
    rescale_c_with_u();
  }
  write_data();
  if(quiet!=1)print_options();
  cleanup();
  return 0;
}

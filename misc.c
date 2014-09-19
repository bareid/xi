#include "header.h"

extern gsl_rng *rngen;
extern int rngensetup;


real my_pow_2(real x)  {
    return x*x;
    }

real my_pow_3(real x)  {
    return x*x*x;
    }

real my_pow_4(real x)  {
    return my_pow_2(x*x);
    }

real my_pow_6(real x)  {
    return my_pow_2(x*x*x);
    }

real fixperiodicsep(real sepval, real Lbox)  {
  if(fabs(sepval + Lbox) < fabs(sepval))  {
    sepval = sepval + Lbox;
    }
  if(fabs(sepval - Lbox) < fabs(sepval))  {
    sepval = sepval - Lbox;
    }
  return sepval;
  }

real fixbox(real coord,real Lbox)  {
  while((coord >= Lbox) || (coord < 0.))  {
    if(coord >= Lbox)  {
      coord = coord - Lbox;
      }
    if(coord < 0.)  {
      coord = coord + Lbox;
      }
    }
  assert((coord >= 0.) && (coord < Lbox));
  return coord;
  }

FILE *open_file_write(char *filename) {
  FILE *ofp;

  if(!(ofp=fopen(filename,"w")))  {
    fprintf(stderr,"Can't open %s\n",filename);
    return NULL;
    }
  else {
    return ofp;
    }
  }

FILE *open_file_append(char *filename) {
  FILE *ofp;

  if(!(ofp=fopen(filename,"a")))  {
    fprintf(stderr,"Can't open %s\n",filename);
    return NULL;
    }
  else {
    return ofp;
    }
  }
  
FILE *open_file_read(char *filename) {
  FILE *ifp;

  if(!(ifp=fopen(filename,"r")))  {
    fprintf(stderr,"Can't open %s\n",filename);
    return NULL;
    }
  else {
    return ifp;
    }
  }
  
int get_file_length(FILE *ifp, int *headercount, int *linecount) {
  char line[MAXLINELEN];
  
  *headercount = 0;
  *linecount = 0;
      
    while(fgets(line,MAXLINELEN,ifp))  {
        if(!feof(ifp))  {
            (*linecount)++;
            }
    else {
      break;
      }
    if(strstr(line,"#"))  {
      (*headercount)++;
      }
        }
  (void) rewind(ifp);
  return 0;
  }


real legendrep(int ell, real x) {

  real xsqr;
  switch(ell) {
    case(0):
  return (1.);
  break;
    case(1):
  return (x);
  break;
    case(2):
  return (-0.5+1.5*x*x);
  break;
    case(3):
  return (-1.5*x+2.5*x*x*x);
  break;
    case(4):
  xsqr = x*x;
  return (0.375 - 3.75*xsqr + 4.375*xsqr*xsqr);
  break;
    case(6):
  xsqr = x*x;
  return (-0.3125+6.5625*xsqr - 19.6875*xsqr*xsqr + 14.4375*xsqr*xsqr*xsqr);
  break;
    default:
  exit(1);
  break;
    }
  }

long double longlegendrep(int ell, long double x) {

  long double xsqr;
  switch(ell) {
    case(0):
  return (1.);
  break;
    case(1):
  return (x);
  break;
    case(2):
  return (-0.5+1.5*x*x);
  break;
    case(3):
  return (-1.5*x+2.5*x*x*x);
  break;
    case(4):
  xsqr = x*x;
  return (0.375 - 3.75*xsqr + 4.375*xsqr*xsqr);
  break;
    case(6):
  xsqr = x*x;
  return (-0.3125+6.5625*xsqr - 19.6875*xsqr*xsqr + 14.4375*xsqr*xsqr*xsqr);
  break;
    default:
  exit(1);
  break;
    }
  }



//Random number generator stuff.

unsigned int devrand(void)  {
  int fn;
  unsigned int r2=1;
  fn = open("/dev/urandom",O_RDONLY);
  if(fn == -1)  {
    fprintf(stderr,"couldn't open /dev/urandom file.\n");
    exit(1);
    }
  //unsigned long int
  if(read(fn,&r2,sizeof(unsigned int)) != sizeof(unsigned int))  {
    fprintf(stderr,"bad /dev/urandom read of %lu bytes\n",sizeof(unsigned int));
    exit(1);
    }
  //printf("%lu\n",r2);
  //crashes when i have this; don't know why.
  close(fn);
  return r2;
  }

// This initializes the random number type rngen
int initialize_rngen(unsigned int unsignedintseed)  {
  if(rngensetup == 1) {
    return 0;
    }

  rngen = gsl_rng_alloc(gsl_rng_mt19937);
/*
  if(fixrngenseed == 1) {
    unsignedintseed = 408606142;
    }
*/

  gsl_rng_set(rngen,unsignedintseed);

  #ifdef REALLYVERBOSE
  printf("Initialized random number generator mt19937 with %u\n",unsignedintseed);
  #endif

  double testval = gsl_rng_uniform(rngen);
  assert(testval >= 0.);
  assert(testval <= 1.);
  testval = gsl_rng_uniform(rngen);
  assert(testval >= 0.);
  assert(testval <= 1.);

  //end setup GSL RNG instead
  rngensetup = 1;
  return 0;
  }

int free_rngen() {
  if(rngensetup == 1)  {
    gsl_rng_free(rngen);
    }
  rngensetup = 0;
  return 0;
  }

//end Random number generator stuff.

//copying these from cosmoxi2d.
double **malloc2ddouble(int l1, int l2) {
  double **m;
  int i;
  m = malloc(l1*sizeof(double *));
  m[0] = malloc(l1*l2*sizeof(double));
  for(i=0;i<l1; i++) {
    m[i] = m[0]+i*l2;
    }
  return m;
  }

double malloc2dfree(double **m) {
  if(m) {
    if(m[0]) {
      free(m[0]);
      }
    free(m);
    }
  }

int **malloc2dint(int l1, int l2) {
  int **m;
  int i;
  m = malloc(l1*sizeof(int *));
  m[0] = malloc(l1*l2*sizeof(int));
  for(i=0;i<l1; i++) {
    m[i] = m[0]+i*l2;
    }
  return m;
  }

int malloc2dfreeint(int **m) {
  if(m) {
    if(m[0]) {
      free(m[0]);
      }
    free(m);
    }
  }


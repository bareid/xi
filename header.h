//basic stuff needed by everyone.
//supported command line options: REALFLOAT, VOPT.

#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_errno.h>

#ifndef NOOMP
#include "omp.h"
#endif

//for drawing random numbers to downsample the random catalog.
#include <fcntl.h>            // this is to get the O_RDONLY to work in the random seed.
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//testing!!
#define SANITYCHECKS
#define ASSERT_ON
//end testing.

#ifdef REALFLOAT
#define real float
#else
#define real double
#endif

#define MAXLINELEN 1000
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

//remove later!
#define REALLYVERBOSE
#define SANITYCHECKS
#define ASSERT_ON
//end remove later!

typedef struct {
  double omegam;
  double omegal;
  double omegak;
  double h;
  } cosmo_params;

typedef struct {
  real pos[3];
#ifdef VOPT
  real vel[3];
#endif
  float wgt;
  int cell;
  int DorR;  //this can be data/random or galaxy/mass.
  } cntparticle;  //let's absorb cattype into this datatype as well.

typedef struct {
//Yu will hardcode these to be double?  think about this later.
    real pos[3];    
    real vel[3]; //leave blank for ra/dec/z counts.
    real weight;
} particle;

//from bethalexie code.
typedef struct {
  int bintype;
  int logxopt;
  int nx;
  double minx;
  double dx;
  int logyopt;
  int ny;
  double miny;
  double dy;
  int realorzspace;
  int periodicopt;
  int nbins2d;
  //wp only:
  double rpimax;
  int zspaceaxis;
  //if your data points are not simple log or linear-spaced, these parameters specify how to compute the theory in the data bins by combining points measured in the simple binning.
  //delete everything below here, do in python?
  int rebinopt;
  int maskopt;
  int *binxiell;
  int *mask;
  double rperpcut;  //masks points below this perpendicular separation; not used except for a sanity check, since the mask is computed outside of this program.
  int nr0;
  int nr2;
  int nr4;
  int ellmaxdata; //this is 0 for lensing, currently 2 (could be 4 if we wanted) for RSD
  int ndata;
  } xibindat;

typedef struct {
  int ntheta;
  real *thetacen, *thetamin, *thetamax, *angweight;
  real dlog10theta;
  real log10thetamin;
  } angwgt;

typedef struct {
  real RSEPMAX,RSEPMAXAP,rminsqr,rmaxsqr;  // min and max scales included in pair counter.
  real APscale[3];
  int Ncell;
  int NFAC;
  real Lbox;
  real originpos[3];
  } cntparams;

//global variables, defined in header.c
//no more global variables!
/*
extern real Lbox;
extern real originpos[3];
extern int Ncell;
extern int NFAC;
*/

//functions from misc.c
real my_pow_2(real x);
real my_pow_3(real x);
real my_pow_4(real x);
real my_pow_6(real x);
real fixperiodicsep(real sepval, real Lbox);
real fixbox(real coord, real Lbox);
FILE *open_file_write(char *filename);
FILE *open_file_append(char *filename);
FILE *open_file_read(char *filename);
int get_file_length(FILE *ifp, int *headercount, int *linecount);
real legendrep(int ell, real x);
long double longlegendrep(int ell, long double x);
//rngen stuff.
unsigned int devrand(void);
int initialize_rngen(unsigned int unsignedintseed);
int free_rngen();
double **malloc2ddouble(int l1, int l2);
double malloc2dfree(double **m);
int **malloc2dint(int l1, int l2);
int malloc2dfreeint(int **m);
//end functions from misc.c

//from cosmoutils.c
int docosmosetup();
int freecosmo();
double comoving_disthinvMpc_lcdm(double z,cosmo_params cosmopfid);
double comoving_distMpc_lcdm(double z, cosmo_params cosmopfid);
//end from cosmoutils.c

//from xibin.c
xibindat readbinfile(char *frebin);
int testmaskcutsmallscale(xibindat b);
//end from xibin.c

//functions from readnbody.c
particle *readsimcat(char *fname, int whichcat, real Lbox, int zspaceaxis, real vscale, float Mmin, float Mmax, int DorR, int *ngals);
//end functions from readnbody.c

//functions from readradecz.c
void radecztopos(double ra, double dec, float z, int unitsMpc, int angopt, cosmo_params cosmopfid, real pos[3], real *chi);
particle *readcat(char *infilename, int ftype, int unitsMpc, int angopt, cosmo_params cosmopfid, float zmin, float zmax, int DorR, float ndownRR, int *ntot, long double *wgttot, real *maxdist);
//end functions from readradecz.c

//functions from computesepnbody.c
int addpairperiodic(const real pos1[], const real pos2[], const real vel1[], const real vel2[], xibindat b, cntparams *cp, int *bin2d, long double *vr, long double *v2perp, long double *v2par);
//end functions from computesepnbody.c

//functions from computesepradecz.c
real getrsep3dmaxforang(xibindat b, const real originpos[]);
int addpairsky2d(const real pos1[], const real pos2[], const real originpos[], xibindat b, int *bin2d);
int addpairsky3d(const real pos1[], const real pos2[], xibindat b, cntparams *cp, int *bin2d, real *angsep);
//end functions from computesepradecz.c

//functions from angupweight.c
angwgt readangwgt(char *filename1pw);
real fbangweight(real asep, angwgt awgt);
void freeangwgt(angwgt awgt);
//end functions from angupweight.c

//functions from xigrid.c
int i2n(const int ipos[],int Ncell);
void n2i(int cellnum, int Ncell, int *ipos);
int chksortindxmatch(int icell, int jcell, int Ncell);
int wrapgridperiodic(int ii, int Ncell);
void nbrsprecompute(int *closenbrs, size_t nnbrsmax, int NFAC);
void getnbrs(int cellnum, int NFAC, int Ncell, int *nbrs, int *nnbrs, int *closenbrs, int periodicopt);
//end functions from xigrid.c

//functions from countpairs.c
//need to import these first!
int countpairsradecz(particle *p1, int np1, particle *p2, int np2, real maxdist, xibindat b, angwgt awgt, long double *Npairsfinal);
int countpairssim(particle *p1, int np1, particle *p2, int np2, real Lbox, xibindat b, real APperp, real APpar, long double *Npairsfinal);
//don't expose countpairs.c to the outside world.

//functions from output.c
//BEWARE, THESE ARE UNTESTED!!
int xiellperiodic(long double *Npairsfinal, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi);
int vstatsperiodic(long double *Npairsfinal, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **vstats);
int xiellperiodicrebin(long double *Npairsfinal,xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi);
int xiellNprebin(long double *Npairsfinal, xibindat b, long double **Npout);
void printxiell(char *outfname, xibindat b, double **xi);
void printxiellrebin(char *outfname, xibindat b, double *r0D, double *r2D, double *r4D, double *xi);
void printxiellwvol(char *outfname, xibindat b, double **xi);
int xiellperiodicrebincutsmallscale(long double *Npairsfinal,xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi);
int xiellNprebincutsmallscale(long double *Npairsfinal, xibindat b, long double **Npout);
int xigridperiodic(long double *Npairsfinal, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi);
void printxigrid(char *outfname, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi);
void printxigridwNp(char *outfname, char *outfnameNp, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi, long double *Npairsfinal);
int wpperiodic(long double *Npairsfinal, xibindat b, int n_gal, real Lbox, real APperp, real APpar, double *wp);
int wpNprebin(long double *Npairsfinal, xibindat b, long double **Npout);
void printwp(char *outfname, xibindat b, double *wp);
void printNpairsgeneric(char *foutbase, long double *Npairsfinal,xibindat b,int DRopt,int n1, long double n1wgt, long double n2wgt,char *binfname, double omfid, double hfid);
void printNpairssim(int n1, int n2, int autoorcross, xibindat b, real Lbox, real APperp, real APpar, char *foutbase, long double *Npairsfinal);
//end functions from output.c






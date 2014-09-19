#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>


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
  } xibindat;

typedef struct {
  double *px;  //these are pointers to 1d arrays for simplicity; sending a 2d array seems like a nightmare.
  double *py;
  double *pz;
  double *vx;
  double *vy;
  double *vz;
  double *wgt;
  int np;
} dataset;





#include "test.h"

int printit(xibindat b) {
  printf("bintype: %d\n",b.bintype);

  printf("logxopt: %d\n",b.logxopt);
  printf("nx: %d\n",b.nx);
  printf("minx: %e\n",b.minx);
  printf("dx: %e\n",b.dx);

  printf("logyopt: %d\n",b.logyopt);
  printf("ny: %d\n",b.ny);
  printf("miny: %e\n",b.miny);
  printf("dy: %e\n",b.dy);

  printf("realorzspace: %d\n",b.realorzspace);
  printf("periodicopt: %d\n",b.periodicopt);
  printf("nbins2d: %d\n",b.nbins2d);
  return 0;
  }

int hellop(dataset d1) {

  printf("hey beth: %d\n",d1.np);

  int i,j;
  int mymax = d1.np;
  if(d1.np > 10) {
    mymax = 10;
    }
  for(i=0;i<mymax;i++) {
    printf("%e ",d1.px[i]);
    printf("%e ",d1.py[i]);
    printf("%e ",d1.pz[i]);
    printf("%e ",d1.vx[i]);
    printf("%e ",d1.vy[i]);
    printf("%e ",d1.vz[i]);
    printf("%e\n",d1.wgt[i]);
    }
  printf("hey beth: %d\n",d1.np);
  return 0;
  }

int testsingle(double *px, int np) {
  int i,j;
  int mymax = np;
  if(np > 10) {
    mymax = 10;
    }
  for(i=0;i<mymax;i++) {
    printf("%e\n",px[i]);
    }
  printf("hey beth: %d\n",np);
  return 0;
  }

#include "header.h"

real mydotprod(real *r1, real *r2) {
  return (r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]);
  }

//velocity stuff only if bintype == 2.
int addpairperiodic(const real pos1[], const real pos2[], const real vel1[], const real vel2[], xibindat b, cntparams *cp, int *bin2d, long double *vr, long double *v2perp, long double *v2par) {
  real rsep;
  real rvec[3];
  int xbin, ybin;
  real mymu2;

//variables needed for velocities.
  real rnorm[3];
  real deltav[3];
  real dvsqr;
  int vi; 

  int i1,i2,i3;
  double rpi, rsig;

  rvec[0] = fixperiodicsep((pos2[0]-pos1[0]),cp->Lbox)/cp->APscale[0];
  //printf("ib 0 %e %e %e\n",pos1[0], pos2[0], rvec[0]);
  //fflush(stdout);
  if(fabs(rvec[0]) >= cp->RSEPMAX) {
    (*bin2d) = -1;
    return 0;
    }
  rvec[1] = fixperiodicsep((pos2[1]-pos1[1]),cp->Lbox)/cp->APscale[1];
  //printf("ib 1 %e %e %e\n",pos1[1], pos2[1], rvec[1]);
 // fflush(stdout);
  if(fabs(rvec[1]) >= cp->RSEPMAX) {
    (*bin2d) = -1;
    return 0;
    }
  rvec[2] = fixperiodicsep((pos2[2]-pos1[2]),cp->Lbox)/cp->APscale[2];
//  printf("ib 2 %e %e %e\n",pos1[2], pos2[2], rvec[2]);
//  fflush(stdout);
  if(fabs(rvec[2]) >= cp->RSEPMAX) {
    (*bin2d) = -1;
    return 0;
    }

  switch(b.bintype) {
    case(0): //xi_ell
  rsep = rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2];
  if(rsep >=cp->rmaxsqr) {
    (*bin2d) = -1;
    return 0;
    }
  if(rsep <cp->rminsqr) {
    (*bin2d) = -1;
    return 0;
    }
  //rsep is still rsep**2.  compute mymu2.
  mymu2 = rvec[b.zspaceaxis]*rvec[b.zspaceaxis]/rsep;
  //fix edge cases.
  (mymu2) = min((mymu2),1-1.0e-9);
  (mymu2) = max((mymu2),1.0e-9);
  ybin = min(((int) floor((sqrt(mymu2)*b.ny))),b.ny-1);

  rsep = sqrt(rsep);
  if(b.logxopt == 0) {
    xbin = min((int) floor((rsep-b.minx)/b.dx),b.nx-1);
    }
  else {
    xbin = min((int) floor((log10(rsep)-b.minx)/b.dx),b.nx-1);
    }
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(xbin < b.nx && xbin >= 0);
  assert(mymu2 >= 0. && mymu2 <= 1.);
  assert(ybin < b.ny && ybin >= 0);
  #endif
  if(xbin == 1) {
    assert(b.zspaceaxis == 2);
    printf("%e %e %e %e %e %e %e %e %e %e %e %d %d\n",pos1[0],pos2[0],rvec[0],pos1[1],pos2[1],rvec[1],pos1[2],pos2[2],rvec[2],rsep,mymu2,xbin,ybin);
    }

  break;    
    case(1): //xi_grid.

  switch(b.zspaceaxis) {
    case(0):
  i1 = 0; i2 = 1; i3 = 2;
  break;
    case(1):
  i1 = 1; i2 = 2; i3 = 0;
  break;
    case(2):
  i1 = 2; i2 = 0; i3 = 1;
  break;
    default:
  exit(1);
  } //end switch zspaceaxis.

  rpi = fabs(rvec[i1]);
  if(b.logyopt == 1 && rpi > 0.) {
    rpi = log10(rpi);
    }
  if(rpi >= b.miny + b.ny*b.dy || rpi < b.miny) {
    (*bin2d) = -1;
    return 0;
    }
  rsig = sqrt(rvec[i2]*rvec[i2]+rvec[i3]*rvec[i3]);

//deal with pairs at same position, rsig <= 0.
/*
  if(b.logxopt == 1 && rsig > 0.) {
    rsig = log10(rsig);
    }
*/

  if(b.logxopt == 1) {
    if(rsig <= 0.) {
      (*bin2d) = -1;
      return 0;
      }
    rsig = log10(rsig);
    }

  if(rsig >= b.minx + b.nx*b.dx || rsig < b.minx) {
    (*bin2d) = -1;
    return 0;
    }

  ybin = min((int) floor((rpi-b.miny)/b.dy),b.ny-1);
  xbin = min((int) floor((rsig-b.minx)/b.dx),b.nx-1);
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(ybin >= 0 && ybin < b.ny);
  assert(xbin >= 0 && xbin < b.nx);
  #endif
  break; //end xigrid.
    case(2):  //velocity option.
  rsep = rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2];
  if(rsep >=cp->rmaxsqr) {
    (*bin2d) = -1;
    return 0;
    }
  if(rsep <cp->rminsqr) {
    (*bin2d) = -1;
    return 0;
    }
  rsep = sqrt(rsep);
  if(b.logxopt == 0) {
    xbin = min((int) floor((rsep-b.minx)/b.dx),b.nx-1);
    }
  else {
    xbin = min((int) floor((log10(rsep)-b.minx)/b.dx),b.nx-1);
    }
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(xbin < b.nx && xbin >= 0);
  #endif
  ybin = 0;

  //now compute velocity statistics.
  dvsqr = 0.;
  for(vi=0;vi<=2;vi++) {
    rnorm[vi] = rvec[vi]/rsep;
    deltav[vi] = (vel2[vi]-vel1[vi]);
    dvsqr += deltav[vi]*deltav[vi];
    }
  *vr = ((long double) mydotprod(deltav,rnorm));
  *v2par = (*vr)*(*vr);
  *v2perp = ((long double) ((dvsqr) - (*v2par)));   
  break;
    default:
  exit(1);
  } //end switch over bintypes.
  (*bin2d) = xbin*b.ny + ybin;
  return 1;
  } //end addpairperiodic.

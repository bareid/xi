#include "header.h"

//these are all defined/setup in countpairs.c
/*
extern real RSEPMAX,rminsqr,rmaxsqr;
extern real originpos[3];
*/


real angsep(real ra1, real dec1, real ra2, real dec2) {
  real tmp1 = sin(dec1)*sin(dec2)+cos(ra1-ra2)*cos(dec1)*cos(dec2);
  if(tmp1>1.0) tmp1=1.0;
  if(tmp1<-1.0) tmp1=-1.0;
  #ifdef SANITYCHECKS
  assert(tmp1 <= 1.0 && tmp1 >= -1.);
  #endif
  return(acos(tmp1)); //in radians
  }

void radectounitpos(double ra,double dec, const real originpos[], real pos[3]) {
  //this is also defined in readradecz.c/radecztopos; separate version here so don't need to input cosmo_params.
  (pos)[0] = cos(dec)*cos(ra)+originpos[0];
  (pos)[1] = cos(dec)*sin(ra)+originpos[1];
  (pos)[2] = sin(dec)+originpos[2];
  }

//BR note: take out long double (replace with "real" = float/double) for angular separation calculation?
real angsepfromunitvec(const real pos1[], const real pos2[], const real originpos[]) {
  long double ndot = 0.;
  int ii;
  for(ii=0;ii<=2;ii++) {
    ndot += (pos1[ii]-originpos[ii])*(pos2[ii]-originpos[ii]);
    }
  if(ndot > 1.0) {
    ndot = 1.0;
    }
  if(ndot < -1.0) {
    ndot = -1.0;
    }
  return ((real) acos(((real) ndot)));
  }

real getrsep3dmaxforang(xibindat b, const real originpos[], real mindist) {
  real angsepmax;
  if(b.logxopt == 1) { //log binning.
    angsepmax = pow(10.,b.minx + ((real) b.nx)*b.dx);
    }
  else { //linear binning.
    angsepmax = (b.minx + ((real) b.nx)*b.dx);
    }
  if(b.bintype == 4) {
    //angsepmax is currently in rp units.
    angsepmax = 2.*asin(0.5*angsepmax/mindist);
    }

  //make sure angsepmax is in radians!
  assert(angsepmax < M_PI);
  assert(angsepmax > 0.);

  real pos1t[3];
  real pos2t[3];
  real pos3t[3];

  real diff, diffbig;
  real ratmp2, dectmp2;
  real ratmp1, dectmp1;
  real ratmp3, dectmp3;
  long double rsqr, rsqrbig, n1, n2, ndot;
  int ii;

  ratmp1 = 0.2;
  dectmp1 = 0.01;
  radectounitpos(ratmp1, dectmp1, originpos, pos1t);
  dectmp2 = dectmp1 + angsepmax;
  dectmp3 = dectmp1 + angsepmax*1.1;
  ratmp2 = ratmp1;
  ratmp3 = ratmp1;
  radectounitpos(ratmp2, dectmp2, originpos, pos2t);
  radectounitpos(ratmp3, dectmp3, originpos, pos3t);
  rsqr = 0.;
  rsqrbig = 0.;
  n1 = 0.;
  n2 = 0.;
  ndot = 0.;
  for(ii=0;ii<=2;ii++) {
    diff = (pos2t[ii] - pos1t[ii]);
    diffbig = (pos3t[ii] - pos1t[ii]);
    rsqr += diff*diff;
    rsqrbig += diffbig*diffbig;
    n1 += (pos1t[ii]-originpos[ii])*(pos1t[ii]-originpos[ii]);
    n2 += (pos2t[ii]-originpos[ii])*(pos2t[ii]-originpos[ii]);
    ndot += (pos1t[ii]-originpos[ii])*(pos2t[ii]-originpos[ii]);
    }
  real a1 = angsep(ratmp1,dectmp1,ratmp2,dectmp2);
  real a2 = ((real) (ndot/sqrt((real) (n1*n2))));
  if(a2>1.0) a2=1.0;
  if(a2<-1.0) a2=-1.0;
  a2 = acos(a2); 
  assert(fabs(a1-a2)/a1 < 2.0e-6);  //fractional error on angular separation should be small!
  assert(fabs(a1 - angsepmax) < 2.0e-6);
  real a3d = sqrt((real) rsqr);
//  bethdidimean rsqrbig to return?  I guess not, this is exact!
  return a3d;
  }

int addpairsky2d(const real pos1[], const real pos2[], const real originpos[], xibindat b, int *bin2d) {

  real asep = angsepfromunitvec(pos1,pos2,originpos);
  int xbin;

  if(b.logxopt == 1) {
    if(asep <= 0.) {
      (*bin2d) = -1;
      return 0;
      }
    asep = log10(asep);
    }

  if(asep >= b.minx + b.nx*b.dx || asep < b.minx) {
    (*bin2d) = -1;
    return 0;
    }

  xbin = min((int) floor((asep-b.minx)/b.dx),b.nx-1); 
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(xbin >= 0 && xbin < b.nx);
  #endif
  (*bin2d) = xbin;
  return 1;
  }  //end addpairsky2d.

int addpairsky2dhogg(const real pos1[], const real pos2[], const real originpos[], xibindat b, int *bin2d, real chi) {

  real asep = angsepfromunitvec(pos1,pos2,originpos);
  real rsig = 2.*chi*sin(0.5*asep);

  //tmp!
/*
  printf("p1 %e %e %e\n",pos1[0],pos1[1],pos1[2]);
  printf("p2 %e %e %e\n",pos2[0],pos2[1],pos2[2]);
  printf("op,chi %e %e %e %e %e %e %e\n",originpos[0],originpos[1],originpos[2],chi,asep,rsig,asep*chi);
  fflush(stdout);
  assert(chi >= 570. && chi <= 1100.);
*/

  int xbin;
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

  xbin = min((int) floor((rsig-b.minx)/b.dx),b.nx-1);
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(b.bintype == 4);
  assert(xbin >= 0 && xbin < b.nx);
  #endif
  (*bin2d) = xbin;
  return 1;
  } //end addpairsky2dhogg.

//BR, make addpairsky faster!!
int addpairsky3d(const real pos1[], const real pos2[], xibindat b, cntparams *cp, int *bin2d, real *angsep) {

  (*angsep) = -1000.;

  real rsep;
  real losvec[3];
  real rvec[3];
  real rpi, rsig;

  int xbin, ybin;

//rest of this function is for xi_ell/xi_grid (3d correlation fxn). 

  rvec[0] = (pos2[0]-pos1[0]);
  if(fabs(rvec[0]) >= cp->RSEPMAX) {
    (*bin2d) = -1;
    return 0;
    }
  rvec[1] = (pos2[1]-pos1[1]);
  if(fabs(rvec[1]) >= cp->RSEPMAX) {
    (*bin2d) = -1;
    return 0;
    }
  rvec[2] = (pos2[2]-pos1[2]);
  if(fabs(rvec[2]) >= cp->RSEPMAX) {
    (*bin2d) = -1;
    return 0;
    }
  rsep = rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2];
  if(rsep >=cp->rmaxsqr) {
    (*bin2d) = -1;
    return 0;
    }
  if(rsep <=cp->rminsqr) {
    (*bin2d) = -1;
    return 0;
    }
  rsep = sqrt(rsep);

  //probably not necessary to be long double!!

  long double losnorm = 0.;
  long double rnorm = 0.;
  long double mymu2 = 0.;
  long double angsep1 = 0.;
  int i;
  long double r1norm, r2norm;
  r1norm = 0.;
  r2norm = 0.;
  for(i=0;i<=2;i++) {
    r1norm += (pos1[i]-cp->originpos[i])*(pos1[i]-cp->originpos[i]);
    r2norm += (pos2[i]-cp->originpos[i])*(pos2[i]-cp->originpos[i]);
    angsep1 += (pos1[i]-cp->originpos[i])*(pos2[i]-cp->originpos[i]);
    }
  r1norm = sqrt(r1norm);
  r2norm = sqrt(r2norm);
  angsep1 = angsep1/(r1norm*r2norm);  //BR Aug 15.  This is the cosine of the angular separation.  IS that what I want?  Yes, acos taken at end of this script.

  for(i=0;i<=2;i++) {
    losvec[i] = (pos1[i]-cp->originpos[i]) + r1norm/r2norm*(pos2[i]-cp->originpos[i]);
    losnorm += losvec[i]*losvec[i];
    rnorm += rvec[i]*rvec[i];
    mymu2 += rvec[i]*losvec[i];
    }

  mymu2 = mymu2*mymu2/(rnorm*losnorm);
  //fix boundary cases.
  mymu2 = min(mymu2,1-1.0e-9);
  mymu2 = max(mymu2,1.0e-9);

  if(angsep1>1.0) angsep1=1.0;
  if(angsep1<-1.0) angsep1=-1.0;
  (*angsep) = ((real) (acos(angsep1)));

  switch(b.bintype) {
    case(0): //xi_ell
  if(b.logxopt == 0) {
    xbin = min((int) floor((rsep-b.minx)/b.dx),b.nx-1);
    }
  else {
    xbin = min((int) floor((log10(rsep)-b.minx)/b.dx),b.nx-1);
    }
  ybin = min(((int) floor((sqrt(mymu2)*b.ny))),b.ny-1);
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(xbin < b.nx && xbin >= 0);
  assert(mymu2 >= 0. && mymu2 <= 1.);
  assert(ybin < b.ny && ybin >= 0);
  #endif
/*
  if(xbin == 0) {
    printf("%e %e %e %e %e %e %e %e %e %e %Le %d %d\n",pos1[0],pos2[0],rvec[0],pos1[1],pos2[1],rvec[1],pos1[2],pos2[2],rvec[2],rsep,mymu2,xbin,ybin);
    }
*/
  break;
    case(1): //xi_grid
  rpi = rsep*sqrt(mymu2);
  if(b.logyopt == 1) {
    if(rpi <= 0.) {
      (*bin2d) = -1;
      return 0;
      }
    rpi = log10(rpi);
    }
  if(rpi >= b.miny + b.ny*b.dy || rpi < b.miny) {
    (*bin2d) = -1;
    return 0;
    }

  rsig = rsep*sqrt(1.-(mymu2));
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
  break; //end of xigrid case.
    default:
  exit(1);
  }  //end switch bintype.

  (*bin2d) = xbin*b.ny + ybin;
  return 1;
  } //end addpairsky3d


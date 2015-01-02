#include "header.h"


//sims output stuff.
int xiellperiodic(long double *Npairsfinal, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi) {
  int i,j,ell;
  long double rmin, rmax, rcen;
  double muval;
  long double *xis;
  long double mynvoltimesNdmu;
  xis = (long double *) malloc(sizeof(long double)*(b.ellmaxdata/2+1));
  //i checked by doing elliptical integrals -- this scaling is exact, doesn't matter what the geometry of the bins is (log/linear r, mu or xigrid, etc.)
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar;

  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
//    assert(fabs(b.dy-1.) < 1.0e-8); //better have set b.dy = 1.
    }
  long double xival;
  for(i=0;i<b.nx;i++) {
    //sum over mu bins; do more intelligent integration later?
    if(b.logxopt == 1) {
      rmin = ((long double) pow(10.,b.minx+i*b.dx));
      rmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      rmin = ((long double) b.minx+i*b.dx);
      rmax = ((long double) b.minx+(i+1)*b.dx);
      rcen = ((long double) b.minx+(i+0.5)*b.dx);
      }
    mynvoltimesNdmu = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin)*mynbar*((long double) n_gal)*((long double) b.dy);
    for(ell=0;ell<=b.ellmaxdata;ell+=2) {
      xis[ell/2] = 0.;
      }
    for(j=0;j<b.ny;j++) {
      muval = (j+0.5)*b.dy;
      if(autoorcross == 0) {
        xival = Npairsfinal[i*b.ny+j]/mynvoltimesNdmu*2.0-1.0;
        }
      else {  //no factor of 2 on pair counts for cross correlation.  We counted them all!
        xival = Npairsfinal[i*b.ny+j]/mynvoltimesNdmu-1.0;
        }
      for(ell=0;ell<=b.ellmaxdata;ell+=2) {
        xis[ell/2] += xival*legendrep(ell,muval);
        }
      }
    for(ell=0;ell<=b.ellmaxdata;ell+=2) {
      xis[ell/2] = xis[ell/2]*b.dy*(2.*ell+1.);
      xi[ell/2][i] = xis[ell/2];
//      printf("%Le %Le %Le %Le %Le %e\n",rmin,rcen,rmax,mynvoltimesNdmu,xival,xi[ell/2][i]);
      }
    }
  free(xis);
  return 0;
  } 

int vstatsperiodic(long double *Npairsfinal, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **vstats) {
  assert(b.ny == 4); //xi, vr, sig2par, sig2perp. 
  assert(fabs(b.dy - 1.) < 1.0e-6);
  int i,j;
  long double rmin, rmax, rcen;
  long double mynvoltimesNdmu;
  long double vrmean, v2parmean, v2perpmean;

  //i checked by doing elliptical integrals -- this scaling is exact, doesn't matter what the geometry of the bins is (log/linear r, mu or xigrid, etc.)
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar;

  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
    assert(fabs(b.dy-1.) < 1.0e-8); //better have set b.dy = 1.
    }
  for(i=0;i<b.nx;i++) {
    //sum over mu bins; do more intelligent integration later?
    if(b.logxopt == 1) {
      rmin = ((long double) pow(10.,b.minx+i*b.dx));
      rmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      rmin = ((long double) b.minx+i*b.dx);
      rmax = ((long double) b.minx+(i+1)*b.dx);
      rcen = ((long double) b.minx+(i+0.5)*b.dx);
      }
    mynvoltimesNdmu = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin)*mynbar*((long double) n_gal)*((long double) b.dy);
    if(autoorcross == 0) {
      vstats[0][i] = Npairsfinal[i*b.ny]/mynvoltimesNdmu*2.0-1.0;
      }
    else {
      vstats[0][i] = Npairsfinal[i*b.ny]/mynvoltimesNdmu-1.0;
      }
    for(j=1;j<=3;j++) {
      vstats[j][i] = Npairsfinal[i*b.ny+j]/Npairsfinal[i*b.ny];
      }
    //subtract the mean infall from the par dispersion, and convert the perp dispersion to 1d.
    vstats[2][i] = vstats[2][i] - vstats[1][i]*vstats[1][i];
    vstats[3][i] = vstats[3][i]*0.5;
    } //end i loop over x bins.
  return 0;
  } 


/* Same as xiellperiodic but also does some rebinning into some
 * funky bins used for the actual data */
int xiellperiodicrebin(long double *Npairsfinal,xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi) {
  int i,j,ell;
  long double rmin, rmax, rcen;
  double muval;
  long double *xis; 
  long double *Npsum;
  long double mynvoltimesNdmu;
  xis = (long double *) malloc(sizeof(long double)*(b.ellmaxdata/2+1));
  Npsum = (long double *) malloc(sizeof(long double)*(b.ny));
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar;

  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
    }

  long double xival;
  int bi;
  int cumi;
  int nmaxlist[10];
  for(bi=0;bi<10;bi++) {
    nmaxlist[bi] = 0;
    }
  nmaxlist[0] = b.nr0;
  nmaxlist[2] = b.nr2;
  nmaxlist[4] = b.nr4;
  int cumbi = 0;
  for(ell=0;ell<=b.ellmaxdata;ell+=2) { //here, let's do this!
    cumi = 0;
    for(bi=0;bi<nmaxlist[ell];bi++) {
      if(b.logxopt == 1) {
        rmin = ((long double) pow(10.,b.minx+(cumi)*b.dx));
        rmax = ((long double) pow(10.,b.minx+(cumi+b.binxiell[bi])*b.dx));
        rcen = ((long double) pow(10.,b.minx+(cumi+b.binxiell[bi]/2.0)*b.dx));
        }
      else {
        rmin = ((long double) b.minx+(cumi)*b.dx);
        rmax = ((long double) b.minx+(cumi+b.binxiell[bi])*b.dx);
        rcen = ((long double) b.minx+(cumi+b.binxiell[bi]/2.0)*b.dx);
        }
      mynvoltimesNdmu = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin)*mynbar*((long double) n_gal)*((long double) b.dy);
      xis[ell/2] = 0.;
      for(j=0;j<b.ny;j++) {
        Npsum[j] = 0;
        muval = (j+0.5)*b.dy;
        for(i=cumi;i<cumi+b.binxiell[bi];i++) {
          Npsum[j] += Npairsfinal[i*b.ny+j];
          }
        if(autoorcross == 0) {
          xival = Npsum[j]/mynvoltimesNdmu*2.0-1.0;
          }
        else {  //no factor of 2 on pair counts for cross correlation.  We counted them all!
          xival = Npsum[j]/mynvoltimesNdmu-1.0;
          }
        //tmp AP testing.
/*
        printf("%d %d %Le %e %Le %Le %Le %e\n",ell,j,rcen,muval,mynvoltimesNdmu,Npsum[j],xival,legendrep(ell,muval));
*/
        xis[ell/2] += xival*legendrep(ell,muval);
        }
      xis[ell/2] = xis[ell/2]*b.dy*(2.*ell+1.);
      (*xi)[cumbi+bi] = xis[ell/2];
      printf("beth: xi %d %e\n",cumbi+bi,(*xi)[cumbi+bi]);
      cumi = cumi+b.binxiell[bi];
      }
    cumbi += nmaxlist[ell];
    }

  free(xis);
  free(Npsum);
  return 0;
  } 

/* Function to convert the countpairs binning of Np to the xiell
 * binning (this reduces total # of numbers we're storing). The
 * result is stored in Npout. Note that Npout is still in counts, and
 * is not normalized to xi.*/
int xiellNprebin(long double *Npairsfinal, xibindat b, long double **Npout) {
  int i,j,ell;
  long double dy = 1./((long double) b.ny);
  assert(fabs(dy-b.dy) < 2.0e-6);
  long double muval;
  long double *Npell; 
  long double *Npsum;
  Npell = (long double *) malloc(sizeof(long double)*(b.ellmaxdata/2+1));
  Npsum = (long double *) malloc(sizeof(long double)*(b.ny));
  
  int bi;
  int cumi;
  int nmaxlist[10];
  for(bi=0;bi<10;bi++) {
    nmaxlist[bi] = 0;
  }
  nmaxlist[0] = b.nr0;
  nmaxlist[2] = b.nr2;
  nmaxlist[4] = b.nr4;
  int cumbi = 0;
  for(ell=0;ell<=b.ellmaxdata;ell+=2) { //here, let's do this!
    cumi = 0;
    for(bi=0;bi<nmaxlist[ell];bi++) {
      Npell[ell/2] = 0.;
      for(j=0;j<b.ny;j++) {
        Npsum[j] = 0;
        muval = ((long double) (j+0.5))*dy;
        for(i=cumi;i<cumi+b.binxiell[bi];i++) {
          Npsum[j] += Npairsfinal[i*b.ny+j];
	}
        Npell[ell/2] += Npsum[j]*longlegendrep(ell,muval);
      }
      Npell[ell/2] = Npell[ell/2]*dy*(2.*ell+1.);
      (*Npout)[cumbi+bi] = Npell[ell/2];
      cumi = cumi+b.binxiell[bi];
    }
    cumbi += nmaxlist[ell];
  }
  free(Npell);
  free(Npsum);
  return 0;
} 

void printxiell(char *outfname, xibindat b, double **xi) {
  FILE *ofp;
  ofp = open_file_write(outfname);
  int i,ell;
  long double rmin, rmax, rcen,myvol;
  for(i=0;i<b.nx;i++) {
    if(b.logxopt == 1) {
      rmin = ((long double) pow(10.,b.minx+i*b.dx));
      rmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      rmin = ((long double) (b.minx+i*b.dx));
      rmax = ((long double) (b.minx+(i+1)*b.dx));
      rcen = ((long double) (b.minx+(i+0.5)*b.dx));
      }
    myvol = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin);
    //which do i want?
    //fprintf(ofp,"%le",((double) rcen),myvol);
    fprintf(ofp,"%le",((double) rcen));
    for(ell=0;ell<=b.ellmaxdata;ell+=2) {
      fprintf(ofp," %e",xi[ell/2][i]);
      }
    fprintf(ofp,"\n");
    }  
  fclose(ofp);
  }

void printxiellrebin(char *outfname, xibindat b, double *r0D, double *r2D, double *r4D, double *xi) {
  FILE *ofp;
  ofp = open_file_write(outfname);
  int i,ell;
  long double rmin, rmax, rcen,myvol;
  fprintf(ofp,"# nr0 = %d\n",b.nr0);
  fprintf(ofp,"# nr2 = %d\n",b.nr2);
  fprintf(ofp,"# nr4 = %d\n",b.nr4);
  int cumi = 0;
  for(i=0;i<b.nr0;i++) {
    fprintf(ofp,"%e %e\n",r0D[i],xi[cumi+i]);
    }
  cumi += b.nr0;
  for(i=0;i<b.nr2;i++) {
    fprintf(ofp,"%e %e\n",r2D[i],xi[cumi+i]);
    }
  cumi += b.nr2;
  for(i=0;i<b.nr4;i++) {
    fprintf(ofp,"%e %e\n",r4D[i],xi[cumi+i]);
    }
  cumi += b.nr4;
  fclose(ofp);
  }

void printxiellwvol(char *outfname, xibindat b, double **xi) {
  FILE *ofp;
  ofp = open_file_write(outfname);
  int i,ell;
  long double rmin, rmax, rcen,myvol;
  for(i=0;i<b.nx;i++) {
    if(b.logxopt == 1) {
      rmin = ((long double) pow(10.,b.minx+i*b.dx));
      rmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      rmin = ((long double) (b.minx+i*b.dx));
      rmax = ((long double) (b.minx+(i+1)*b.dx));
      rcen = ((long double) (b.minx+(i+0.5)*b.dx));
      }
    myvol = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin);
    fprintf(ofp,"%le %.9Le",((double) rcen),myvol);
    for(ell=0;ell<=b.ellmaxdata;ell+=2) {
      fprintf(ofp," %e",xi[ell/2][i]);
      }
    fprintf(ofp,"\n");
    }  
  fclose(ofp);
  }

//start functions to cut out the small-mu bins.
int xiellperiodicrebincutsmallscale(long double *Npairsfinal,xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi) {
//still need to set up smallscalemask.
//int xiellperiodicrebincutsmallscale(long double *Npairsfinal,xibindat b, int n_gal, int n2, int autoorcross, double **xi) {
  int i,j,ell;
  long double rmin, rmax, rcen;
  long double rminbin, rmaxbin;
  double muval;
  long double *xis; 
  long double *Npsum;
  long double mynvoltimesNdmu,mynvoltimesNdmuuncut;
  xis = (long double *) malloc(sizeof(long double)*(b.ellmaxdata/2+1));
  Npsum = (long double *) malloc(sizeof(long double)*(b.ny));
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar;

  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
    }

  long double xival;
  int bi;
  int cumi;
  int nmaxlist[10];
  for(bi=0;bi<10;bi++) {
    nmaxlist[bi] = 0;
    }
  nmaxlist[0] = b.nr0;
  nmaxlist[2] = b.nr2;
  nmaxlist[4] = b.nr4;
  int cumbi = 0;
  int nobinflag;
  

  for(ell=0;ell<=b.ellmaxdata;ell+=2) { //here, let's do this!
    cumi = 0;
    for(bi=0;bi<nmaxlist[ell];bi++) {
      if(b.logxopt == 1) {
        rmin = ((long double) pow(10.,b.minx+(cumi)*b.dx));
        rmax = ((long double) pow(10.,b.minx+(cumi+b.binxiell[bi])*b.dx));
        rcen = ((long double) pow(10.,b.minx+(cumi+b.binxiell[bi]/2.0)*b.dx));
        }
      else {
        rmin = ((long double) b.minx+(cumi)*b.dx);
        rmax = ((long double) b.minx+(cumi+b.binxiell[bi])*b.dx);
        rcen = ((long double) b.minx+(cumi+b.binxiell[bi]/2.0)*b.dx);
        }
      mynvoltimesNdmuuncut = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin)*mynbar*((long double) n_gal)*((long double) b.dy);
      xis[ell/2] = 0.;
      for(j=0;j<b.ny;j++) {
        Npsum[j] = 0;
        mynvoltimesNdmu = 0.;
        muval = (j+0.5)*b.dy;
        nobinflag = 1; //check if any bins are included; if not, set xi to 0.
        for(i=cumi;i<cumi+b.binxiell[bi];i++) {
          if(b.logxopt == 1) {
            rminbin = ((long double) pow(10.,b.minx+(i)*b.dx));
            rmaxbin = ((long double) pow(10.,b.minx+(i+1)*b.dx));
            }
          else {
            rminbin = ((long double) b.minx+(i)*b.dx);
            rmaxbin = ((long double) b.minx+(i+1)*b.dx);
            }

          if(b.mask[i*b.ny+j] == 0) {
            nobinflag = 0;
            Npsum[j] += Npairsfinal[i*b.ny+j];
            mynvoltimesNdmu += 4./3.*M_PI*(rmaxbin*rmaxbin*rmaxbin-rminbin*rminbin*rminbin)*mynbar*((long double) n_gal)*((long double) b.dy);
            }
          }
        if(nobinflag == 1) {
          xival = 0.;
          }
        else {
          if(autoorcross == 0) {
            //this is temporary cross-check -- passed with no mask case.
/*
            printf("%Le %Le %Le\n",mynvoltimesNdmu,mynvoltimesNdmuuncut,fabs(mynvoltimesNdmu - mynvoltimesNdmuuncut)/mynvoltimesNdmuuncut);
            assert(fabs(mynvoltimesNdmu - mynvoltimesNdmuuncut)/mynvoltimesNdmuuncut < 2.0e-6);
*/
            xival = Npsum[j]/mynvoltimesNdmu*2.0-1.0;
            }
          else {  //no factor of 2 on pair counts for cross correlation.  We counted them all!
            xival = Npsum[j]/mynvoltimesNdmu-1.0;
            }
          }
        xis[ell/2] += xival*legendrep(ell,muval);
        }
      xis[ell/2] = xis[ell/2]*b.dy*(2.*ell+1.);
      (*xi)[cumbi+bi] = xis[ell/2];
      cumi = cumi+b.binxiell[bi];
      }
    cumbi += nmaxlist[ell];
    }

  free(xis);
  free(Npsum);
  return 0;
  } 

/* Function to convert the countpairs binning of Np to the xiell
 * binning (this reduces total # of numbers we're storing). The
 * result is stored in Npout. Note that Npout is still in counts, and
 * is not normalized to xi.*/
int xiellNprebincutsmallscale(long double *Npairsfinal, xibindat b, long double **Npout) {
  printf("this function needs to be written!\n");
  exit(1);
  int i,j,ell;
  long double dy = 1./((long double) b.ny);
  assert(fabs(dy-b.dy) < 2.0e-6);
  long double muval;
  long double *Npell; 
  long double *Npsum;
  Npell = (long double *) malloc(sizeof(long double)*(b.ellmaxdata/2+1));
  Npsum = (long double *) malloc(sizeof(long double)*(b.ny));
  
  int bi;
  int cumi;
  int nmaxlist[10];
  for(bi=0;bi<10;bi++) {
    nmaxlist[bi] = 0;
  }
  nmaxlist[0] = b.nr0;
  nmaxlist[2] = b.nr2;
  nmaxlist[4] = b.nr4;
  int cumbi = 0;
  for(ell=0;ell<=b.ellmaxdata;ell+=2) { //here, let's do this!
    cumi = 0;
    for(bi=0;bi<nmaxlist[ell];bi++) {
      Npell[ell/2] = 0.;
      for(j=0;j<b.ny;j++) {
        Npsum[j] = 0;
        muval = ((long double) (j+0.5))*dy;
        for(i=cumi;i<cumi+b.binxiell[bi];i++) {
          Npsum[j] += Npairsfinal[i*b.ny+j];
	}
        Npell[ell/2] += Npsum[j]*longlegendrep(ell,muval);
      }
      Npell[ell/2] = Npell[ell/2]*dy*(2.*ell+1.);
      (*Npout)[cumbi+bi] = Npell[ell/2];
      cumi = cumi+b.binxiell[bi];
    }
    cumbi += nmaxlist[ell];
  }
  free(Npell);
  free(Npsum);
  return 0;
}

int xigridperiodic(long double *Npairsfinal, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi) {
  int i,j;
  long double rsigmin, rsigmax, rsigcen, rpimin, rpimax, rpicen, mynvoltimesN;
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar;
  long double xival;

  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
    }
  
  for(i=0;i<b.nx;i++) {
    if(b.logxopt == 1) {
      rsigmin = ((long double) pow(10.,b.minx+i*b.dx));
      rsigmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rsigcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      // Bins in Rp
      rsigmin = ((long double) b.minx+i*b.dx);
      rsigmax = ((long double) b.minx+(i+1)*b.dx);
      rsigcen = ((long double) b.minx+(i+0.5)*b.dx);
      }
    for(j=0;j<b.ny;j++) {
      if(b.logyopt == 1) {
        rpimin = ((long double) pow(10.,b.miny+j*b.dy));
        rpimax = ((long double) pow(10.,b.miny+(j+1)*b.dy));
        rpicen = ((long double) pow(10.,b.miny+(j+0.5)*b.dy));
        }
      else {
        // Bins in Rp
        rpimin = ((long double) b.miny+j*b.dy);
        rpimax = ((long double) b.miny+(j+1)*b.dy);
        rpicen = ((long double) b.miny+(j+0.5)*b.dy);
        }



      // Volume of annulus*height
      mynvoltimesN = M_PI*(rsigmax*rsigmax-rsigmin*rsigmin)*mynbar*((long double) n_gal)*((long double) (2.*(rpimax-rpimin)));
      if(autoorcross == 0) {
        xival = (Npairsfinal[i*b.ny+j]/mynvoltimesN)*2.0-1.0;
        }
      else {
        xival = (Npairsfinal[i*b.ny+j]/mynvoltimesN)-1.0;
        }
      xi[i][j] = ((double) (xival));
      } //end j loop
    } //end i loop.
  return 0;
  } //end xigridperiodic.

void printxigrid(char *outfname, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi) {
  int i,j;
  long double rsigmin, rsigmax, rsigcen, rpimin, rpimax, rpicen, mynvoltimesN;
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar; 

  FILE *ofp;
  ofp = open_file_write(outfname);
  fprintf(ofp,"# nbar1 = %.6Le\n",((long double) n_gal)/myboxvol);
  if(autoorcross == 0) {
    fprintf(ofp,"# nbar2 = %.6Le\n",((long double) 0.0));
    }
  else {
    fprintf(ofp,"# nbar2 = %.6Le\n",((long double) n2)/myboxvol);
    }
  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
    }

  for(i=0;i<b.nx;i++) {
    if(b.logxopt == 1) {
      rsigmin = ((long double) pow(10.,b.minx+i*b.dx));
      rsigmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rsigcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      // Bins in Rp
      rsigmin = ((long double) b.minx+i*b.dx);
      rsigmax = ((long double) b.minx+(i+1)*b.dx);
      rsigcen = ((long double) b.minx+(i+0.5)*b.dx);
      }
    for(j=0;j<b.ny;j++) {
      if(b.logyopt == 1) {
        rpimin = ((long double) pow(10.,b.miny+j*b.dy));
        rpimax = ((long double) pow(10.,b.miny+(j+1)*b.dy));
        rpicen = ((long double) pow(10.,b.miny+(j+0.5)*b.dy));
        }
      else {
        // Bins in Rp
        rpimin = ((long double) b.miny+j*b.dy);
        rpimax = ((long double) b.miny+(j+1)*b.dy);
        rpicen = ((long double) b.miny+(j+0.5)*b.dy);
        }
      // Volume of annulus*height
      mynvoltimesN = M_PI*(rsigmax*rsigmax-rsigmin*rsigmin)*mynbar*((long double) n_gal)*((long double) (2.*(rpimax-rpimin)));
      fprintf(ofp,"%le %le %le\n",((double) rsigcen), ((double) rpicen), xi[i][j]);
      } //end j loop
    } //end i loop.
  fclose(ofp);
  }

void printxigridwNp(char *outfname, char *outfnameNp, xibindat b, int n_gal, int n2, int autoorcross, real Lbox, real APperp, real APpar, double **xi, long double *Npairsfinal) {
  int i,j;
  long double rsigmin, rsigmax, rsigcen, rpimin, rpimax, rpicen, mynvoltimesN;
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar; 

  FILE *ofp;
  ofp = open_file_write(outfname);
  FILE *ofpNp;
  ofpNp = open_file_write(outfnameNp);

  fprintf(ofp,"# nbar1 = %.6Le\n",((long double) n_gal)/myboxvol);
  fprintf(ofpNp,"# nbar1 = %.6Le\n",((long double) n_gal)/myboxvol);
  if(autoorcross == 0) {
    fprintf(ofp,"# nbar2 = %.6Le\n",((long double) 0.0));
    fprintf(ofpNp,"# nbar2 = %.6Le\n",((long double) 0.0));
    }
  else {
    fprintf(ofp,"# nbar2 = %.6Le\n",((long double) n2)/myboxvol);
    fprintf(ofpNp,"# nbar2 = %.6Le\n",((long double) n2)/myboxvol);
    }

  if(autoorcross == 0) {
    mynbar = ((long double) n_gal)/myboxvol;
    }
  else {
    mynbar = ((long double) n2)/myboxvol;
    }

  for(i=0;i<b.nx;i++) {
    if(b.logxopt == 1) {
      rsigmin = ((long double) pow(10.,b.minx+i*b.dx));
      rsigmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rsigcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
      }
    else {
      // Bins in Rp
      rsigmin = ((long double) b.minx+i*b.dx);
      rsigmax = ((long double) b.minx+(i+1)*b.dx);
      rsigcen = ((long double) b.minx+(i+0.5)*b.dx);
      }
    for(j=0;j<b.ny;j++) {
      if(b.logyopt == 1) {
        rpimin = ((long double) pow(10.,b.miny+j*b.dy));
        rpimax = ((long double) pow(10.,b.miny+(j+1)*b.dy));
        rpicen = ((long double) pow(10.,b.miny+(j+0.5)*b.dy));
        }
      else {
        // Bins in Rp
        rpimin = ((long double) b.miny+j*b.dy);
        rpimax = ((long double) b.miny+(j+1)*b.dy);
        rpicen = ((long double) b.miny+(j+0.5)*b.dy);
        }
      // Volume of annulus*height
      mynvoltimesN = M_PI*(rsigmax*rsigmax-rsigmin*rsigmin)*mynbar*((long double) n_gal)*((long double) (2.*(rpimax-rpimin)));
      fprintf(ofp,"%le %le %le\n",((double) rsigcen), ((double) rpicen), xi[i][j]);
      if(autoorcross == 0) {
          fprintf(ofpNp,"%le %le %.12Le %.12Le\n",((double) rsigcen), ((double) rpicen), mynvoltimesN, Npairsfinal[i*b.ny+j]*2.0);
          }
      else {
        fprintf(ofpNp,"%le %le %.12Le %.12Le\n",((double) rsigcen), ((double) rpicen), mynvoltimesN, Npairsfinal[i*b.ny+j]);
        }
      } //end j loop
    } //end i loop.
  fclose(ofp);
  fclose(ofpNp);
  }


/* Function to integrate counts over LOS axis and return wp. 
 * wp = integral xi(rp, rpi) dy - some constant.
 * xi = Npairs x some factor.
 * The result is filled into wp. Note that wp is already
 * allocated with nx items */
int wpperiodic(long double *Npairsfinal, xibindat b, int n_gal, real Lbox, real APperp, real APpar, double *wp) {
  int i,j;
  long double rsigmin, rsigmax, rsigcen,rpiedge,mynvoltimesN;
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar = ((long double) n_gal)/myboxvol;
  long double xival;
  
  //make sure y-binning correct using my assumptions; relax later?
  assert(fabs(b.miny) < 2.0e-6);
  assert(b.miny+b.ny*b.dy >= b.rpimax);
  assert(b.logyopt == 0);
  for(i=0;i<b.nx;i++) {
    //sum over mu bins; do more intelligent integration later?
    if(b.logxopt == 1) {
      rsigmin = ((long double) pow(10.,b.minx+i*b.dx));
      rsigmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
      rsigcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
    }
    else {
      // Bins in Rp
      rsigmin = ((long double) b.minx+i*b.dx);
      rsigmax = ((long double) b.minx+(i+1)*b.dx);
      rsigcen = ((long double) b.minx+(i+0.5)*b.dx);
    }
    wp[i] = 0.;
    for(j=0;j<b.ny;j++) {
      rpiedge = ((long double) b.miny+(j+1)*b.dy);
      // Integrate out to rpimax.
      if(rpiedge <= b.rpimax*1.0001) {  //keep this bin in the calculation
	// Volume of annulus*height
        mynvoltimesN = M_PI*(rsigmax*rsigmax-rsigmin*rsigmin)*mynbar*((long double) n_gal)*((long double) (2.*b.dy));
        xival = (Npairsfinal[i*b.ny+j]/mynvoltimesN)*2.0-1.0;
        wp[i] += xival;
      }
    }
    wp[i] = wp[i]*b.dy*2.; // Factor of 2 to go from -pimax to pimax instead of 0 to pimax
  }
  return 0;
} 

/* Function to convert the countpairs binning of Np to the wp
 * binning. Do the integration over the LOS axis, and only store the
 * ndata rp bins for each mass bin pair.Don't normalize to wp */
int wpNprebin(long double *Npairsfinal, xibindat b, long double **Npout) {
  int i,j;
  long double rpiedge;  
  long double Npsum;  

  //make sure y-binning correct using my assumptions; relax later?
  assert(fabs(b.miny) < 2.0e-6);
  assert(b.miny+b.ny*b.dy >= b.rpimax);
  assert(b.logyopt == 0);

  // Integrate over LOS
  for(i=0;i<b.nx;i++) {

    //sum over bins; do more intelligent integration later?
    Npsum = 0.;
    for(j=0;j<b.ny;j++) {
      rpiedge = ((long double) b.miny+(j+1)*b.dy);
      if(rpiedge <= b.rpimax*1.0001) {  //keep this bin in the calculation
        Npsum += Npairsfinal[i*b.ny+j]; // sum over bins along LOS
      }
    } // end j
    (*Npout)[i] = Npsum;

  } // end i

  return 0;
} 

void printwp(char *outfname, xibindat b, double *wp) {
  FILE *ofp;
  ofp = open_file_write(outfname);
  int i;
  for(i=0;i<b.nx;i++) {

    if(b.logxopt == 1) {
      fprintf(ofp,"%e %e\n",pow(10.,b.minx+(i+0.5)*b.dx),wp[i]);
      fflush(stdout);
    }
    else {
      fprintf(ofp,"%e %e\n",b.minx+(i+0.5)*b.dx,wp[i]);
      fflush(stdout);
    }
  }
  fclose(ofp);
}

//end sims output stuff.

void printNpairsgeneric(char *foutbase, long double *Npairsfinal,xibindat b,int DRopt,int n1, long double n1wgt, long double n2wgt,char *binfname, double omfid, double hfid) {

  double xcen,ycen;
  char outfnameNp[MAXLINELEN];
  int i,j;
  if(DRopt == -1) {
    sprintf(outfnameNp,"%s.Np",foutbase);
    }
  else {
    sprintf(outfnameNp,"%s.Np.DRopt%d",foutbase,DRopt);
    }

  if(b.bintype == 3) { //append with .ang to avoid overwriting 3d xi counts.
    sprintf(outfnameNp,"%s.ang",outfnameNp);
    }

  FILE *ofpNp;
  ofpNp = open_file_write(outfnameNp);
  if(DRopt != 3) {
    fprintf(ofpNp,"# nobj: %d\n",n1);
    fprintf(ofpNp,"# weight sums: %.12Le, %.12Le\n",n1wgt,n2wgt);
    fprintf(ofpNp,"# binfile: %s\n",binfname);
    }
  else {
    fprintf(ofpNp,"# nobj: %d\n",n1);
    fprintf(ofpNp,"# weight sums: %.12Le, %.12Le\n",n2wgt,n1wgt);
    fprintf(ofpNp,"# binfile: %s\n",binfname);
    }
  fprintf(ofpNp,"# omfid: %f\n", omfid);
  fprintf(ofpNp,"# hfid: %f\n", hfid);

  for(i=0;i<b.nx;i++) {
    if(b.logxopt == 0) {
      xcen = b.minx + (i+0.5)*b.dx;
      }
    else {
      xcen = pow(10.,b.minx + (i+0.5)*b.dx);
      }
    for(j=0;j<b.ny;j++) {
      if(b.logyopt == 0) {
        ycen = b.miny + (j+0.5)*b.dy;
        }
      else {
        ycen = pow(10.,b.miny + (j+0.5)*b.dy);
        }
      if(b.ny > 1) {
        fprintf(ofpNp,"%e %e %.12Le\n",xcen,ycen,Npairsfinal[i*b.ny+j]);
        }
      else {
        fprintf(ofpNp,"%e %.12Le\n",xcen,Npairsfinal[i*b.ny+j]);
        }

      }  //end j loop.
    } //end loop over b.nx
  fclose(ofpNp);
  } //end printNpairsradecz


void printNpairssim(int n1, int n2, int autoorcross, xibindat b, real Lbox, real APperp, real APpar, char *foutbase, long double *Npairsfinal) {

  char outfname[MAXLINELEN];
  char outfnameNp[MAXLINELEN];

  int i,j,ell;

  sprintf(outfnameNp,"%s.Np",foutbase);
  if(b.bintype == 0) {
    sprintf(outfname,"%s.xiell",foutbase);
    }
  if(b.bintype == 1) {
    sprintf(outfname,"%s.xigrid",foutbase);
    }
  if(b.bintype == 2) {
    sprintf(outfname,"%s.xiv",foutbase);
    }

  FILE *ofp, *ofpNp;

  double **xi;
  double *xiellrebin;
  long double rmin, rmax, rcen;
  long double APrescale = ((long double) (APperp*APperp*APpar));
  long double myboxvol = ((long double) Lbox)*((long double) Lbox)*((long double) Lbox)/APrescale;
  long double mynbar;
  long double mynbar2;
  long double mynvoltimesNdmu;

  mynbar = ((long double) n1)/myboxvol;
  mynbar2 = ((long double) n2)/myboxvol;

  double muval;
  int myn2 = n2;
  if(myn2 < 0) {
    myn2 = n1;
    }
  if(b.rebinopt == 0) {
    if(b.bintype == 0 || b.bintype == 2) {
      ofp = open_file_write(outfname);
      ofpNp = open_file_write(outfnameNp);
      if(b.bintype == 0) {
        xi = malloc2ddouble(b.ellmaxdata/2+1,b.nx);
        }
      else {
        xi = malloc2ddouble(b.ny,b.nx);
        }  
      for(i=0;i<b.nx;i++) {
        //sum over mu bins; do more intelligent integration later?
        if(b.logxopt == 1) {
          rmin = ((long double) pow(10.,b.minx+i*b.dx));
          rmax = ((long double) pow(10.,b.minx+(i+1)*b.dx));
          rcen = ((long double) pow(10.,b.minx+(i+0.5)*b.dx));
          }
        else {
          rmin = ((long double) b.minx+i*b.dx);
          rmax = ((long double) b.minx+(i+1)*b.dx);
          rcen = ((long double) b.minx+(i+0.5)*b.dx);
          }
        mynvoltimesNdmu = 4./3.*M_PI*(rmax*rmax*rmax-rmin*rmin*rmin)*mynbar*((long double) myn2)*((long double) b.dy);
        for(j=0;j<b.ny;j++) {
          muval = (j+0.5)*b.dy;
          if(autoorcross == 0) {
            fprintf(ofpNp,"%Le %le %.12Le %.12Le\n",rcen,muval,mynvoltimesNdmu,Npairsfinal[i*b.ny+j]*2.0);
            }
          else {
            fprintf(ofpNp,"%Le %le %.12Le %.12Le\n",rcen,muval,mynvoltimesNdmu,Npairsfinal[i*b.ny+j]);
            }
          }
        } //end Np print loop over nx.
      fclose(ofpNp);
      if(b.bintype == 0) {
        xiellperiodic(Npairsfinal,b,n1,n2,autoorcross,Lbox,APperp,APpar,xi);
        }
      else {
        vstatsperiodic(Npairsfinal,b,n1,n2,autoorcross,Lbox,APperp,APpar,xi); //xi here is really xiand velocity stats.
        }
// print nbar at file header for simplicity.
      fprintf(ofp,"# nbar1 = %.6Le\n",mynbar);
      assert(autoorcross == 0);
      if(autoorcross == 0) {
        fprintf(ofp,"# nbar2 = %.6Le\n",((long double) 0.0));
        }
      else {
        fprintf(ofp,"# nbar2 = %.6Le\n",((long double) n2)/myboxvol);
        }

      for(i=0;i<b.nx;i++) {
        if(b.logxopt == 1) {
          fprintf(ofp,"%e",pow(10.,b.minx+(i+0.5)*b.dx));
          }
        else {
          fprintf(ofp,"%e",b.minx+(i+0.5)*b.dx);
          }
        if(b.bintype == 0) {
          for(ell=0;ell<=b.ellmaxdata;ell+=2) {
            fprintf(ofp," %e",xi[ell/2][i]);
            }
          fprintf(ofp,"\n");
          } //bintype = 0
        if(b.bintype == 2) {
          for(j=0;j<b.ny;j++) {
            fprintf(ofp," %e",xi[j][i]);
            }
          fprintf(ofp,"\n");
          } //bintype = 2
        } //end loop over x bins.
      malloc2dfree(xi);
      fclose(ofp);
      } //end bintype = 0/2 print xi
    if(b.bintype == 1) {  //bintype == 1; this is used for butterfly plots.
      xi = malloc2ddouble(b.nx, b.ny);
      xigridperiodic(Npairsfinal,b,n1,0,0,Lbox,APperp,APpar,xi);
      printf("these are fnames %s %s\n",outfname,outfnameNp);
      printxigridwNp(outfname,outfnameNp,b, n1,0,0, Lbox, APperp, APpar, xi,Npairsfinal);
      }
    }
  else { //rebinopt == 1!!
    fprintf(stderr,"nevermind, this is not coded yet!\n");
/*
    xiellrebin = (double *) malloc(sizeof(double)*(b.ndata));
    if(b.maskopt == 0) {
      xiellperiodicrebin(Npairsfinal,b,n1,0,0,Lbox,APperp,APpar,&xiellrebin);
      printxiellrebin(outfname,b,r0D,r2D,r4D,xiellrebin);
      }
    else {  //mask file!
      xiellperiodicrebincutsmallscale(Npairsfinal,b,n1,0,0,&xiellrebin);
      printxiellrebin(outfname,b,r0D,r2D,r4D,xiellrebin);
      }
    free(xiellrebin);
*/
    }
  } //end printNpairssim

#include "header.h"

extern gsl_rng *rngen;

//ra and dec input must be in radians already!!!
//not shifting by originpos anymore!!!

void radecztopos(double ra, double dec, float z, int unitsMpc, int angopt, cosmo_params cosmopfid, real pos[3], real *chi) {
  if(angopt == 1) {
    (*chi) = 1.;
/*
    (pos)[0] = cos(dec)*cos(ra)+originpos[0];
    (pos)[1] = cos(dec)*sin(ra)+originpos[1];
    (pos)[2] = sin(dec)+originpos[2];
*/
    (pos)[0] = cos(dec)*cos(ra);
    (pos)[1] = cos(dec)*sin(ra);
    (pos)[2] = sin(dec);

    } //angular clustering option.
  else {
    if(unitsMpc == 1) {
      (*chi) = comoving_distMpc_lcdm((double) z, cosmopfid);
      }
    else {
      (*chi) = comoving_disthinvMpc_lcdm((double) z, cosmopfid);
      }
/*
    (pos)[0] = (*chi)*cos(dec)*cos(ra)+originpos[0];
    (pos)[1] = (*chi)*cos(dec)*sin(ra)+originpos[1];
    (pos)[2] = (*chi)*sin(dec)+originpos[2];
*/
    (pos)[0] = (*chi)*cos(dec)*cos(ra);
    (pos)[1] = (*chi)*cos(dec)*sin(ra);
    (pos)[2] = (*chi)*sin(dec);
    } //end angopt = 0 (3d correlation option)
  } //end radecztopos.

float getcatzmax(char *infilename, int ftype) {
  FILE *ifp;
  if(!(ifp = open_file_read(infilename))) {
    fprintf(stderr,"Exiting!\n");
    exit(1);
    }
  int header, linestot;
  get_file_length(ifp,&header,&linestot);
  double ra, dec;
  float z,finalwgt;
  int i;
  char line[MAXLINELEN];

  float myzmax;

  z = -1.;  //value if no input z.

  for(i=0;i<header;i++) {
    fgets(line,MAXLINELEN,ifp);
    }
  for(i=0;i<linestot-header;i++) {
    switch(ftype) {
      case(1):
    fscanf(ifp,"%lf %lf %f %f\n",&ra,&dec,&z,&finalwgt);
    break;
      case(2):
    fscanf(ifp,"%lf %lf %f\n",&ra,&dec,&z);
    finalwgt = 1.;
    break;
      case(3):
    fscanf(ifp,"%lf %lf %f\n",&ra,&dec,&finalwgt);
//    z = 0.5*(zmin + zmax);
    break;
      case(4):
    fscanf(ifp,"%lf %lf\n",&ra,&dec);
//    z = 0.5*(zmin + zmax);
    finalwgt = 1.;
    break;
      default:
    fprintf(stderr,"your ftype is not supported.  Exiting!\n");
    exit(1);
    break;
    } //end ftype switch.
    if(i==0) {
      myzmax = z;
      }
    else {
      myzmax = max(z,myzmax);
      }
    } //end loop to count gals and weights in redshift range.
  return myzmax;
  }  //end getcatzmax

particle *readcat(char *infilename, int ftype, int unitsMpc, int angopt, cosmo_params cosmopfid, float zmin, float zmax, int DorR, float ndownRR, int *ntot, long double *wgttot, real *maxdist) {
  FILE *ifp;
  ifp = open_file_read(infilename);
  int header, linestot;
  get_file_length(ifp,&header,&linestot);
  double ra, dec;
  float z,finalwgt;
  (*wgttot) = 0.;
  int i,ii;
  char line[MAXLINELEN];
  int ipos[3];
  float catzmax = -1;
  particle *gallist;

  int dosubsample = 1;
  if(fabs(ndownRR - 1.0) < 1.0e-3 || ndownRR < 0.) {
    dosubsample = 0;
    }

  double ndownthresh = 1./ndownRR;
  int nRmax;

  (*ntot) = 0;
  (*wgttot) = 0.;

  for(i=0;i<header;i++) {
    fgets(line,MAXLINELEN,ifp);
    }
  for(i=0;i<linestot-header;i++) {
    switch(ftype) {
      case(1):
    fscanf(ifp,"%lf %lf %f %f\n",&ra,&dec,&z,&finalwgt);
    break;
      case(2):
    fscanf(ifp,"%lf %lf %f\n",&ra,&dec,&z);
    finalwgt = 1.;
    break;
      case(3):
    fscanf(ifp,"%lf %lf %f\n",&ra,&dec,&finalwgt);
    z = 0.5*(zmin + zmax);
    break;
      case(4):
    fscanf(ifp,"%lf %lf\n",&ra,&dec);
    z = 0.5*(zmin + zmax);
    finalwgt = 1.;
    break;
      default:
    fprintf(stderr,"your ftype is not supported.  Exiting!\n");
    exit(1);
    break;
    } //end ftype switch.
    catzmax = max(catzmax,z);
    if(z >= zmin && z <= zmax) {
      (*ntot) += 1;
      }
    } //end loop to count gals and weights in redshift range.
  rewind(ifp);

  if(dosubsample == 0) {
    gallist = (particle *) malloc(sizeof(particle)*(*ntot));
    }
  else {
    double mybuffer = 1.+30./sqrt((*ntot));
    nRmax = ((int) floor(((*ntot)/ndownRR*mybuffer))) + 1;
    gallist = (particle *) malloc(sizeof(particle)*nRmax);
    }

  if((zmax - catzmax) > 0.1*max(zmax,catzmax)) {
    fprintf(stderr,"You should update your zmax to be closer to the zmax of the input catalog, %f\n",catzmax);
    exit(1);
    }
  real posmin[3];
  real posmax[3];
  for(ii=0;ii<=2;ii++) {
    if(unitsMpc == 1) {
      posmin[ii] = comoving_distMpc_lcdm((double) zmax, cosmopfid);
      posmax[ii] = -comoving_distMpc_lcdm((double) zmax, cosmopfid);
      }
    else {
      posmin[ii] = comoving_disthinvMpc_lcdm((double) zmax, cosmopfid);
      posmax[ii] = -comoving_disthinvMpc_lcdm((double) zmax, cosmopfid);
      }
    }

  real chi;
  int gindx = 0;
  for(i=0;i<header;i++) {
    fgets(line,MAXLINELEN,ifp);
    }
  for(i=0;i<linestot-header;i++) {
switch(ftype) {
      case(1):
    fscanf(ifp,"%lf %lf %f %f\n",&ra,&dec,&z,&finalwgt);
    break;
      case(2):
    fscanf(ifp,"%lf %lf %f\n",&ra,&dec,&z);
    finalwgt = 1.;
    break;
      case(3):
    fscanf(ifp,"%lf %lf %f\n",&ra,&dec,&finalwgt);
    z = 0.5*(zmin + zmax);
    break;
      case(4):
    fscanf(ifp,"%lf %lf\n",&ra,&dec);
    z = 0.5*(zmin + zmax);
    finalwgt = 1.;
    break;
      default:
    fprintf(stderr,"your ftype is not supported.  Exiting!\n");
    exit(1);
    break;
    } //end ftype switch.

    if(z >= zmin && z <= zmax) {
      if(dosubsample == 1) {
        //draw a random uniform number.
        if(gsl_rng_uniform(rngen) > ndownthresh) {
          continue;
          }
        }
      ra = ra*M_PI/180.;
      dec = dec*M_PI/180.;
      radecztopos(ra,dec,z, unitsMpc, angopt, cosmopfid, gallist[gindx].pos, &chi);  //chi used to be saved in gallist, now its not.
//void radecztopos(double ra, double dec, float z, int unitsMpc, int angopt, cosmo_params cosmopfid, real pos[3], real *chi) {
      for(ii=0;ii<=2;ii++) {
        //ipos[ii] = ((int) floor(gallist[gindx].pos[ii]*Ncell/Lbox));
        posmin[ii] = min(posmin[ii],gallist[gindx].pos[ii]);
        posmax[ii] = max(posmax[ii],gallist[gindx].pos[ii]);
        }
      gallist[gindx].weight = finalwgt;
      (*wgttot) += finalwgt;
      //this is now set in the countpairs part.
      //gallist[gindx].DorR = DorR;
      //gallist[gindx].cell = i2n(ipos);  //going to assign this in the countpairs section.
      gindx += 1;
      if(gindx == nRmax+1) {
        printf("nR alloc error for subsample case!  exiting!!\n");
        exit(1);
        }
      }
    } //end loop over gals.
  if(dosubsample == 0) {
    assert(gindx == (*ntot));
    }
  else {
    assert(gindx < nRmax);
    #ifdef REALLYVERBOSE
    printf("downsample produced this many, this expected: %d %d %e\n",(*ntot), gindx, ((float) (*ntot))/ndownRR); 
    #endif
    (*ntot) = gindx;
    gallist = (particle *) realloc(gallist,sizeof(particle)*(*ntot));
    }
  fclose(ifp);
/*
  for(ii=0;ii<=2;ii++) {
    assert(posmin[ii] >= 0. && posmin[ii] <= Lbox);
    }
*/
  #ifdef REALLYVERBOSE
  printf("catalog mins: %e %e %e %e %e\n",comoving_distMpc_lcdm((double) zmax, cosmopfid),comoving_disthinvMpc_lcdm((double) zmax, cosmopfid),posmin[0],posmin[1],posmin[2]);
  printf("catalog maxs: %e %e %e %e %e\n",comoving_distMpc_lcdm((double) zmax, cosmopfid),comoving_disthinvMpc_lcdm((double) zmax, cosmopfid),posmax[0],posmax[1],posmax[2]);
  #endif

  //compute the maximum distance away to set Lbox for this data set.
  (*maxdist) = -1000.;
  for(ii=0;ii<=2;ii++) {
    (*maxdist) = max(fabs(posmin[ii]),(*maxdist));
    (*maxdist) = max(fabs(posmax[ii]),(*maxdist));
     }
  return gallist;
  } //end readcat

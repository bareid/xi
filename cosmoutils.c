#include "header.h"
#include <gsl/gsl_integration.h>

//DH = c/H0 in h^-1 Mpc units = 2.997e8/1e5
#define DH 2997.92458
#define OGH2FIX 2.469e-5 //got this numbe from page 14 of WMAP7 Komatsu et al.

static int cosmosetup = 0;

size_t limitI=20000;
double epsabs=1.0e-8;
double epsrel=1.0e-5;  //adjust later.
static gsl_integration_workspace *workspace;

int freecosmo() {
  if(cosmosetup == 1) {
    gsl_integration_workspace_free(workspace);
    return 0;
    }
  return 0;
  }

int docosmosetup() {
  if(cosmosetup == 1) {
    return 0;
    }

  workspace = gsl_integration_workspace_alloc(limitI);
  cosmosetup = 1;
  return 0;
  }

//stuff from olcdm.c in cosmoxi2dmcmcv3/
//H(a)/H0
double normH(double a, void *params) {
  cosmo_params *myp = (cosmo_params *) params;
  double omegam = myp -> omegam;
  double omegal = myp -> omegal;
  double omegak = myp -> omegak;
  return(sqrt(omegam/gsl_pow_3(a) + omegal + omegak/gsl_pow_2(a)));
  }

double G_integrand(double ap, void *params) {
  return(1./gsl_pow_3(ap*normH(ap,params)));
  }

double Dofa(double aval, void *params) {
  docosmosetup();
  gsl_function AA;
  AA.function = &G_integrand;
  AA.params = params;

  cosmo_params *myp = (cosmo_params *) params;
  double omegam = myp -> omegam;
  //double omegal = myp -> omegal;
  //double omegak = myp -> omegak;

  int intstatus;
  double i1, i1err;
  intstatus = gsl_integration_qag(&AA,1.0e-6,aval,epsabs,epsrel,limitI,GSL_INTEG_GAUSS31,workspace,&i1,&i1err);
  return(5./2.*omegam*normH(aval,params)*i1);
  }

double dlnDdlna(double aval, void *params) {
  cosmo_params *myp = (cosmo_params *) params;
  double omegam = myp -> omegam;
  //double omegal = myp -> omegal;
  double omegak = myp -> omegak;

  return(1./gsl_pow_2(normH(aval,params))*(-3./2.*omegam/gsl_pow_3(aval)-omegak/gsl_pow_2(aval)+5./2.*aval/Dofa(aval,params)*omegam/gsl_pow_3(aval)));
  }

double cod_integrand(double ap, void *params) {
  return(1./(ap*ap*normH(ap,params)));
  }

//units are Mpc, not h^-1 Mpc!!
double comoving_distMpc(double aval, void *params) {
  docosmosetup();
  gsl_function BB;
  BB.function = &cod_integrand;
  BB.params = params;
  int intstatus;
  double i1, i1err;

  cosmo_params *myp = (cosmo_params *) params;
  double hval = myp -> h;
  double omegak = myp -> omegak;

  intstatus = gsl_integration_qag(&BB,aval,1.0,epsabs,epsrel,limitI,GSL_INTEG_GAUSS31,workspace,&i1,&i1err);
  if(fabs(omegak) > 1.0e-5) {
    if(omegak > 0) {
      return (DH/hval/sqrt(omegak)*sinh(sqrt(omegak)*i1));
      }
    else { //omegak < 0
      return (DH/hval/sqrt(fabs(omegak))*sin(sqrt(fabs(omegak))*i1));
      }
    } //end |omegak| > 0
  else { //omegak == 0
    return (i1*DH/hval);
    }
  //shouldn't be here.
  exit(1);
  }

//c./H(z) gives the distance between pairs separated by dz along the line of sight.
double HinvMpc(double aval, void *params) {
  cosmo_params *myp = (cosmo_params *) params;
  double hval = myp -> h;
  return(DH/hval/normH(aval,params));
  }

//as defined by Eisenstein et al 2005 Eqn 2
double DVMpc(double aval, void *params) {
  double DM = comoving_distMpc(aval, params);
  double z = 1./aval - 1.;
  return(pow(DM*DM*HinvMpc(aval,params)*z,1./3.));
  }

double comoving_disthinvMpc_lcdm(double z, cosmo_params cosmopfid) {
  double aval = 1./(1.+z);
  return(comoving_distMpc(aval,&cosmopfid)*cosmopfid.h);
  }

double comoving_distMpc_lcdm(double z, cosmo_params cosmopfid) {
  double aval = 1./(1.+z);
  return(comoving_distMpc(aval,&cosmopfid));
  }



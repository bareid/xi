#include "header.h"


//at the moment only set up for log-binned angular weights.
angwgt readangwgt(char *filename1pw) {

  angwgt awgt;

  FILE *ifpw;
  int headertot, linestot,i;
  if(!(ifpw = open_file_read(filename1pw))) {
    exit(1);
    //return NULL;
    }
  get_file_length(ifpw,&headertot,&linestot);
  assert(headertot==0);

  awgt.ntheta = linestot;
  awgt.thetacen = (real *) malloc(sizeof(real)*awgt.ntheta);
  awgt.thetamin = (real *) malloc(sizeof(real)*awgt.ntheta);
  awgt.thetamax = (real *) malloc(sizeof(real)*awgt.ntheta);
  awgt.angweight = (real *) malloc(sizeof(real)*awgt.ntheta);

  for(i=0;i<awgt.ntheta;i++) {
    #ifdef REALFLOAT
    fscanf(ifpw,"%e %e %e %e\n",&(awgt.thetacen[i]), &(awgt.thetamin[i]), &(awgt.thetamax[i]), &(awgt.angweight[i]));
    #else
    fscanf(ifpw,"%le %le %le %le\n",&(awgt.thetacen[i]), &(awgt.thetamin[i]), &(awgt.thetamax[i]), &(awgt.angweight[i]));
    #endif
    if(i==0) {
      awgt.dlog10theta = log10(awgt.thetamax[i]/awgt.thetamin[i]);
      }
    else {
      assert(fabs(log10(awgt.thetamax[i]/awgt.thetamin[i]) - awgt.dlog10theta) < 2.0e-6);
      assert(fabs(log10(awgt.thetacen[i]/awgt.thetacen[i-1]) - awgt.dlog10theta) < 2.0e-6);
      }
    }
  fclose(ifpw);
  awgt.log10thetamin = log10(awgt.thetamin[0]);
  return awgt;
  }

real fbangweight(real asep, angwgt awgt) {
  if(asep < awgt.thetamin[0]) {
    return(awgt.angweight[0]);
    }
  if(asep > awgt.thetamax[awgt.ntheta-1]) {
    //no angular weighting at large angular separations.
    return(1.);
    }
  int tbin = min(((int) floor((log10(asep) - awgt.log10thetamin)/awgt.dlog10theta)),awgt.ntheta-1);
  #ifdef SANITYCHECKS
  if(!(asep >= awgt.thetamin[tbin]*0.99999 && asep <= awgt.thetamax[tbin]*1.00001)) {
    printf("%e %e %e %d\n",asep,awgt.thetamin[tbin],awgt.thetamax[tbin],tbin);
    fflush(stdout);
    }
  assert(asep >= awgt.thetamin[tbin]*0.99999 && asep <= awgt.thetamax[tbin]*1.00001);
  #endif
  //rbin = min(((int) floor((log10rsep-LOG10RSEPMIN)/DELTALOG10R)),NRBINS-1);
  return (awgt.angweight[tbin]);
  }


void freeangwgt(angwgt awgt) {
  free(awgt.thetacen);
  free(awgt.thetamin);
  free(awgt.thetamax);
  free(awgt.angweight);
  }




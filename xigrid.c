#include "header.h"

extern const int SORTINDX;

//stuff for working on the grid.
//copying from ~/boss/correlationfxnMASTERv6/countpairs.c and ~/boss/zdistvXlogbinsompcleverLSsmallscalev2bughunt/xiLS.c

int i2n(const int ipos[], int Ncell) {
  return(Ncell*Ncell*ipos[0]+Ncell*ipos[1]+ipos[2]);
  }

void n2i(int cellnum, int Ncell, int *ipos) {
  (ipos)[2] = cellnum%(Ncell);
  (ipos)[1] = ((cellnum - (ipos)[2])%(Ncell*Ncell))/(Ncell);
  (ipos)[0] = cellnum/(Ncell*Ncell);
  }


int chksortindxmatch(int icell, int jcell, int Ncell) {
  int iposi[3];
  int iposj[3];
  n2i(icell,Ncell,iposi);
  n2i(jcell,Ncell,iposj);
  int iii;
  for(iii=0;iii<=2;iii++) {
    }
  if(iposi[SORTINDX] == iposj[SORTINDX]) {
    return 1;
    }
  else {
    return 0;
    }
  }

int wrapgridperiodic(int ii, int Ncell) {
  if(ii < 0) return (ii+Ncell);
  if(ii >= Ncell) return (ii-Ncell);
  return ii;
  }

void nbrsprecompute(int *closenbrs, size_t nnbrsmax, int NFAC) {
  int ix, iy, iz;
  int dsqrx, dsqry, dsqrz;
  int cnt = 0;
  int nskip = 0;
  for(ix=-NFAC;ix<=NFAC;ix++) {
    dsqrx = min(min(ix*ix,(ix-1)*(ix-1)),(ix+1)*(ix+1));
    for(iy=-NFAC;iy<=NFAC;iy++) {
      dsqry = min(min(iy*iy,(iy-1)*(iy-1)),(iy+1)*(iy+1));
      for(iz=-NFAC;iz<=NFAC;iz++) {
        dsqrz = min(min(iz*iz,(iz-1)*(iz-1)),(iz+1)*(iz+1));
        if(dsqrx + dsqry + dsqrz < NFAC*NFAC) {
          closenbrs[cnt] = 1;
          }
        else {
          closenbrs[cnt] = 0;
          nskip += 1;
          }
        cnt += 1;
        }
      }
    }
  #ifdef ASSERT_ON
  assert(cnt == nnbrsmax);
  #endif
  #ifdef REALLYVERBOSE
  printf("skipped %d of %d\n",nskip,cnt);
  #endif
  }

//put periodicopt == 1 for periodic boundary conditions (i.e., in a box)
void getnbrs(int cellnum, int NFAC, int Ncell, int *nbrs, int *nnbrs, int *closenbrs, int periodicopt) {
  int iposcen[3];
  int ipos[3];
  //tmp!!
  int iposchk[3];
  n2i(cellnum,Ncell,iposcen);
  int ix, iy, iz;
  int nbrindx = 0;
  int closecnt = 0;
  for(ix=iposcen[0]-NFAC;ix<=iposcen[0]+NFAC;ix++) {
    for(iy=iposcen[1]-NFAC;iy<=iposcen[1]+NFAC;iy++) {
      for(iz=iposcen[2]-NFAC;iz<=iposcen[2]+NFAC;iz++) {
        if(closenbrs[closecnt] == 0) {
          closecnt += 1;
          continue;
          }
        else {
          closecnt += 1;
          }
        if(periodicopt == 0) {
          if(ix < 0 || ix >= Ncell || iy < 0 || iy >= Ncell || iz < 0 || iz >= Ncell) {
            continue;
            }
          ipos[0] = ix;
          ipos[1] = iy;
          ipos[2] = iz;
          }
        else { //periodicopt = 1
          ipos[0] = wrapgridperiodic(ix,Ncell);
          ipos[1] = wrapgridperiodic(iy,Ncell);
          ipos[2] = wrapgridperiodic(iz,Ncell);
          }
        nbrs[nbrindx] = i2n(ipos,Ncell);
        #ifdef SANITYCHECKS
        n2i(nbrs[nbrindx],Ncell,iposchk);
        if(iposchk[0] != ix && iposchk[0] != ix-Ncell && iposchk[0] != ix+Ncell) {
          printf("x err!\n");
          exit(1);
          }
        if(iposchk[1] != iy && iposchk[1] != iy-Ncell && iposchk[1] != iy+Ncell) {
          printf("y err!\n");
          exit(1);
          }
        if(iposchk[2] != iz && iposchk[2] != iz-Ncell && iposchk[2] != iz+Ncell) {
          printf("z err!\n");
          exit(1);
          }
        #endif
        nbrindx += 1;
        }
      }
    }
  (*nnbrs) = nbrindx;
  #if defined(SANITYCHECKS) && defined(ASSERT_ON)
  assert(closecnt == (2*NFAC+1)*(2*NFAC+1)*(2*NFAC+1));
  #endif
  }

//end stuff for working on the grid.




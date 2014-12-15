#include "header.h"

//need to fill in cell field here rather than the reads.  let's make everything self contained in this function!

//read in APscale from countpairs function.
//get AP stuff correct!  

//make this a struct to pass around.
/*
real RSEPMAX,RSEPMAXAP,rminsqr,rmaxsqr;
real APscale[3];
int Ncell;
int NFAC;
real Lbox;
real originpos[3];
*/

const int SORTINDX = 2;  //this shows up in the compare_sort, I don't know how to code it if not a global variable.


//set up bounds for rmin/rmax.
int setpairrminmax(xibindat b, int radeczorsim, real APperp, real APpar, real mindist, cntparams *cp) {

  double maxx, maxy;
  double maxxAP, maxyAP;
  double RSEPMIN;

  real APmax = max(APperp,APpar);

  if(b.bintype == 3) {
    cp->RSEPMAXAP = getrsep3dmaxforang(b,cp->originpos,-1.);
//we don't need RSEPMIN,RSEPMAX,rminsqr,rmaxsqr
    }

  if(b.bintype == 4) {
    cp->RSEPMAXAP = getrsep3dmaxforang(b,cp->originpos,mindist);
    if(b.logxopt == 1) {
      maxx = pow(10.,b.minx + ((double) b.nx)*b.dx);
      cp->rminsqr = pow(10.,b.minx)*pow(10.,b.minx);
      }
    else {
      maxx = (b.minx + ((double) b.nx)*b.dx);
      cp->rminsqr = (b.minx*b.minx);
      }
    cp->RSEPMAX = maxx;
    }

  if(b.bintype == 0 || b.bintype == 2) {  //xi_ell
    assert(b.logyopt == 0);
    if(b.logxopt == 1) {
      cp->RSEPMAX = pow(10.,b.minx + ((double) b.nx)*b.dx);
      RSEPMIN = pow(10.,b.minx);
      cp->rminsqr = RSEPMIN*RSEPMIN;
      }
    else {  //linear binning.
      cp->RSEPMAX = (b.minx + ((double) b.nx)*b.dx);
      RSEPMIN = (b.minx); 
      cp->rminsqr = RSEPMIN*RSEPMIN;
      }
    if(radeczorsim == 0) {
      cp->RSEPMAXAP = cp->RSEPMAX;
      }
    else { //sims with AP effect.
      cp->RSEPMAXAP = cp->RSEPMAX*APmax;
      }
    }  //xi_ell

  if(b.bintype == 1) { //xi_grid
    if(b.logxopt == 1) {
      maxx = pow(10.,b.minx + ((double) b.nx)*b.dx);
      cp->rminsqr = pow(10.,b.minx)*pow(10.,b.minx);
      }
    else {
      maxx = (b.minx + ((double) b.nx)*b.dx);
      cp->rminsqr = (b.minx*b.minx);
      }
    maxxAP = maxx*APperp;

    if(b.logyopt == 1) {
      maxy = pow(10.,b.miny + ((double) b.ny)*b.dy);
      cp->rminsqr = min(cp->rminsqr,pow(10.,b.miny)*pow(10.,b.miny));
      }
    else {
      maxy = (b.miny + ((double) b.ny)*b.dy);
      cp->rminsqr = min(cp->rminsqr,b.miny*b.miny);
      }
    maxyAP = maxy*APpar;

    cp->RSEPMAX = sqrt(maxx*maxx+maxy*maxy);
    if(radeczorsim == 0) {
      cp->RSEPMAXAP = cp->RSEPMAX;
      }
    else { //sims with AP effect.
      cp->RSEPMAXAP = sqrt(maxxAP*maxxAP+maxyAP*maxyAP);
      }
    }  //xi_grid

  cp->rmaxsqr = cp->RSEPMAX*cp->RSEPMAX;

  #ifdef REALLYVERBOSE
  printf("bintype: %d\n",b.bintype);
  printf("%e %e %e %e\n",cp->RSEPMAXAP,cp->rminsqr,cp->rmaxsqr,cp->RSEPMAX);
  #endif

  return 0;
  }


int compare_cntparticle_sort(const void *a, const void *b)  {
  const cntparticle *dat1 = (const cntparticle *) a;
  const cntparticle *dat2 = (const cntparticle *) b;
  if(dat1->cell != dat2->cell) {
    return (dat1->cell > dat2->cell) - (dat1->cell < dat2->cell);
    }
  else {
    return (dat1->pos[SORTINDX] > dat2->pos[SORTINDX]) - (dat1->pos[SORTINDX] < dat2->pos[SORTINDX]);
    }  
  }

/*
sepopt = 0: regular 3d correlation function on radecz.
sepopt = 1: angular 2d correlation function on radecz.
sepopt = 2: regular periodic box; does velocities or regular pair counting according to bintype.
sepopt = 3: hogg option.
*/

//python needs to know particle, xibindat.
int countpairs(cntparticle *plist,int np, int autoorcross, int sepopt, xibindat b, angwgt awgt, cntparams *cp, long double *Npairsfinal) {

  int angwgtopt = 0;
  if(awgt.ntheta != -1) {
    angwgtopt = 1;
    assert(b.bintype == 3);
    }

  if(b.logxopt == 1) {
    if(!(b.minx + b.dx*b.nx < 3.0)) {
      printf("b.minx needs to be log10!!!\n");
      exit(1);
      }
    assert(b.minx + b.dx*b.nx < 3.0);
    } //b.logxopt
  if(b.logyopt == 1) {
    if(!(b.miny + b.dy*b.ny < 3.0)) {
      printf("b.miny needs to be log10!!!\n");
      exit(1);
      }
    assert(b.miny + b.dy*b.ny < 3.0);
    } //b.logyopt

  int i,j,k;

  int bin2dtot=b.nx*b.ny;
  int nthreads;
  #ifdef NOOMP
  nthreads = 1;
  #else
  nthreads = omp_get_max_threads();
  #endif
  #ifdef REALLYVERBOSE
  printf("nthreads = %d\n",nthreads);
  #endif
  long double *Npairs;

  //round down to nearest integer to put one particle per box on average.
  cp->NFAC = ((int) floor(pow(((float) np),1./3.)/cp->Lbox*cp->RSEPMAXAP));

  cp->NFAC = min(cp->NFAC,16);
  cp->NFAC = max(cp->NFAC,1);
  assert(cp->NFAC >= 1 && cp->NFAC <= 16);
  cp->Ncell = ((int) floor(cp->Lbox/cp->RSEPMAXAP))*cp->NFAC;
  //let's try to set cp->NFAC automatically given RSEPMAX and np.
  if(cp->Ncell > pow(0.5*INT_MAX,1./3.)) {
    cp->NFAC = ((int) floor(pow(0.5*INT_MAX,1./3.)/floor(cp->Lbox/cp->RSEPMAXAP)));
    cp->Ncell = ((int) floor(cp->Lbox/cp->RSEPMAXAP))*cp->NFAC;
    }
  assert(cp->Ncell < pow(INT_MAX,1./3.));
  #ifdef REALLYVERBOSE
  printf("BR: remember to experiment with the NFAC determination!!\n");
  printf("here is NFAC and particles per cell: %d %d %f %e\n",cp->NFAC,cp->Ncell,cp->Lbox,((float) np)/(cp->Ncell*cp->Ncell*cp->Ncell));
  printf("setting up grid with Ncell = %d\n",cp->Ncell);
  #endif  

  //set up the grid.
  int ipos[3];
  int cellnum;

  size_t nnbrsmax = ((size_t) ((2*cp->NFAC+1)*(2*cp->NFAC+1)*(2*cp->NFAC+1)));
  int **nbrs, *closenbrs;
  // doesn't work!  need to malloc some space for each thread and fill accordingly.
  //nbrs = (int *) malloc(sizeof(int)*nnbrsmax);
  nbrs = malloc2dint(nthreads,nnbrsmax);

  //precompute which cells are kept and which aren't.
  closenbrs = (int *) malloc(sizeof(int)*nnbrsmax);

  //stuff to keep track of cell-> list map.
  int cellmax = cp->Ncell*cp->Ncell*cp->Ncell;
  int *cellstart, *cellend, *occupiedcells;
  int noccupied;
  cellstart = (int *) malloc(sizeof(int)*cellmax);
  cellend = (int *) malloc(sizeof(int)*cellmax);
  occupiedcells = (int *) malloc(sizeof(int)*cellmax);

/*
  int nbrs[((2*NFAC+1)*(2*NFAC+1)*(2*NFAC+1))];
  int closenbrs[((2*NFAC+1)*(2*NFAC+1)*(2*NFAC+1))];
*/

  nbrsprecompute(closenbrs,nnbrsmax,cp->NFAC);

  real vtmp1[3];
  real vtmp2[3];

  int nnbrs;  //in NONPERIODIC case, neighbors that require periodic wrapping are invalid.
  #ifdef SANITYCHECKS
  for(cellnum=0;cellnum<cp->Ncell*cp->Ncell*cp->Ncell;cellnum++) {
    n2i(cellnum,cp->Ncell,ipos);
    if(i2n(ipos,cp->Ncell) != cellnum) {
      printf("fail\n");
      exit(1);
      }
    }
  #endif
 
  //assign cell number to every object in the list.
  for(i=0;i<np;i++) {
    for(j=0;j<=2;j++) {
      ipos[j] = ((int) floor(plist[i].pos[j]*cp->Ncell/cp->Lbox));
      #ifdef SANITYCHECKS
      if(!(ipos[j] >= 0 && ipos[j] < cp->Ncell)) {
        printf("%d %e %e %d %d\n",i,plist[i].pos[j],plist[i].pos[j]/cp->Lbox,cp->Ncell,ipos[j]);
        }
      assert(ipos[j] >= 0 && ipos[j] < cp->Ncell);
      #endif
      }
    plist[i].cell = i2n(ipos,cp->Ncell);
    }
  qsort(plist,np,sizeof(cntparticle),compare_cntparticle_sort); 


  //this indicates no objects in the cells.
  for(i=0;i<cellmax;i++) {
    cellstart[i] = -1;
    cellend[i] = -1;
    occupiedcells[i] = -1;
    }

  //check what the sorted list looks like.
  assert(plist[0].cell >= 0);
  int currcell = plist[0].cell;
  cellstart[currcell] = 0;
  assert(plist[np-1].cell < cellmax);

  noccupied = 0;
  for(i=0;i<np;i++) {
    if(plist[i].cell != currcell) {
      cellend[currcell] = i-1;
      occupiedcells[noccupied] = currcell;
      noccupied += 1;
      currcell = plist[i].cell;
      cellstart[currcell] = i;
      }
    }
  printf("yoyo3\n");
  #ifdef REALLYVERBOSE
  printf("%d of %d cells occupied\n",noccupied,cellmax);
  #endif
  cellend[currcell] = np-1;
  if(occupiedcells[noccupied-1] != currcell) {
    occupiedcells[noccupied] = currcell;
    noccupied += 1;
    }

  //zero the Npairs counter in every bin.
  Npairs = (long double *) malloc(sizeof(long double)*(nthreads*bin2dtot));  
  for(i=0;i<b.nx;i++) {
    for(k=0;k<b.ny;k++) {
      for(j=0;j<nthreads;j++) {
        Npairs[bin2dtot*j+i*b.ny+k] = 0.0;  
        }
      }
    }

 
  //now loop over all neighboring cells and count pairs!

  //private variables inside pragma.
  int oo,xa,ya,nn,icell,jcell,jstart,jend,sgn,inbin,bin2d;
  float offset,tmpdiff;
  long double pairwgt;
  long double vr, v2perp, v2par;
  vr = 0.;
  v2perp = 0.;
  v2par = 0.;
  real angsep = 1.0; //need every nthread to have this.  set up option for angular upweighting?
  int mythread;

  //let's write a test loop.
/*
  #pragma omp parallel for default(none) \
    private(mythread,oo,icell,nnbrs,nn) shared(noccupied,occupiedcells,cp,closenbrs,nbrs,b,nnbrsmax)
  for(oo=0;oo<noccupied;oo++) {
    mythread = omp_get_thread_num();
    icell = occupiedcells[oo];
    getnbrs(icell,cp->NFAC,cp->Ncell,nbrs[omp_get_thread_num()],&nnbrs,closenbrs,b.periodicopt);
    if(oo%10 == 0) {
      printf("neighbors of %d\n",oo);
      for(nn=0;nn<nnbrsmax;nn++) {
        printf("%d %d %d %d\n",mythread,omp_get_thread_num(),nn,nbrs[omp_get_thread_num()][nn]);
        }
      }
    }
  exit(1);
*/

  #ifndef NOOMP
  #pragma omp parallel for default(none) \
              private(oo,xa,ya,nn,icell,jcell,nnbrs,jstart,jend,sgn,tmpdiff,offset,inbin,bin2d,pairwgt,vr,v2perp,v2par, angsep, mythread) shared(noccupied,occupiedcells,cellmax,closenbrs,nbrs,cellstart,cellend,b,cp,plist,autoorcross,sepopt,vtmp1,vtmp2,angwgtopt,awgt,bin2dtot,Npairs)
  #endif
  for(oo=0;oo<noccupied;oo++) {
    #ifdef REALLYVERBOSE
    if(oo%10000 == 0) {
      printf("working on %d of %d\n",oo,noccupied);
      }
    #endif
    icell = occupiedcells[oo];
    #if defined(SANITYCHECKS) && defined(ASSERT_ON)
    assert(icell >= 0 && icell < cellmax);
    #endif
    mythread = omp_get_thread_num();
    getnbrs(icell,cp->NFAC,cp->Ncell,nbrs[mythread],&nnbrs,closenbrs,b.periodicopt);
    for(nn=0;nn<nnbrs;nn++) {
//tmp!!
      assert(mythread == omp_get_thread_num());
      jcell = nbrs[mythread][nn];
      if(jcell < icell) {
        continue;
        }
      if(cellstart[jcell] == -1) {
        continue;
        }
      if(icell == jcell) {
        sgn = -2;
        }
      else { //different cells, so count all pairs, worry about Lbox.  as long as you're not going to sizes of the same order as Lbox/2, just define one offset of Lbox.
        //getnbrs doesn't wrap around if periodicopt == 0.
        tmpdiff = (plist[cellstart[jcell]].pos[SORTINDX] - plist[cellstart[icell]].pos[SORTINDX]);
        offset = 0.;
        sgn = 0;
        if(fabs(tmpdiff + cp->Lbox) < fabs(tmpdiff)) {
          offset = cp->Lbox;
          sgn = 1; //means jcell is to right of icell
          }
        if(fabs(tmpdiff - cp->Lbox) < fabs(tmpdiff)) {
          offset = -cp->Lbox;
          sgn = -1; //jcell to left of icell (after periodic wrapping)
          }
        if(fabs(offset) < 0.001*cp->Lbox) { //it's still 0.
          if(tmpdiff > 0.) {
            sgn = 1;
            }
          else {
            sgn = -1;
            }
          }
        //not offset in SORTINDX direction.
        if(chksortindxmatch(icell,jcell,cp->Ncell) == 1) {
          sgn = 0;
          }
        } //else icell == jcell
      //BR IMPROVEMENT -- I THINK THIS COULD BE SPED UP WITH A BINARY SEARCH INSTEAD OF A COMPARISON AGAINST EVERY OBJECT!!  Come back later.
      for(xa=cellstart[icell];xa<=cellend[icell];xa++) {
        switch(sgn) {
          case(-2):
        jend = cellend[jcell];
        jstart = xa+1;
        break;
          case(-1):
        jstart = cellend[jcell];
        jend = cellend[jcell];
        //check sign on x.  is diff pos or neg? also, need to add Lbox or not?
        while((jstart > cellstart[jcell]) && ((-plist[jstart].pos[SORTINDX] - offset + plist[xa].pos[SORTINDX]) < cp->RSEPMAXAP)) {
          jstart--;
          }
        while((jstart < cellend[jcell]) && ((-plist[jstart].pos[SORTINDX] - offset + plist[xa].pos[SORTINDX]) > cp->RSEPMAXAP)) {
          jstart++;
          }
        jstart--;
        if(jstart < cellstart[jcell]) {
          jstart = cellstart[jcell];
          }
        break;
          case(0):
        jstart = cellstart[jcell];
        jend = cellend[jcell];
        while((jstart < cellend[jcell]) && ((-plist[jstart].pos[SORTINDX] + plist[xa].pos[SORTINDX]) > cp->RSEPMAXAP)) {
          jstart++;
          }
        jstart--;
  
  //should we do the same while loop on jend??
        if(jend > cellend[jcell]) {
          jend = cellend[jcell];
          }
        if(jstart < cellstart[jcell]) {
          jstart = cellstart[jcell];
          }
        break;
          case(1):
        jstart = cellstart[jcell];
        jend = cellstart[jcell];
        //check sign on x.  is diff pos or neg? also, need to add Lbox or not?
        while((jend < cellend[jcell]) && ((plist[jend].pos[SORTINDX] + offset - plist[xa].pos[SORTINDX]) < cp->RSEPMAXAP)) {
          jend++;
          }
        while((jend > cellstart[jcell]) && ((plist[jend].pos[SORTINDX] + offset - plist[xa].pos[SORTINDX]) > cp->RSEPMAXAP)) {
          jend--;
          }
        jend++;
  
        //new!!!!!!
        if(jend > cellend[jcell]) {
          jend = cellend[jcell];
          }
        break;
          default:
        exit(1);
        }  //end switch(sgn)
        for(ya=jstart;ya<=jend;ya++) {
          if((autoorcross == 1) && (plist[xa].DorR == plist[ya].DorR)) {
            continue;
            }
          switch(sepopt) {
            case(0): //sepopt = 0: regular 3d correlation function on radecz.
          inbin = addpairsky3d(plist[xa].pos, plist[ya].pos, b, cp, &bin2d, &angsep);
          //old sanity check revive if I am bug-hunting again.
/*
              #ifdef SANITYCHECKS
              inbinchk = calcsep(gallist[xa], gallist[ya], &rbinchk, &mymu2chk);
              assert(inbinchk == inbin);
              assert(rbin == rbinchk);
              //printf("%d %Le %Le %Le\n",rbin,mymu2,mymu2chk,mymu2-mymu2chk);
              assert(fabs(mymu2chk - mymu2) < 1.0e-9);
              #endif

*/
          break; //end sepopt = 0
            case(1):  //sepopt = 1: angular 2d correlation function on radecz.
          inbin = addpairsky2d(plist[xa].pos, plist[ya].pos, cp->originpos, b, &bin2d);
          break; //end sepopt = 1
            case(2): //sepopt = 2: regular periodic box; optional velocity stuff.
          #ifdef VOPT
          inbin = addpairperiodic(plist[xa].pos, plist[ya].pos, plist[xa].vel, plist[ya].vel, b, cp, &bin2d, &vr, &v2perp, &v2par);
          #else
          inbin = addpairperiodic(plist[xa].pos, plist[ya].pos, vtmp1, vtmp2, b, cp, &bin2d, &vr, &v2perp, &v2par);
          #endif
          break; //end sepopt = 2
            case(3):
          if(plist[xa].DorR == 0) {
            inbin = addpairsky2dhogg(plist[xa].pos, plist[ya].pos, cp->originpos, b, &bin2d, plist[xa].chi);
            }
          else {
            inbin = addpairsky2dhogg(plist[xa].pos, plist[ya].pos, cp->originpos, b, &bin2d, plist[ya].chi);
            }
          break;
            default:
          exit(1);
          } //end switch on sepopt (type of pair counter).
          if(inbin == 1) {
            if(angwgtopt == 1) {
              pairwgt = ((long double) (plist[xa].wgt*plist[ya].wgt*fbangweight(angsep,awgt)));
              }
            else {
              pairwgt = ((long double) (plist[xa].wgt*plist[ya].wgt));
              }
            #ifdef NOOMP
            Npairs[bin2d] += pairwgt;
            #else
            Npairs[bin2dtot*omp_get_thread_num()+bin2d] += pairwgt;
            #endif

            if(b.bintype == 2)  {
              #ifdef NOOMP
              Npairs[bin2d+1] += pairwgt*vr;
              Npairs[bin2d+2] += pairwgt*v2par;
              Npairs[bin2d+3] += pairwgt*v2perp;
              #else
              Npairs[bin2dtot*omp_get_thread_num()+bin2d+1] += pairwgt*vr;
              Npairs[bin2dtot*omp_get_thread_num()+bin2d+2] += pairwgt*v2par;
              Npairs[bin2dtot*omp_get_thread_num()+bin2d+3] += pairwgt*v2perp;
              #endif
              }
            } //end if inbin == 1
          } //end ya loop.
        } //end xa loop.
      } //end nn loop over icell neighbors.
    } //end icell loop

  //now sum pair counts.  Put all counts in the first NRBINS.
  long double bigsum = 0.;
  for(i=0;i<b.nx;i++) {
    for(k=0;k<b.ny;k++) {
      for(j=1;j<nthreads;j++) {
        Npairs[i*b.ny+k] += Npairs[bin2dtot*j+i*b.ny+k];
        }
      Npairsfinal[i*b.ny+k] = Npairs[i*b.ny+k];
      bigsum += Npairsfinal[i*b.ny+k];
//      printf("%d %d %Le\n",i,k,Npairsfinal[i*b.ny+k]);
      }
    }
  malloc2dfreeint(nbrs);
  free(closenbrs);
  free(cellstart);
  free(cellend);
  free(occupiedcells);
  return 0;
  } //end countpairs.

/*
sepopt = 0: regular 3d correlation function on radecz.
sepopt = 1: angular 2d correlation function on radecz.
sepopt = 2: regular periodic box.
*/

//int countpairssim(particle *plist, int np, int autoorcross, real Lbox, xibindat b, real APperp, real APpar, long double *Npairsfinal) {
int countpairssim(particle *p1, int np1, particle *p2, int np2, real Lbox, xibindat b, real APperp, real APpar, long double *Npairsfinal) {

  assert(b.zspaceaxis >= 0 && b.zspaceaxis <= 2); //if originally -1, xi.c sets it to 0 and realorzspace = 0.
  angwgt awgt;
  awgt.ntheta = -1; //this thing should never get used for sims!

  cntparams cp;

  int i,j;
  int sepopt;
  int autoorcross = 0;
  int ntot;
  cntparticle *ctot;

  if(np2 > 0) {
    assert(p2 != NULL);
    autoorcross = 1;
    ntot = np1 + np2;
    }
  else {
    ntot = np1;
    }

  ctot = (cntparticle *) malloc(sizeof(cntparticle)*(ntot));


  //set global variables for Lbox, originpos, and APscale.
  cp.Lbox = Lbox;
  for(i=0;i<=2;i++) {
    cp.originpos[i] = 0.;
    }

  for(i=0;i<=2;i++) {
    cp.APscale[i] = APperp;
    }
  if(b.realorzspace == 0) {
    assert(APperp == APpar);
    }
  else { //calculation in redshift space.
    cp.APscale[b.zspaceaxis] = APpar;
    }
  setpairrminmax(b, 1, APperp, APpar, -1, &cp);
  sepopt = 2; // will do velocity or regular pair counts.

  //change type to cntparticle, put in one list and sort.
  for(i=0;i<np1;i++) {
    for(j=0;j<=2;j++) {
      ctot[i].pos[j] = p1[i].pos[j]; //no originpos
      #ifdef VOPT
      ctot[i].vel[j] = p1[i].vel[j];
      #endif
      ctot[i].wgt = p1[i].weight;
      ctot[i].DorR = 0;
      }
    }
  for(i=0;i<np2;i++) {
    for(j=0;j<=2;j++) {
      ctot[i+np1].pos[j] = p2[i].pos[j];
      #ifdef VOPT
      ctot[i+np1].vel[j] = p2[i].vel[j];
      #endif
      ctot[i+np1].wgt = p2[i].weight;
      ctot[i+np1].DorR = 1;
      }
    }

  countpairs(ctot,ntot,autoorcross,sepopt,b,awgt,&cp,Npairsfinal);
  free(ctot);
  return 0;
  } //end countpairssims

//functions that can be called from python instead of c.
//doesn't work yet.
/*
int countpairssimP(real *p1, real *v1, real *w1, int np1, real *p2, real *v2, real *w2, int np2, real Lbox, xibindat b, real APperp, real APpar, long double *Npairsfinal) {

  assert(b.zspaceaxis >= 0 && b.zspaceaxis <= 2); //if originally -1, xi.c sets it to 0 and realorzspace = 0.
  angwgt awgt;
  awgt.ntheta = -1; //this thing should never get used for sims!

  cntparams cp;

  int i,j;
  int sepopt;
  int autoorcross = 0;

  if(np2 > 0) {
    assert(p2 != NULL);
    autoorcross = 1;
    }

  cntparticle *ctot;
  int ntot = np1 + np2;
  ctot = (cntparticle *) malloc(sizeof(cntparticle)*(ntot));


  //set global variables for Lbox, originpos, and APscale.
  cp.Lbox = Lbox;
  for(i=0;i<=2;i++) {
    cp.originpos[i] = 0.;
    }

  for(i=0;i<=2;i++) {
    cp.APscale[i] = APperp;
    }
  if(b.realorzspace == 0) {
    assert(APperp == APpar);
    }
  else { //calculation in redshift space.
    cp.APscale[b.zspaceaxis] = APpar;
    }
  setpairrminmax(b, 1, APperp, APpar, &cp);
  sepopt = 2; // will do velocity or regular pair counts.

  //change type to cntparticle, put in one list and sort.
  for(i=0;i<np1;i++) {
    for(j=0;j<=2;j++) {
      ctot[i].pos[j] = p1[i].pos[j]; //no originpos
      #ifdef VOPT
      ctot[i].vel[j] = p1[i].vel[j];
      #endif
      ctot[i].wgt = p1[i].weight;
      ctot[i].DorR = 0;
      }
    }
  for(i=0;i<np2;i++) {
    for(j=0;j<=2;j++) {
      ctot[i+np1].pos[j] = p2[i].pos[j];
      #ifdef VOPT
      ctot[i+np1].vel[j] = p2[i].vel[j];
      #endif
      ctot[i+np1].wgt = p2[i].weight;
      ctot[i+np1].DorR = 1;
      }
    }

  countpairs(ctot,ntot,autoorcross,sepopt,b,awgt,&cp,Npairsfinal);
  free(ctot);
  return 0;
  } //end countpairssims
*/

//BR: generalize to read in ra,dec and do the spherical coordinate transformation here?  LATER.
int countpairsradecz(particle *p1, int np1, particle *p2, int np2, real mindist, real maxdist, xibindat b, angwgt awgt, long double *Npairsfinal) {
//  printf("yo beth, %e %e\n",mindist,maxdist);
//this assert goes off for funny locations of bootstrap regions, so let's skip this for now!
/*
  if(mindist >= maxdist) {
    assert(maxdist <= 1.0); //make sure it's an angopt Hogg situation.
    }
*/

  int angopt = 0;
  if(b.bintype == 3 || b.bintype == 4) {
    angopt = 1;
    assert(b.ny == 1);
    }


  int i,j;
  int autoorcross = 0;
  int sepopt;
  cntparams cp;

  int ntot;
  cntparticle *ctot;

  if(np2 > 0) {
    assert(p2 != NULL);
    autoorcross = 1;
    ntot = np1 + np2;
    }
  else {
    ntot = np1;
    }
  ctot = (cntparticle *) malloc(sizeof(cntparticle)*(ntot));

  //determine Lbox given maxdist from the origin in the catalog.
  if(angopt == 0) { //computing a 3d correlation fxn with sky coordinates.
    cp.Lbox = 2.*maxdist*1.01; //avoid floating point rounding errors.
    }
  else { //computing a 2d (angular) correlation fxn with sky coordinates.
    assert(b.bintype == 3 || b.bintype == 4);
    cp.Lbox = 2.;
    }
  for(i=0;i<=2;i++) {
    cp.originpos[i] = (cp.Lbox)*0.5;
    }

  for(i=0;i<np1;i++) {
    for(j=0;j<=2;j++) {
      ctot[i].pos[j] = p1[i].pos[j] + cp.originpos[j];
      #ifdef VOPT
      ctot[i].vel[j] = p1[i].vel[j];
      #endif
      ctot[i].wgt = p1[i].weight;
      ctot[i].chi = p1[i].chi;
      ctot[i].DorR = 0;
      }
    }
  for(i=0;i<np2;i++) {
    for(j=0;j<=2;j++) {
      ctot[i+np1].pos[j] = p2[i].pos[j] + cp.originpos[j];
      #ifdef VOPT
      ctot[i+np1].vel[j] = p2[i].vel[j];
      #endif
      ctot[i+np1].wgt = p2[i].weight;
      ctot[i+np1].DorR = 1;
      }
    }

  setpairrminmax(b, 0, 1.0, 1.0,mindist,&cp);

  if(angopt == 0) {
    sepopt = 0;
    }
  else {
    sepopt = 1;
    }
  if(b.bintype == 4) {
    sepopt = 3;
    }
  countpairs(ctot,ntot,autoorcross,sepopt,b,awgt,&cp,Npairsfinal);
  free(ctot);
  return 0;
  } //end countpairsradecz.




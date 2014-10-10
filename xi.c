#include "header.h"

typedef struct {
  int radeczorsim;
  int unitsMpc;
  double omfid;
  double hfid;
  char binfname[MAXLINELEN];
  char foutbase[MAXLINELEN];
  int Dftype;
  char Dfilename[MAXLINELEN];
  int Rftype;
  char Rfilename[MAXLINELEN];
  float zmin;
  float zmax;
  int DRopt;
  float ndownRR;
  float ndownDD;
  unsigned int unsignedintseed; 
  char fname1pw[MAXLINELEN];
  int angupweight;
  real Lbox;
  float abox;
  int zspaceaxis;
  float lg10Mmin;
  float lg10Mmax;
  float lg10Mmin2;
  float lg10Mmax2;
  int D2ftype;
  char D2filename[MAXLINELEN];
  float APperp;
  float APpar;
  } runparams;

const int nrunparams = 200;  //when finished coding, make this the actual final number of parameters in the parameter file.

void printrundat(runparams p) {
  printf("Common variables:\n");
  printf("radeczorsim: %d\n",p.radeczorsim);
  printf("unitsMpc: %d\n",p.unitsMpc);
  printf("omfid: %f\n",p.omfid);
  printf("hfid: %f\n",p.hfid);
  printf("binfname: %s\n",p.binfname);
  printf("foutbase: %s\n",p.foutbase);
  printf("Dftype: %d\n",p.Dftype);
  printf("Dfilename: %s\n",p.Dfilename);
  if(p.radeczorsim == 0) {
    printf("\n radecz specific variables:\n");
    printf("Rftype: %d\n",p.Rftype);
    printf("Rfilename: %s\n",p.Rfilename);
    printf("zmin: %f\n",p.zmin);
    printf("zmax: %f\n",p.zmax);
    printf("DRopt: %d\n",p.DRopt);
    printf("ndownRR: %f\n",p.ndownRR);
    printf("ndownDD: %f\n",p.ndownDD);
    printf("unsignedintseed: %u\n",p.unsignedintseed);
    }
  else {
    assert(p.radeczorsim == 1);
    printf("\nsim specific variables:\n");
    printf("Lbox: %f\n",p.Lbox);
    printf("abox: %f\n",p.abox);
    printf("zspaceaxis: %d\n",p.zspaceaxis);
    printf("lg10Mmin: %f\n",p.lg10Mmin);
    printf("lg10Mmax: %f\n",p.lg10Mmax);
    printf("lg10Mmin2: %f\n",p.lg10Mmin2);
    printf("lg10Mmax2: %f\n",p.lg10Mmax2);
    printf("D2ftype: %d\n",p.D2ftype);
    if(p.D2ftype != -1) {
      printf("D2filename: %s\n",p.D2filename);
      } 
    printf("APperp: %f\n",p.APperp);
    printf("APpar: %f\n",p.APpar);
    }
  } //end printrundat.

runparams parserunfile(char *infilename) {
  runparams myparams;
  FILE *ifpparams;
  ifpparams = open_file_read(infilename);
  char line[MAXLINELEN];
  char *r;

  int gotlist[nrunparams];
  int i;

  myparams.DRopt = -1; //for sims.

  for(i=0;i<nrunparams;i++) {
    gotlist[i] = 0;
    }

  while(fgets(line,MAXLINELEN,ifpparams))  {
    if(strstr(line,"#"))  {
      continue;
      }
    if(!feof(ifpparams))  {
      if(strstr(line,"radeczorsim")) {
        r = strchr(line,'=');
        assert(r);
        myparams.radeczorsim = atoi(r+1);
        assert(myparams.radeczorsim == 0 || myparams.radeczorsim == 1);
        gotlist[0] = 1;
        continue;
        }
      if(strstr(line,"unitsMpc")) {
        r = strchr(line,'=');
        assert(r);
        myparams.unitsMpc = atoi(r+1);
        assert(myparams.unitsMpc == 0 || myparams.unitsMpc == 1);
        gotlist[1] = 1;
        continue;
        }
      if(strstr(line,"omfid")) {
        r = strchr(line,'=');
        assert(r);
        myparams.omfid = atof(r+1);
        assert(myparams.omfid <= 1. && myparams.omfid >= 0.);
        gotlist[2] = 1;
        continue;
        }
      if(strstr(line,"hfid")) {
        r = strchr(line,'=');
        assert(r);
        myparams.hfid = atof(r+1);
        assert(myparams.hfid > 0. && myparams.hfid < 3.);  //make sure its in the right units, 100 h km/s/Mpc.
        gotlist[3] = 1;
        continue;
        }
      if(strstr(line,"binfname")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.binfname);
        gotlist[4] = 1;
        continue;
        }
      if(strstr(line,"foutbase")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.foutbase);
        gotlist[5] = 1;
        continue;
        }
      if(strstr(line,"Dftype")) {
        r = strchr(line,'=');
        assert(r);
        myparams.Dftype = atoi(r+1);
        gotlist[6] = 1;
        continue;
        }
      if(strstr(line,"Dfilename")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.Dfilename);
        gotlist[7] = 1;
        continue;
        }
      //end parameters common to sims and data.
      if(strstr(line,"Rftype")) {
        r = strchr(line,'=');
        assert(r);
        myparams.Rftype = atoi(r+1);
        gotlist[8] = 1;
        continue;
        }
      if(strstr(line,"Rfilename")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.Rfilename);
        gotlist[9] = 1;
        continue;
        }
      if(strstr(line,"zmin")) {
        r = strchr(line,'=');
        assert(r);
        myparams.zmin = atof(r+1);
        gotlist[10] = 1;
        continue;
        }
      if(strstr(line,"zmax")) {
        r = strchr(line,'=');
        assert(r);
        myparams.zmax = atof(r+1);
        gotlist[11] = 1;
        continue;
        }
      if(strstr(line,"DRopt")) {
        r = strchr(line,'=');
        assert(r);
        myparams.DRopt = atoi(r+1);
        gotlist[12] = 1;
        continue;
        }
      if(strstr(line,"ndownRR")) {  //optional for data.
        r = strchr(line,'=');
        assert(r);
        myparams.ndownRR = atof(r+1);
        gotlist[13] = 1;
        continue;
        }
      if(strstr(line,"ndownDD")) {  //optional for data.
        r = strchr(line,'=');
        assert(r);
        myparams.ndownDD = atof(r+1);
        gotlist[14] = 1;
        continue;
        }
      if(strstr(line,"rngenseed")) { //optional for data.
        r = strchr(line,'=');
        assert(r);
        myparams.unsignedintseed = ((unsigned int) (atoi(r+1)));
        gotlist[15] = 1;
        assert(myparams.unsignedintseed > 0);
        continue;
        }
      if(strstr(line,"1pwfname")) { //optional for data.
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.fname1pw);
        myparams.angupweight = 1;
        gotlist[16] = 1;
        continue;
        }
      if(strstr(line,"Lbox")) {
        r = strchr(line,'=');
        assert(r);
        myparams.Lbox = atof(r+1);
        gotlist[17] = 1;
        continue;
        }
      if(strstr(line,"abox")) {
        r = strchr(line,'=');
        assert(r);
        myparams.abox = atof(r+1);
        gotlist[18] = 1;
        continue;
        }
      if(strstr(line,"zspaceaxis")) {
        r = strchr(line,'=');
        assert(r);
        myparams.zspaceaxis = atoi(r+1);
        assert(myparams.zspaceaxis >= -1 && myparams.zspaceaxis <= 2);
        gotlist[19] = 1;
        continue;
        }
      if(strstr(line,"lg10Mmin")) {
        r = strchr(line,'=');
        assert(r);
        myparams.lg10Mmin = atof(r+1);
        gotlist[20] = 1;
        continue;
        }
      if(strstr(line,"lg10Mmax")) {
        r = strchr(line,'=');
        assert(r);
        myparams.lg10Mmax = atof(r+1);
        gotlist[21] = 1;
        continue;
        }
      if(strstr(line,"lg10Mmin2")) {
        r = strchr(line,'=');
        assert(r);
        myparams.lg10Mmin2 = atof(r+1);
        gotlist[22] = 1;
        continue;
        }
      if(strstr(line,"lg10Mmax2")) {
        r = strchr(line,'=');
        assert(r);
        myparams.lg10Mmax2 = atof(r+1);
        gotlist[23] = 1;
        continue;
        }
      if(strstr(line,"D2ftype")) {
        r = strchr(line,'=');
        assert(r);
        myparams.D2ftype = atoi(r+1);
        gotlist[24] = 1;
        continue;
        }
      if(strstr(line,"D2filename")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.D2filename);
        gotlist[25] = 1;
        continue;
        }
      if(strstr(line,"APperp")) {
        r = strchr(line,'=');
        assert(r);
        myparams.APperp = atof(r+1);
        gotlist[26] = 1;
        continue;
        }
      if(strstr(line,"APpar")) {
        r = strchr(line,'=');
        assert(r);
        myparams.APpar = atof(r+1);
        gotlist[27] = 1;
        continue;
        }
      } //end !feof
    else {
      fprintf(stderr,"non-key word line: %s\n",line);
      fprintf(stderr,"Quitting!\n");
      exit(1);
      }
    } //end while

  assert(gotlist[0] == 1); //got radeczorsim.
  switch(myparams.radeczorsim) {
    case(0):
  for(i=0;i<=12;i++) {
    if(gotlist[i] != 1) {
      fprintf(stderr,"missing elt %d from radecz runparams, aborting!\n",i);
      exit(1);
      }
    }
  //fill in optional input parameters with defaults.
  if(gotlist[13] == 0) {  //didn't get ndown.
    myparams.ndownRR = 1.0;
    }
  if(gotlist[14] == 0) {  //didn't get ndown.
    myparams.ndownDD = 1.0;
    }
  if(gotlist[15] == 0) {
    myparams.unsignedintseed = devrand();
    if(myparams.unsignedintseed == 1)  {
      fprintf(stderr,"initialize_rngen failed!  Exiting\n");
      exit(1);
      }
    }
  if(gotlist[16] == 0) {
    myparams.angupweight = 0;
    }
  break;  
    case(1):
  for(i=0;i<=7;i++) { 
    if(gotlist[i] != 1) {
      fprintf(stderr,"missing elt %d from sim runparams, aborting!\n",i);
      exit(1);
      }
    }
  for(i=17;i<=19;i++) { 
    if(gotlist[i] != 1) {
      fprintf(stderr,"missing elt %d from sim runparams, aborting!\n",i);
      exit(1);
      }
    }
  //specify default behavior for the optional simulation inputs.
  if(gotlist[20] == 0) { 
    myparams.lg10Mmin = -1.;
    }
  if(gotlist[21] == 0) {  
    myparams.lg10Mmax = 17.;
    }
  if(gotlist[22] == 0) {  
    myparams.lg10Mmin2 = myparams.lg10Mmin;
    }
  if(gotlist[23] == 0) {  
    myparams.lg10Mmax2 = myparams.lg10Mmax;
    }
  if(gotlist[24] == 0) {  
    myparams.D2ftype = -1;
    }
  if(gotlist[26] == 0) {  
    myparams.APperp = 1.0;
    }
  if(gotlist[27] == 0) {  
    myparams.APpar = 1.0;
    }
  //no need to specify default for D2filename.
  break;
    default:
  fprintf(stderr,"Illegal radeczorsim.  Quitting!\n");
  exit(1);
  }
  return myparams;
  } //end parserunfile.

/*
void set_Lbox_origin(runparams runp, cosmo_params cosmopfid, int angopt, real originpos[3], real *Lbox) {
  float zmaxD, zmaxR,myzmax;
  double maxchi;
  int i;
  if(runp.radeczorsim == 0 && angopt == 0) { //computing a 3d correlation fxn with sky coordinates.
//this is too slow, we'll just put in a check/assert in radecz.

//    zmaxD = getcatzmax(runp.Dfilename,runp.Dftype);
//    zmaxR = getcatzmax(runp.Rfilename,runp.Rftype);
//    myzmax = max(zmaxD,zmaxR);
//    if(runp.zmax < myzmax) {
//      myzmax = runp.zmax;
//      }
    myzmax = runp.zmax;
    if(runp.unitsMpc == 1) {
      maxchi = comoving_distMpc_lcdm((double) myzmax, cosmopfid);
      }
    else {
      maxchi = comoving_disthinvMpc_lcdm((double) myzmax, cosmopfid);
      }

    (*Lbox) = ((real) (2.*maxchi));
    for(i=0;i<=2;i++) {
      originpos[i] = (*Lbox)*0.5;
      }
    }
  if(runp.radeczorsim == 0 && angopt == 1) { //computing a 2d (angular) correlation fxn with sky coordinates.
    (*Lbox) = 2.;
    for(i=0;i<=2;i++) {
      originpos[i] = (*Lbox)*0.5;
      }
    }
  if(runp.radeczorsim == 1) {
    (*Lbox) = runp.Lbox; //read Lbox in from parameter file.
    for(i=0;i<=2;i++) {
      originpos[i] = 0.;
      }
    } 
  } //end set_Lbox_origin
*/

int main(int argc, char *argv[]) {
  if(argc != 2) {
    fprintf(stderr,"Usage ./xi pdata.params\n");
    fprintf(stderr,"Usage ./xi psims.params\n");
    exit(1);
    }

//variables needed by everyone.
  xibindat b;
  angwgt awgt;


//variables needed for sims.
  real vscale;
//  real APscale[3];

  runparams runp = parserunfile(argv[1]);
  printrundat(runp);
  docosmosetup(); //allocates integration memory.
  cosmo_params cosmopfid = (cosmo_params) {runp.omfid,1.-runp.omfid,0.,runp.hfid};

  b = readbinfile(runp.binfname); 
  int angopt = 0;
  if(b.bintype == 3 || b.bintype == 4) {
    angopt = 1;
    assert(b.ny == 1);
    }

  if(runp.angupweight == 1) {
    awgt = readangwgt(runp.fname1pw);
    }
  else {
    awgt.ntheta = -1; //signal to countpairs that we're not using angular upweighting!
    }

  //this function needs to know b to find out if it's angular correlation fxn or 
  //now we do this in countpairs!
  //set_Lbox_origin(runp,cosmopfid, angopt, originpos, &Lbox);

  //only needed if ndownRR != 1, but we'll initialize it anyways.
  initialize_rngen(runp.unsignedintseed);


  time_t time0, time1, time2;
  #ifdef REALLYVERBOSE
  printf("starting read\n");
  #endif
  (void) time(&time0);

  int DorR = 0;
  int autoorcross = 0;
  float ndownRRnone = -1;  //indicates do not do catalog downsampling.

  int n1, n2;
  long double n1wgt, n2wgt;
  long double n1wgtSI, n2wgtSI;
  particle *c1tmp, *c2tmp;

  int i;

  real maxdist1 = -1000.;
  real maxdist2 = -1000.;
  real maxdist; //this will be maxdist between the two catalogs.

  //this is needed for the Hogg method.  Minimum distance scale sets maximum angular scale you need to look for pairs at.
  real mindist1 = 200000.;
  real mindist2 = 200000.;
  real mindist; //this will be mindist between the two catalogs.

  n1wgt = 0.;
  n2wgt = 0.;
  n1wgtSI = 0.;
  n2wgtSI = 0.;

  //for 
  c2tmp = NULL;
  assert(c2tmp == NULL);
  n2 = -1;

  if(runp.radeczorsim == 0) {
    #ifdef REALLYVERBOSE
    printf("initializing RNGEN with %u\n",runp.unsignedintseed);
    #endif
    b.periodicopt = 0;
    b.realorzspace = 1;

    assert(runp.Dftype >= 1 && runp.Dftype <= 4);
    assert(runp.Rftype >= 1 && runp.Rftype <= 4); 

//BR next: put in read in data and/or random catalogs.
    assert((runp.DRopt >= 1 && runp.DRopt <= 3) || (runp.DRopt >= 11 && runp.DRopt <= 14));
    if(runp.DRopt == 1) {
      DorR = 0;
      autoorcross = 0;
      c1tmp = readcat(runp.Dfilename, runp.Dftype, runp.unitsMpc, angopt, cosmopfid, runp.zmin, runp.zmax, DorR, ndownRRnone, &n1, &n1wgt, &n1wgtSI, &mindist1, &maxdist1);
      }
    if(runp.DRopt == 3) {
      DorR = 1;
      autoorcross = 0;
      c1tmp = readcat(runp.Rfilename, runp.Rftype, runp.unitsMpc, angopt, cosmopfid, runp.zmin, runp.zmax, DorR, runp.ndownRR, &n1, &n1wgt, &n1wgtSI, &mindist1, &maxdist1);
      }
    if(runp.DRopt == 2) {  //need to read in both D and R catalogs and condense them.
      //I think we want to enforce no downsampling for DR counts.  It's not necessary in principle
      //but I think most likely it's a typo if I set ndownRR for DRopt == 2.
      autoorcross = 1;
      assert(fabs(runp.ndownRR - 1.0) < 2.0e-5);
      DorR = 0;
      c1tmp = readcat(runp.Dfilename, runp.Dftype, runp.unitsMpc, angopt, cosmopfid, runp.zmin, runp.zmax, DorR, ndownRRnone, &n1, &n1wgt, &n1wgtSI, &mindist1,&maxdist1);
      DorR = 1;
      c2tmp = readcat(runp.Rfilename, runp.Rftype, runp.unitsMpc, angopt, cosmopfid, runp.zmin, runp.zmax, DorR, runp.ndownRR, &n2, &n2wgt, &n2wgtSI, &mindist2,&maxdist2);
      }
    if(runp.DRopt >= 11 && runp.DRopt <= 14) {  //now "random" catalog is imaging, "data" catalog is spectroscopy.
      autoorcross = 1;
      if(runp.DRopt <= 12) {
        assert(fabs(runp.ndownDD - 1.0) < 2.0e-5);
        }
      if(runp.DRopt == 11 || runp.DRopt == 13) {
        assert(fabs(runp.ndownRR - 1.0) < 2.0e-5);
        }
      DorR = 0;
      c1tmp = readcat(runp.Dfilename, runp.Dftype, runp.unitsMpc, angopt, cosmopfid, runp.zmin, runp.zmax, DorR, runp.ndownDD, &n1, &n1wgt, &n1wgtSI, &mindist1,&maxdist1);
      DorR = 1;
      c2tmp = readcat(runp.Rfilename, runp.Rftype, runp.unitsMpc, angopt, cosmopfid, runp.zmin, runp.zmax, DorR, runp.ndownRR, &n2, &n2wgt, &n2wgtSI, &mindist2,&maxdist2);
      } //end DRopt == 11-14

    maxdist = max(maxdist1,maxdist2);
//    mindist = min(mindist1,mindist2);
    mindist = mindist1; //use the spectroscopic catalog to get mindist.
#ifdef REALLYVERBOSE
    printf("Using maxdist = %e\n",maxdist);
    printf("Using mindist = %e\n",mindist);
#endif
    } //end ra,dec,z preliminaries..

  else { //sim preliminaries.
    if(runp.zspaceaxis == -1) {
      runp.zspaceaxis = 0;
      b.realorzspace = 0;
      assert(runp.APperp == runp.APpar);
      }
    else {
      b.realorzspace = 1;
      }  
    b.periodicopt = 1;  //not fixing this took hours of bug hunting!! gaa!!
    b.zspaceaxis = runp.zspaceaxis;
    b.ellmaxdata = 4;
    assert(runp.Dftype >= 0 && runp.Dftype <= 5);
    if(runp.Dftype == 3) {
      if(runp.D2ftype != -1) {
        assert(runp.D2ftype == 3);
        }
      vscale = 100.*sqrt(runp.omfid/(runp.abox*runp.abox*runp.abox) + (1-runp.omfid))*runp.abox;
      vscale = 1./vscale;
      #ifdef REALLYVERBOSE
      printf("this is the vscale I calculated: %e\n",vscale);
      #endif
      }
    else { 
      vscale = 1.;
      }
    autoorcross = 1; //assume cross-correlation unless they're equal.
    assert(runp.Dftype >= 0 && runp.Dftype <= 5);

    if(runp.D2ftype != -1) {  //cross-correlation with dark matter particles.
      assert(runp.D2ftype >= 0 && runp.D2ftype <= 5);
      runp.lg10Mmin2 = -1.;
      runp.lg10Mmax2 = 17.;
      DorR = 0;
      c1tmp = readsimcat(runp.Dfilename,runp.Dftype,runp.Lbox,runp.zspaceaxis,vscale,runp.lg10Mmin,runp.lg10Mmax,DorR,&n1);
      DorR = 1;
      c2tmp = readsimcat(runp.D2filename,runp.D2ftype,runp.Lbox,runp.zspaceaxis,vscale,runp.lg10Mmin2,runp.lg10Mmax2,DorR,&n2);
      }
    else {
      if(fabs(runp.lg10Mmin-runp.lg10Mmin2) < 0.001 && fabs(runp.lg10Mmax-runp.lg10Mmax2) < 0.001) {
        autoorcross = 0;
        #ifdef REALLYVERBOSE
        printf("Mass bins are assumed equal, setting autoorcross = 0\n");
        #endif
        DorR = 0;
        c1tmp = readsimcat(runp.Dfilename,runp.Dftype,runp.Lbox,runp.zspaceaxis,vscale,runp.lg10Mmin,runp.lg10Mmax,DorR,&n1);
        }
      else { //don't allow overlap between mass bins.
        assert(runp.lg10Mmin < runp.lg10Mmax);
        assert(runp.lg10Mmin2 < runp.lg10Mmax2);
        assert((runp.lg10Mmax <= runp.lg10Mmin2) || (runp.lg10Mmax2 <= runp.lg10Mmin));
        DorR = 0;
        c1tmp = readsimcat(runp.Dfilename,runp.Dftype,runp.Lbox,runp.zspaceaxis,vscale,runp.lg10Mmin,runp.lg10Mmax,DorR,&n1);
        DorR = 1;
        c2tmp = readsimcat(runp.Dfilename,runp.Dftype,runp.Lbox,runp.zspaceaxis,vscale,runp.lg10Mmin2,runp.lg10Mmax2,DorR,&n2);
        }
      } //read halos from the same file.
    
    } //end sim preliminaries.

  (void) time(&time1);
  #ifdef REALLYVERBOSE
  printf("read took %e\n",difftime(time1,time0)); 
  #endif

//don't do this anymore here!
/*
  if(autoorcross == 1) { //merge two catalogs into a single catalog.
    ntot = n1 + n2;
    ctot = (particle *) malloc(sizeof(particle)*(ntot));
    for(i=0;i<n1;i++) {
      ctot[i] = c1tmp[i];
      }
    for(i=0;i<n2;i++) {
      ctot[i+n1] = c2tmp[i];
      }
    free(c1tmp);
    free(c2tmp);
    }
*/

  //remove this eventually (?)
  if(b.maskopt == 1) {
    testmaskcutsmallscale(b);
    printf("passed testmaskcutsmallscale!!\n");
    }

  long double *Npairsfinal;
  Npairsfinal = (long double *) malloc(sizeof(long double)*(b.nbins2d));
  if(runp.radeczorsim == 0) {
    countpairsradecz(c1tmp,n1,c2tmp,n2,mindist,maxdist,b,awgt,Npairsfinal);
    }
  else { //sim.
    //no need to pass originpos for a sim because it will be set to 0.  Just need Lbox.
    countpairssim(c1tmp,n1,c2tmp,n2,runp.Lbox,b,runp.APperp,runp.APpar,Npairsfinal);
    }
  (void) time(&time2);
  #ifdef REALLYVERBOSE
  printf("counts took %e\n",difftime(time2,time1)); 
  #endif

  if(runp.radeczorsim == 0 || b.bintype == 1) {
    printNpairsgeneric(runp.foutbase,Npairsfinal,b,runp.DRopt,n1,n1wgt,n2wgt,runp.binfname,runp.omfid,runp.hfid);
    }
  if(runp.radeczorsim == 1 && b.bintype != 1) { //sims.
    //still need to debug this.
    assert(autoorcross==0);
    printNpairssim(n1,n2,autoorcross,b,runp.Lbox,runp.APperp,runp.APpar,runp.foutbase,Npairsfinal);
    }

  if(runp.radeczorsim == 0 && ((runp.DRopt == 11) || (runp.DRopt == 13))) {
    printNpairsgeneric(runp.foutbase,Npairsfinal,b,runp.DRopt,n1,n1wgt,n2wgt,runp.binfname,runp.omfid,runp.hfid); 
    }
  if(runp.radeczorsim == 0 && ((runp.DRopt == 12) || (runp.DRopt == 14))) {
    printNpairsgeneric(runp.foutbase,Npairsfinal,b,runp.DRopt,n1,n1wgtSI,n2wgt,runp.binfname,runp.omfid,runp.hfid); 
    }

  if(runp.angupweight == 1) {
    freeangwgt(awgt);
    }
  free_rngen();
  freecosmo();
  return 0;
  }


#include "header.h"

xibindat readbinfile(char *frebin) {
  xibindat b;

  //common set-up stuff.
  b.rpimax = -1.; //only for wp.
  b.rperpcut = -1.;  //this set from mask file if there is one. 

  FILE *ifpbin, *ifpmask;
  char line[MAXLINELEN];
  char buf1[MAXLINELEN], buf2[MAXLINELEN], buf3[MAXLINELEN];
  char maskfname[MAXLINELEN];
  char maskfnametot[MAXLINELEN];

  char *r;
  int hbin, lbin;
  int hmask, lmask;
  int h1, l1, h2, l2;
  int i1, i2;

  int blah1, blah2;

  double *r0D, *r2D, *r4D;

  if(!(ifpbin = open_file_read(frebin))) {
    fprintf(stderr,"%s file does not exist.  Exiting...\n",frebin);
    exit(1);
    }

  get_file_length(ifpbin,&hbin,&lbin);
  assert(hbin == 4 || hbin == 5);
  if(hbin == 5) { //then there is a maskfile to read in.
    fgets(line,MAXLINELEN,ifpbin);
    assert(strstr(line,"#"));
    assert(strstr(line,"maskfname: "));
    if(sscanf(line,"%s%s%s",buf1,buf2,buf3)<2) {
      printf("this line does not have properly formatted mask file.  Exiting!\n %s\n",line);
      exit(1);
      }
    if(strstr(buf2,"maskfname:")) {  //leading space
      strcpy(maskfname,buf3);
      }
    else {
      strcpy(maskfname,buf2);
      }
    sprintf(maskfnametot,"%s",maskfname);
    b.maskopt = 1;
    }

  //line 1
  fgets(line,MAXLINELEN,ifpbin);
  assert(strstr(line,"bintype:"));
  if(sscanf(line,"# bintype: %d\n",&(b.bintype)) != 1) {
    printf("error in bintype: binfile header format.  try again\n");
    exit(1);
    }
  assert(b.bintype == 0 || b.bintype == 1 || b.bintype == 2);

  //we need bintype == 2 for velocity statistics.

  //line 2
  fgets(line,MAXLINELEN,ifpbin);
  assert(strstr(line,"xbins:"));
  b.logxopt = -1;
  b.nx = -1;
  b.minx = -1.;
  b.dx = -1.;
  if(sscanf(line,"# xbins: %d %d %lf %lf\n",&(b.logxopt),&(b.nx),&(b.minx),&(b.dx)) != 4) {
    printf("testo %d %d %lf %lf\n",(b.logxopt),(b.nx),(b.minx),(b.dx));
    printf("error in xbin: binfile header format. try again!\n");
    exit(1);
    }
  //line 3
  fgets(line,MAXLINELEN,ifpbin);
  assert(strstr(line,"ybins:"));
  if(b.bintype == 0 || b.bintype == 2) {
    if(sscanf(line,"# ybins: %d\n",&(b.ny)) != 1) {
      printf("error in ybin: binfile header format. try again!\n");
      exit(1);
      }
    b.dy = 1./((double) b.ny);
    b.miny = 0.;
    b.logyopt = 0;
    if(b.bintype == 2) {
      assert(b.ny == 1);
      //set it equal to 4 bins, we're going to fill in all of them!
      //bins are counts, vr, sig2par, sig2perp.
      b.ny = 4;
      }
    }
  else { //xi_grid
    if(sscanf(line,"# ybins: %d %d %lf %lf\n",&(b.logyopt),&(b.ny),&(b.miny),&(b.dy)) != 4) {
      printf("testo %d %d %lf %lf\n",(b.logyopt),(b.ny),(b.miny),(b.dy));
      printf("error in ybin: binfile header format. try again!\n");
      exit(1);
      }
    }
  //line 4
  fgets(line,MAXLINELEN,ifpbin);
  assert(strstr(line,"rebinopt:"));
  if(sscanf(line,"# rebinopt: %d\n",&(b.rebinopt)) != 1) {
    printf("error in ybin: binfile header format. try again!\n");
    exit(1);
    }

  b.nbins2d = b.nx*b.ny;

  printf("here are the xiell x settings: %d %d %lf %lf\n",(b.logxopt),(b.nx),(b.minx),(b.dx));
  printf("here are the xiell y settings: %d %d %lf %lf\n",(b.logyopt),(b.ny),(b.miny),(b.dy));
  printf("xiell rebinopt: %d\n",b.rebinopt);
  int nchk,ncnt,mubin;

  if(b.rebinopt != 0) {  //in this case, we expect rebinopt
    //read in the rebin information to follow
    b.ndata = lbin - hbin;
    (b.binxiell) = (int *) malloc(sizeof(int)*(b.ndata));
    //read through the file once, get nr0, b.nr2, b.nr4
    b.nr0 = 0;
    b.nr2 = 0;
    b.nr4 = 0;
    for(i1=0;i1<lbin-hbin;i1++) {
      fscanf(ifpbin,"%d %d %d\n",&nchk,&blah1,&blah2);
      if(nchk == 0) {
        b.nr0 += 1;
        b.ellmaxdata = 0;
        }
      if(nchk == 2) {
        b.nr2 += 1;
        b.ellmaxdata = 2;
        }
      if(nchk == 4) {
        b.nr4 += 1;
        b.ellmaxdata = 4;
        }
      printf("%d %d %d %d %d\n",nchk,b.nr0,b.nr2,b.nr4,blah1);
      }
    (void) rewind(ifpbin);
    for(i1=0;i1<hbin;i1++) {
      fgets(line,MAXLINELEN,ifpbin);
      }

    r0D = (double *) calloc(b.nr0,sizeof(double));
    r2D = (double *) calloc(b.nr2,sizeof(double));
    r4D = (double *) calloc(b.nr4,sizeof(double));

    ncnt = 0;
    printf("%d %d %d %d\n",b.ndata,b.nr0,b.nr2,b.nr4);
    assert((b.ndata) == (b.nr0)+(b.nr2)+(b.nr4));
    for(i1=0;i1<(b.nr0);i1++) {
      fscanf(ifpbin,"%d %d %d\n",&nchk,&((b.binxiell)[i1]),&mubin);
      assert(mubin == 200/b.ny);
      assert(200%mubin == 0);
      assert(200%b.ny == 0);
      assert(nchk == 0);
      ncnt += (b.binxiell)[i1];
      }
    assert(ncnt <= b.nx);
    ncnt = 0;
    for(i1=0;i1<(b.nr2);i1++) {
      fscanf(ifpbin,"%d %d %d\n",&nchk,&((b.binxiell)[(b.nr0)+i1]),&mubin);
      assert(mubin == 200/b.ny);
      assert(200%mubin == 0);
      assert(200%b.ny == 0);
      assert(nchk == 2);
      ncnt += (b.binxiell)[(b.nr0)+i1];
      }
    assert(ncnt <= b.nx);  //all the observed x bins are counted.
    ncnt = 0;
    for(i1=0;i1<(b.nr4);i1++) {
      fscanf(ifpbin,"%d %d %d\n",&nchk,&((b.binxiell)[(b.nr0)+(b.nr2)+i1]),&mubin);
      assert(mubin == 200/b.ny);
      assert(200%mubin == 0);
      assert(200%b.ny == 0);
      assert(nchk == 4);
      ncnt += (b.binxiell)[(b.nr0)+(b.nr2)+i1];
      }
    assert(ncnt <= b.nx);  //all the observed x bins are counted.
    }
  fclose(ifpbin);

  //read in mask file.
  if(b.maskopt == 1) {
    if(!(ifpmask = open_file_read(maskfnametot))) {
      fprintf(stderr,"%s file does not exist.  Exiting...\n",maskfnametot);
      exit(1);
      }
    get_file_length(ifpmask,&hmask,&lmask);
    assert(lmask-hmask == b.nbins2d);
    (b.mask) = (int *) malloc(sizeof(int)*(b.nbins2d));
    assert(hmask == 1);
    for(i1=0;i1<hmask;i1++) {
      fgets(line,MAXLINELEN,ifpmask);
      assert(strstr(line,"#"));
      assert(strstr(line,"rperpcut: "));
      r = strchr(line,':');
      assert(r);
      b.rperpcut = atof(r+1);
      }
    for(i1=0;i1<lmask;i1++) {
      fscanf(ifpmask,"%d\n",&((b.mask)[i1]));
      }
    }
  else {
    b.maskopt = 0;
    (b.mask) = NULL;
    }
  return b;
  }

int testmaskcutsmallscale(xibindat b) {
  int ix,iy,mi;
  double rminbin,mumax,rperpmin;
  mi = 0;
  for(ix=0;ix<b.nx;ix++) {
    rminbin = ((long double) pow(10.,b.minx+(ix)*b.dx));
//    rmaxbin = ((long double) pow(10.,b.minx+(ix+1)*b.dx));
    for(iy=0;iy<b.ny;iy++) {
      mumax = b.dy*(iy+1);
      rperpmin = rminbin*sqrt(1.-mumax*mumax);
      if(rperpmin < b.rperpcut) {
        assert(b.mask[mi] == 1);
        }
      else {
        assert(b.mask[mi] == 0);
        }
      mi += 1;
      }
    }
  return 0;
  }


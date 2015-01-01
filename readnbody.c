#include "header.h"

particle *readsimcat(char *fname, int whichcat, real Lbox, int zspaceaxis, real vscale, float Mmin, float Mmax, int DorR, int *ngals) {

  double minP[3];
  double maxP[3];

  real tmppos[3];
  real tmpvel[3];

  int i1,i2,i3;
  float f1,f2;
  float vel[3];
  real myzspacepos;


  minP[0] = minP[1] = minP[2] = 10000.;
  maxP[0] = maxP[1] = maxP[2] = -10000.;

  int linecount, headercount;
  int i,j;
  FILE *ifp;
  char line[MAXLINELEN];

  int mycnt = 0; //count halos that pass the mass cut.
  float lg10Mfof, lg10M180b;

  ifp = open_file_read(fname);
  get_file_length(ifp,&headercount,&linecount);
  if(whichcat == 0) {
    assert(headercount == 0);
    (*ngals) = linecount;
    }
  if(whichcat == 1 || whichcat == 4 || whichcat == 5) {
    assert(headercount == 0 || headercount == 1 || headercount == 5);
    for(i=0;i<headercount;i++) {
      fgets(line,MAXLINELEN,ifp);
      }
    mycnt = 0;
    for(i=0;i<linecount;i++) {
      #ifdef REALFLOAT
      fscanf(ifp,"%E %E %E %E %E %E\n",&(tmppos[0]),&(tmppos[1]),&(tmppos[2]),&(tmpvel[0]),&(tmpvel[1]),&(tmpvel[2]));
      #else
      fscanf(ifp,"%lE %lE %lE %lE %lE %lE\n",&(tmppos[0]),&(tmppos[1]),&(tmppos[2]),&(tmpvel[0]),&(tmpvel[1]),&(tmpvel[2]));
      #endif
      if(whichcat == 1) {
        fscanf(ifp,"%d %f\n",&i1,&lg10Mfof);
        if(lg10Mfof >= Mmin && lg10Mfof <= Mmax) {
          mycnt += 1;
          }
        }
      if(whichcat == 4 || whichcat == 5) {
        fscanf(ifp,"%d %f %f\n",&i1,&lg10Mfof,&lg10M180b);
        if(whichcat == 4) {
          if(lg10Mfof >= Mmin && lg10Mfof <= Mmax) {
            mycnt += 1;
            }
          }  //end 4
        else { //whichcat 5
          if(lg10M180b >= Mmin && lg10M180b <= Mmax) {
            mycnt += 1;
            }
          } //end whichcat 5
        } //whichcat 4 or 5
      } //end for loop
    (void) rewind(ifp);
    (*ngals) = mycnt;
    } //145
  if(whichcat == 2) {
    assert(headercount == 0);
    (*ngals) = linecount; 
    }
  if(whichcat == 3) {
    assert(headercount == 0);
    (*ngals) = linecount; 
    }
  if(whichcat == 6) {
    assert(headercount == 2);
    (*ngals) = linecount - headercount; 
    }

  particle *c1 = (particle *) malloc(sizeof(particle)*(*ngals));
  #ifdef REALLYVERBOSE 
  printf("ahh beth %d %d %d %d\n",linecount,headercount,*ngals,whichcat);
  #endif
  int mycntchk = 0;

  int iindx = 0;

  if(whichcat == 2) {
    #ifdef VOPT
    fprintf(stderr,"VOPT should not be defined with whichcat == 2!\n");
    exit(1);
    #endif
    }

  for(i=0;i<headercount;i++) {
    fgets(line,MAXLINELEN,ifp);
    }
  for(i=0;i<linecount-headercount;i++) {
    if(i%1000000 == 0) {
      #ifdef REALLYVERBOSE 
      printf("finished %d of %d\n",i,linecount-headercount);
      #endif
      }
    if(whichcat == 2) {
      #ifdef REALFLOAT
      fscanf(ifp,"%E %E %E %E %E %E %d %d %d\n",&(c1[iindx].pos[0]),&(c1[iindx].pos[1]),&(c1[iindx].pos[2]),&(myzspacepos),&f1,&f2,&i1,&i2,&i3);
      #else
      fscanf(ifp,"%lE %lE %lE %lE %E %E %d %d %d\n",&(c1[iindx].pos[0]),&(c1[iindx].pos[1]),&(c1[iindx].pos[2]),&(myzspacepos),&f1,&f2,&i1,&i2,&i3);
      #endif 
      mycntchk += 1;
      if(zspaceaxis >= 0) {
        c1[iindx].pos[zspaceaxis] = fixbox(myzspacepos,Lbox);
        }
      for(j=0;j<=2;j++) {
        assert(c1[iindx].pos[j] >= 0. && c1[iindx].pos[j] <= Lbox);
        minP[j] = min(c1[iindx].pos[j],minP[j]);
        maxP[j] = max(c1[iindx].pos[j],maxP[j]);
        }
      } //2
    else {
      //JT mock challenge catalog.
      if(whichcat == 6) {
        #ifdef REALFLOAT
        fscanf(ifp,"%f %f %f %f %f %f\n",&(c1[iindx].pos[0]),&(c1[iindx].pos[1]),&(c1[iindx].pos[2]),&(tmpvel[0]),&(tmpvel[1]),&(tmpvel[2]));
        #else
        fscanf(ifp,"%lf %lf %lf %lf %lf %lf\n",&(c1[iindx].pos[0]),&(c1[iindx].pos[1]),&(c1[iindx].pos[2]),&(tmpvel[0]),&(tmpvel[1]),&(tmpvel[2]));
        #endif
        for(j=0;j<=2;j++) {
          fixbox(c1[iindx].pos[j],Lbox);
          minP[j] = min(c1[iindx].pos[j],minP[j]);
          maxP[j] = max(c1[iindx].pos[j],maxP[j]);
          }
        mycntchk += 1;
        }
      else {
        #ifdef REALFLOAT
        fscanf(ifp,"%E %E %E %E %E %E\n",&(c1[iindx].pos[0]),&(c1[iindx].pos[1]),&(c1[iindx].pos[2]),&(tmpvel[0]),&(tmpvel[1]),&(tmpvel[2]));
        #else
        fscanf(ifp,"%lE %lE %lE %lE %lE %lE\n",&(c1[iindx].pos[0]),&(c1[iindx].pos[1]),&(c1[iindx].pos[2]),&(tmpvel[0]),&(tmpvel[1]),&(tmpvel[2]));
        #endif
        #ifdef VOPT
        for(ii=0;ii<=2;ii++) {
          c1[iindx].vel[ii] = tmpvel[ii];
          }
        #endif
        if(whichcat == 0) {
          fscanf(ifp,"%d %d %f\n",&i1,&i2,&f1);
          mycntchk += 1;
          } //0
        if(whichcat == 1) {
          fscanf(ifp,"%d %f\n",&i1,&lg10Mfof);
          if(lg10Mfof >= Mmin && lg10Mfof <= Mmax) {
            mycntchk += 1;
            }
          else {
            continue;
            }
          } //1
        if(whichcat == 4) {
          fscanf(ifp,"%d %f %f\n",&i1,&lg10Mfof,&lg10M180b);
          if(lg10Mfof >= Mmin && lg10Mfof <= Mmax) {
            mycntchk += 1;
            }
          else {
            continue;
            }
          } //4
        if(whichcat == 3) {
          fscanf(ifp,"%f %d\n",&f1,&i1);
          mycntchk += 1;
          } // 3
        if(whichcat == 5) {
          fscanf(ifp,"%d %f %f\n",&i1,&lg10Mfof,&lg10M180b);
          if(lg10M180b >= Mmin && lg10M180b <= Mmax) {
            mycntchk += 1;
            }
          else {
            continue;
            }
          } //5
        } //else whichcat != 6
      } //else !2
    if(whichcat == 0 || whichcat == 1 || whichcat == 4 || whichcat == 5) {
      if(zspaceaxis >= 0) {
        //c1[iindx].pos[zspaceaxis] = c1[iindx].pos[zspaceaxis] + c1[iindx].vel[zspaceaxis];
        c1[iindx].pos[zspaceaxis] = c1[iindx].pos[zspaceaxis] + tmpvel[zspaceaxis];
      //  c1[iindx].pos[zspaceaxis] = fixbox(c1[iindx].pos[zspaceaxis]*Lbox,Lbox);
        }
      for(j=0;j<=2;j++) {
        if(whichcat == 0 || whichcat == 1 || whichcat == 4 || whichcat == 5) {
          c1[iindx].pos[j] = fixbox(c1[iindx].pos[j]*Lbox,Lbox);
          #ifdef VOPT
          c1[iindx].vel[j] = c1[iindx].vel[j]*Lbox;
          #endif
          }
        minP[j] = min(c1[iindx].pos[j],minP[j]);
        maxP[j] = max(c1[iindx].pos[j],maxP[j]);
        }
      } //0,1,4,5
    if(whichcat == 3) {
      if(zspaceaxis >= 0) {
        c1[iindx].pos[zspaceaxis] = c1[iindx].pos[zspaceaxis] + tmpvel[zspaceaxis]*vscale;
        c1[iindx].pos[zspaceaxis] = fixbox(c1[iindx].pos[zspaceaxis],Lbox);
        assert(c1[iindx].pos[zspaceaxis] >= 0. && c1[iindx].pos[zspaceaxis] <= Lbox);
        }
      for(j=0;j<=2;j++) {
        assert(c1[iindx].pos[j] >= 0. && c1[iindx].pos[j] <= Lbox);
        #ifdef VOPT
        c1[iindx].vel[j] = c1[iindx].vel[j]*vscale;
        #endif
        minP[j] = min(c1[iindx].pos[j],minP[j]);
        maxP[j] = max(c1[iindx].pos[j],maxP[j]);
        }
      } //3
    c1[iindx].weight = 1.0;
//    c1[iindx].DorR = DorR;
/*
    for(ii=0;ii<=2;ii++) {
      ipos[ii] = ((int) floor(gallist[gindx].pos[ii]*Ncell/Lbox));
      posmin[ii] = min(posmin[ii],gallist[gindx].pos[ii]);
      posmax[ii] = max(posmax[ii],gallist[gindx].pos[ii]);
      }

    c1[iindx].cell = i2n(ipos);
*/
    iindx += 1;
    if(iindx == (*ngals)) {
      break;
      }
    }  //end for loop.
  printf("%d %d\n",mycntchk,iindx);
  assert(mycntchk == iindx);
  fclose(ifp);
  #ifdef REALLYVERBOSE
  printf("this many gals,%d\n",*ngals);
  #endif

  if(whichcat == 6) {
    printf("Jeremy cat min/max\n");
    printf("%f %f %f\n",minP[0],minP[1],minP[2]);
    printf("%f %f %f\n",maxP[0],maxP[1],maxP[2]);
    }

  return c1; 
  }


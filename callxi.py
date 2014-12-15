import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import os

def writeparamfile(pfname,runp,DRopt):
  """
  runp is a dictionary with all hte relevant parameters.
  """

  mykeylist = ['omfid','hfid','binfname','zmin','zmax','foutbase','Dfilename','Rfilename']
  defaultlist = [0.292, 0.69, "xibinfiles/bin_xismu1Mpcohtest.txt", 0.43, 0.7, "outputdr12/testing","tmp",'tmp']

  ofp = open(pfname,'w')
  ## things we're holding fixed for now.
  ofp.write('DRopt = %d\n' % DRopt)
  ofp.write('radeczorsim = 0\n')
  ofp.write('unitsMpc = 0\n')
  ofp.write('Dftype = 1\n') 
  ofp.write('Rftype = 1\n')
 
  for kk in mykeylist:
    try:
      aa = runp[kk]
    except:
      if kk is not 'omfid' and kk is not 'hfid':
        print 'key not found, aborting!'
        sys.exit(1)  

    if type(runp[kk]) is str:
      ofp.write('%s = %s\n',kk,runp[kk])
    if type(runp[kk]) is int:
      ofp.write('%s = %d\n',kk,runp[kk])
    if type(runp[kk]) is float:
      ofp.write('%s = %e\n',kk,runp[kk])

  ofp.close()

def callxi(sample):
  pass

def callxisimcat(pfname,runp,fullpathxi='./xi',runxi=0):

  ofp = open(pfname,'w')
  ## things we're holding fixed.
  ofp.write('radeczorsim = 1\n')

  #mykeylist = ['omfid','hfid','unitsMpc','binfname','foutbase','Dftype','Dfilename','Lbox','abox','zpsaceaxis']
  #defaultlist = [0.292, 0.69, 0, "/home/howdiedoo/xibinfiles/bin_xismu5Mpcoh.txt",\ 
  #               "testo",2,"PB00BF_LOWZv0arg_PB00_0.8000_P0V1_HV1.000_IHV1.000_CENV0.000_testo_SLOW.cat",\
  #               1380.0, 0.8000, 2]

  ddict = {'omfid':0.292, 'hfid':0.69, 'unitsMpc':0, 'binfname': '/home/howdiedoo/gitprojects/xi/xibinfiles/bin_xismu4Mpcoh.txt', 'foutbase':'testo','Dftype':2, 'Dfilename': '/home/howdiedoo/boss/bethalexieDR12/PB00BF_LOWZv0arg_PB00_0.8000_P0V1_HV1.000_IHV1.000_CENV0.000_testo_SLOW.cat', 'Lbox':1380.0, 'abox':0.8000, 'zspaceaxis':2}


  for kk in ddict.keys():
    try:
      aa = runp[kk]
    except:
      if kk is 'omfid' or kk is 'hfid' or kk is 'Dfilename' or kk is 'foutbase':
        print 'key not found, aborting!',kk
        sys.exit(1)  
      ## otherwise, just assume default!
      aa = ddict[kk]

    if type(ddict[kk]) is str:
      ofp.write('%s = %s\n' % (kk,aa))
    if type(ddict[kk]) is int:
      ofp.write('%s = %d\n' % (kk,aa))
    if type(ddict[kk]) is float:
      ofp.write('%s = %e\n' % (kk,aa))

  ## finally, add optional stuff.
  for kk in runp:
    if kk in ddict.keys(): continue
    ofp.write('%s = %s\n' % (kk,runp[kk]))


  ofp.close()
  if runxi == 1:
    os.system('%s %s' % (fullpathxi,pfname))


if __name__ == '__main__':

 ## data directory locations.  Later add riemann and nersc.  copy gitprojects/LSSanalysis/
 ## whichtask
 ## datadir (I should know everything else based on cmass/lowz, N/S, ang/targ, etc.

  runp = {}
  ## these are RunPB defaults.
  runp['omfid'] = 0.292 
  runp['hfid'] = 0.69

  if len(sys.argv) != 1:
    print 'Usage: python paircountwrapper.py sampletag NorS whichtask runopt' 
    print 'ASSUMES OMFID, HFID VALUES INSIDE THIS CODE!  ALTER IF SO DESIRED.'
    print 'omfid = ',runp['omfid']
    print 'hfid = ',runp['hfid']
    sys.exit(1)

  dHH = "/home/howdiedoo/boss/mksamplecatslatestdr12/"
  
  ## test for the sims.
  dfile = '/home/howdiedoo/boss/bethalexieDR12/PB00BF_LOWZv0arg_PB00_0.8000_P0V1_HV1.000_IHV1.000_CENV0.000_testo_SLOW.cat'
  foutbase = dfile + '.testo'
  myd = {'omfid':0.292, 'hfid':0.69, 'Dfilename': dfile, 'foutbase': foutbase, 'Lbox':1380.0, 'abox':0.8000, 'zspaceaxis':2}
  callxisimcat('tmprunxi.params',myd)




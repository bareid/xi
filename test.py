import numpy as np
from ctypes import Structure, POINTER, c_int, c_double, c_void_p, CDLL

# Here's our C structure defined in Python
class xibininfo(Structure):
  _fields_ = ("bintype",c_int), ("logxopt",c_int), ("nx",c_int), ("minx",c_double), ("dx",c_double),\
                                ("logyopt",c_int), ("ny",c_int), ("miny",c_double), ("dy",c_double),\
                                ("realorzspace",c_int), ("periodicopt",c_int), ("nbins2d",c_int)

class dset(Structure):
  _fields_ = ("px",c_void_p), ("py",c_void_p), ("pz",c_void_p),\
             ("vx",c_void_p), ("vy",c_void_p), ("vz",c_void_p),\
             ("wgt",c_void_p), ("np",c_int)

test = CDLL("test.so")

## read in numpy array.
### ack these give me things that are offset.
px,py,pz,w1 = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[0,1,2,3],unpack=True)
vx,vy,vz = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[0,1,2],unpack=True)

print 'well this?'
test.testsingle(px.ctypes.data_as(c_void_p),len(px))

nx = 20
ny = 200
nbins2d = nx*ny
tt = xibininfo(0,0,nx,0.0,1.5,0,ny,0.,1.,0,1,nbins2d)

test.printit(tt)

#d1 = dset(px.ctypes.data_as(c_void_p),py.ctypes.data_as(c_void_p),pz.ctypes.data_as(c_void_p),None,None,None,w1.ctypes.data_as(c_void_p),len(px))

ppx = np.arange(0.,100.,1.,dtype='float64')
ppy = np.arange(0.,100.,1.,dtype='float64') + 1000
ppz = np.arange(0.,100.,1.,dtype='float64') + 2000

vvx = np.arange(0.,100.,1.,dtype='float64') + 3000
vvy = np.arange(0.,100.,1.,dtype='float64') + 4000
vvz = np.arange(0.,100.,1.,dtype='float64') + 5000

ww = np.ones(len(ppx))*-1000 + np.arange(0.,100.,1.,dtype='float64')

d1 = dset(ppx.ctypes.data_as(c_void_p),ppy.ctypes.data_as(c_void_p),ppz.ctypes.data_as(c_void_p),\
          vvx.ctypes.data_as(c_void_p),vvy.ctypes.data_as(c_void_p),vvz.ctypes.data_as(c_void_p),\
          ww.ctypes.data_as(c_void_p),len(ppx))

test.hellop(d1)
ppx[:] = px[:len(ppx)]
ppy[:] = py[:len(ppx)]
ppz[:] = pz[:len(ppx)]

vvx[:] = vx[:len(ppx)]
vvy[:] = vy[:len(ppx)]
vvz[:] = vz[:len(ppx)]

ww[:] = w1[:len(ppx)]

d1 = dset(ppx.ctypes.data_as(c_void_p),ppy.ctypes.data_as(c_void_p),ppz.ctypes.data_as(c_void_p),\
          vvx.ctypes.data_as(c_void_p),vvy.ctypes.data_as(c_void_p),vvz.ctypes.data_as(c_void_p),\
          ww.ctypes.data_as(c_void_p),len(ppx))
print 'does this work better?'

d2 = dset(px.ctypes.data_as(c_void_p),py.ctypes.data_as(c_void_p),pz.ctypes.data_as(c_void_p),\
          vx.ctypes.data_as(c_void_p),vy.ctypes.data_as(c_void_p),vz.ctypes.data_as(c_void_p),\
          w1.ctypes.data_as(c_void_p),len(px))

print 'type stuff.',px.dtype, ppx.dtype
print 'type stuff.',vx.dtype, vvx.dtype

print 'shape stuff.',px.shape, ppx.shape
print 'shape stuff.',vx.shape, vvx.shape

test.hellop(d1)
# native outputs from np.loadtxt are screwed up!
test.hellop(d2)
print 'that is broken, what if we read one column at a time?'

px = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[0])
py = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[1])
pz = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[2])
vx = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[0])
vy = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[1])
vz = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[2])
www = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[3])

print 'does this work better? yes, this seems to work.  now try single numpy array?  I think we need to worry about the memory being contiguous or something.'
d3 = dset(px.ctypes.data_as(c_void_p),py.ctypes.data_as(c_void_p),pz.ctypes.data_as(c_void_p),\
          vx.ctypes.data_as(c_void_p),vy.ctypes.data_as(c_void_p),vz.ctypes.data_as(c_void_p),\
          www.ctypes.data_as(c_void_p),len(px))
test.hellop(d3)



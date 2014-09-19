from ctypes import Structure, POINTER, c_int, c_double, CDLL

# Here's our C structure defined in Python
class xibininfo(Structure):
  _fields_ = ("bintype",c_int), ("logxopt",c_int), ("nx",c_int), ("minx",c_double), ("dx",c_double),\
                                ("logyopt",c_int), ("ny",c_int), ("miny",c_double), ("dy",c_double),\
                                ("realorzspace",c_int), ("periodicopt",c_int), ("nbins2d",c_int)

  nx = 20
  ny = 200
  nbins2d = nx*ny
  tt = xibininfo(0,0,nx,0.0,1.5,0,ny,0.,1.,0,1,nbins2d)

  ## read in numpy array.
  p1 = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[0,1,2])
  w1 = np.loadtxt("/Users/bareid/gitprojects/BIDS/xi/usecase/PB00_0.6452.halos.short",usecols=[3])


  test = CDLL("test.so")
  test.printit(tt)





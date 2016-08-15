import pylab as pl
import numpy as np
pl.rc('legend',fontsize=20)


data=np.loadtxt('2L_180.dat',skiprows=0)

x=data[0,:]
h1= data[1,:]

h2= data[2,:]
h3= data[3,:]
h4= data[4,:]
h5= data[5,:]
h6= data[6,:]
h7= data[7,:]
h8= data[8,:]
#h9= data[9,:]
#h10= data[10,:]
pl.figure(1,figsize=(5,3))
pl.plot(
            x,h1,x,h2,
            x,h3,x,h4,x,h5,
    x,h6,x,h7,
    x,h8,
#    x,
#    h9,x,h10,
    'k-',linewidth=0.8
#    ,label=r'time =  '
    )
#pl.legend(loc=1)
pl.xlabel(r'x',fontsize=14)
pl.ylabel(r'h',fontsize=14)
#pl.xlim(0.0,0.08)
#pl.ylim(0.0,30)
pl.locator_params(axis = 'x', nbins = 20)
pl.locator_params(axis = 'y', nbins = 10)
pl.tight_layout()
#pl.savefig("graph@time=.pdf")
pl.show()

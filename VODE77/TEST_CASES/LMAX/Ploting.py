import pylab as pl
import numpy as np
pl.rc('legend',fontsize=20)


data=np.loadtxt('check.dat',skiprows=0)

x=data[0,:]
h1= data[1,:]

h2= data[2,:]
h3= data[3,:]
h4= data[4,:]
h5= data[5,:]
h6= data[6,:]
h7= data[7,:]
h8= data[8,:]
h9= data[9,:]
h10= data[10,:]
pl.figure(1,figsize=(5,3))
pl.plot(x,h1,x,h2,x,h3,x,h4,x,h5,
    x,h6,x,h7,x,h8,x,h9,x,h10,
    'k-',linewidth=2.5
#    ,label=r'time =  '
    )
#pl.legend(loc=1)
pl.xlabel(r'X',fontsize=20)
pl.ylabel(r'H',fontsize=20)
#pl.xlim(0.0,0.08)
pl.xlim(0.0,9.0)
pl.locator_params(axis = 'x', nbins = 7)
pl.locator_params(axis = 'y', nbins = 5)
pl.tight_layout()
pl.savefig("fig.pdf")
pl.show()

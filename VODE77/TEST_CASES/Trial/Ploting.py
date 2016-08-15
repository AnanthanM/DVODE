import pylab as pl
import numpy as np
pl.rc('legend',fontsize=12)


data=np.loadtxt('sixty.dat',skiprows=0)

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
pl.figure(1,figsize=(6,3))
#pl.plot(x,h1,x,h2,x,h3,x,h4,x,h5,
#    x,h6,x,h7,x,h8,x,h9,x,h10,
#    'k-',linewidth=1.0
#    )
pl.plot(x,h1,'--',linewidth=1.5, label="t=0" )
pl.plot(x,h2,'--',linewidth=1.5, label="t=15" )
pl.plot(x,h3,'--',linewidth=1.5, label="t=20" )
pl.plot(x,h4,'--',linewidth=1.5, label="t=25" )
pl.plot(x,h5,'--',linewidth=1.5, label="t=27" )
pl.plot(x,h6,linewidth=1.5, label="t=28" )
pl.plot(x,h7,linewidth=1.5, label="t=28.1" )
pl.plot(x,h8,linewidth=1.5, label="t=28.16" )
pl.plot(x,h9,linewidth=1.5, label="t=28.2" )
pl.plot(x,h10,linewidth=1.5, label="t=28.5" )
#pl.legend(loc=1)
pl.xlabel(r'X',fontsize=14)
pl.ylabel(r'H',fontsize=14)
pl.xticks(fontsize=12)
pl.yticks(fontsize=12)
pl.xlim(0,9)
#pl.ylim(0.0,30)
leg3 = pl.legend(ncol=2,bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.70 )
pl.locator_params(axis = 'x', nbins = 7)
pl.locator_params(axis = 'y', nbins = 5)
pl.tight_layout()
pl.savefig("2_2pi.pdf", bbox_extra_artists=(leg3,), bbox_inches='tight')
pl.show()

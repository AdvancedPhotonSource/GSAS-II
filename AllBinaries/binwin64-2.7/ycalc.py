import os
import numpy as np
import numpy.ma as ma
#import GSASIIpath # needed if pypowder is in a ./bin... directory
#import GSASIIpwd as G2pwd
import cPickle
import pypowder as pyd

#def getFCJVoigt3(pos,sig,gam,shl,xdata):
#    'needs a doc string'
#    Df = pyd.pypsvfcj(len(xdata),xdata-pos,pos,sig,gam,shl)
#    Df /= np.sum(Df)
#    return Df
loc = os.path.split(os.path.abspath(__file__))[0]
fp = open(os.path.join(loc,'yc_info.CPickle'),'r')
x = cPickle.load(fp)
yc = np.zeros_like(x)
while(True): # loop on phases
    l = cPickle.load(fp)
    if len(l) != 2:
        yc_orig = l
        break
    else:
        im,shl = l
    while(True): # loop on reflections
        l = cPickle.load(fp)
        if l is None: break
        refl,iBeg,iFin = l
        print iBeg,iFin
        #yc[iBeg:iFin] += refl[11+im]*refl[9+im]*getFCJVoigt3(refl[5+im],refl[6+im],refl[7+im],shl,ma.getdata(x[iBeg:iFin]))    #>90% of time spent here
        xdata = ma.getdata(x[iBeg:iFin])-refl[5+im]
        #def getFCJVoigt3(pos,sig,gam,shl,xdata):
        Df = pyd.pypsvfcj(len(xdata),xdata,refl[5+im],refl[6+im],refl[7+im],shl)
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]* Df / np.sum(Df)
        print np.sum(Df)

fp.close()
import pylab
pl = pylab.figure(facecolor='white')
ax = pl.add_subplot(111)
ax.set_xlabel(r'$\mathsf{2\theta}$')
ax.set_ylabel('Intensity')
#ax.set_title(Title)
ax.plot(x,yc,label='calc')
ax.plot(x,yc_orig,label='pickle')
ax.plot(x,yc-yc_orig,label='diff')
ax.legend()
pylab.show()

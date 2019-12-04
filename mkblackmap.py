from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.time import Time
import numpy as np
import time,getopt,sys
import healpy as hp
import matplotlib.pyplot as plt
from pylab import cm

rootname='nov2019' # default rootname

def Usage():
    print("usage: python mkblackmap.py -h -r rootname (string)")

def plotblackmap(rootname):

    expMap = hp.read_map(rootname+'.fits')
    smap=expMap
    smap=np.ma.masked_where(smap==0,smap)
    smap.set_fill_value(np.nan)

    cmap = cm.jet
    cmap.set_under('k')
    cmap.set_bad('xkcd:charcoal')
    plt.style.use('dark_background')
    plt.rcParams.update({'font.size':15})
    lon=np.arange(360)
    lat=np.zeros(360)

    hp.mollview(smap,cmap=cmap,min=0.0,max=3e4,rot=-80,
                flip='geo',coord='C',cbar=False,notext=True,
                norm='linear',xsize=1000,title=' ')
    hp.projplot(lon,lat,'r',lonlat=True,coord='G')
    hp.graticule(dpar=45.,dmer=30.,coord='C')
    pngfile1=rootname+'_black.png'
    plt.savefig(pngfile1,dpi=600)
    plt.close()

def main(rootname):

    plotblackmap(rootname)

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hr:",["help","rootname="])

    except getopt.GetoptError:
        print("Error processing arguments:")

    for opt,arg in opts:
        if opt in ("-h","--help"):
            Usage()
            sys.exit()
        elif opt in ("-r","--rootname"):
            rootname = arg
    
    main(rootname)

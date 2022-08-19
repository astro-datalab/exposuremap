from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.time import Time
import numpy as np
import time,sys,getopt
import healpy as hp
import matplotlib.pyplot as plt
from pylab import cm

rootname = 'aug2022' # default rootname
yearStartEnd = (2005,2022)

def Usage():
    print("usage: python mkyears.py -h -r rootname (string) -y years (tuple)")

def plotyears(rootname,yearStartEnd):

    expTable = Table.read(rootname+'_data.fits')
    indTable = Table.read(rootname+'_index.fits')
    expMap = hp.read_map(rootname+'.fits')

    if 'mjd_obs' in expTable.colnames:
        mjd = expTable['mjd_obs']
    elif 'MJD-OBS' in expTable.colnames:
        mjd = expTable['MJD-OBS']
    else:
        sys.exit('No MJD column')

    keep=np.isfinite(mjd)
    expTable=expTable[keep]
    indTable=indTable[keep]
    mjd=mjd[keep]
    srt=np.argsort(mjd)
    expTable=expTable[srt]
    indTable=indTable[srt]
    mjd=mjd[srt]


    tt = Time(mjd, format='mjd')
    yy = tt.decimalyear

    #plt.style.use('dark_background')
    cmap = cm.jet
    cmap.set_bad('xkcd:charcoal')
    cmap.set_under('w')
    plt.rcParams.update({'font.size':15})
    lon=np.arange(360)
    lat=np.zeros(360)

    years = np.arange(yearStartEnd[0],yearStartEnd[1]+1)
    els = np.arange(len(expTable))

    expMap1 = np.zeros(len(expMap))

    for year in years:
        
        yels = els[(yy < year) & (yy >= (year-1))]
        tobj=tt[yels[-1]]
        title=tobj.datetime.strftime('%Y')
        for el in yels:
            mels=indTable[el]['ind'][0:indTable[el]['np']]
            expMap1[mels]=expMap1[mels]+expTable[el]['exposure']

        smap=expMap1.copy()
        smap=np.ma.masked_where(smap==0,smap)
        smap.set_fill_value(np.nan)
        hp.mollview(np.log10(smap),cmap=cmap,min=0.0,max=6.,rot=-80,
                flip='geo',coord='C',cbar=True,notext=True,
                norm='linear',xsize=1000,title=title,unit='$\mathrm{log_{10} sec}$')
        hp.projplot(lon,lat,'r',lonlat=True,coord='G')
        hp.graticule(dpar=45.,dmer=30.,coord='C')
        pngfile1='log_expmap_'+title+'.png'
        plt.savefig(pngfile1,dpi=300)
        plt.close()


def main(rootname,yearStartEnd):

    plotyears(rootname,yearStartEnd)

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hr:y:",["help","rootname=","years="])

    except getopt.GetoptError:
        print("Error processing arguments:")

    for opt,arg in opts:
        if opt in ("-h","--help"):
            Usage()
            sys.exit()
        elif opt in ("-r","--rootname"):
            rootname = arg
        elif opt in ("-y","--years"):
            yearStartEnd = arg

    
    main(rootname,yearStartEnd)

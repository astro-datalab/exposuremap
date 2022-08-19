from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.time import Time
import numpy as np
import time
import healpy as hp
import matplotlib.pyplot as plt
from pylab import cm

rootname='aug2022'
expTable = Table.read(rootname+'_data.fits')
indTable = Table.read(rootname+'_index.fits')
expMap = hp.read_map(rootname+'.fits')
expMap1 = np.zeros(len(expMap))

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

tstep=7. # days
t0=np.min(mjd)
tf=np.max(mjd)

week_bin=np.trunc((mjd-t0)/tstep)
weeks=(expTable.group_by(week_bin)).groups.indices

plt.style.use('dark_background')
cmap = cm.jet
cmap.set_bad('xkcd:charcoal')
cmap.set_under('k')
plt.rcParams.update({'font.size':15})
lon=np.arange(360)
lat=np.zeros(360)

for i in np.arange(len(weeks)-1):
    el1=weeks[i]
    el2=weeks[i+1]
    els=np.arange(el1,el2)
    t1=mjd[el1]
    tobj=Time(t1,format='mjd')
    title=tobj.datetime.strftime('%Y')
    for el in els:
        mels=indTable[el]['ind'][0:indTable[el]['np']]
        expMap1[mels]=expMap1[mels]+expTable[el]['exposure']

    smap=expMap1.copy()
    smap=np.ma.masked_where(smap==0,smap)
    smap.set_fill_value(np.nan)
    hp.mollview(np.log10(smap),cmap=cmap,min=0.0,max=6.,rot=-80,
                flip='geo',coord='C',cbar=False,notext=True,
                norm='linear',xsize=1000,title=title,bgcolor='black',
                badcolor='gray')
    hp.projplot(lon,lat,'r',lonlat=True,coord='G')
    hp.graticule(dpar=45.,dmer=30.,coord='C')
    pngfile1='frames/'+rootname+'_frame'+str(i)+'.png'
    plt.savefig(pngfile1,dpi=300)
    plt.close()

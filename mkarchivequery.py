#import psycopg2
#import petl as etl
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import time
import healpy as hp
import getopt,sys,os
from pylab import cm
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import astropy.utils as autils
import requests
import json
from pprint import pprint as pp  # pretty print

############
# Globals  #
############


natroot = 'https://astroarchive.noirlab.edu'
fieldlist = ["archive_filename","ra_center","dec_center","proposal","ifilter","MJD-OBS","OBJECT","survey","dateobs_center","instrument","proc_type","exposure","obs_type","prod_type"]
searchlist = [["instrument", "decam", "90prime", "mosaic", "mosaic3", "mosaic_1", "mosaic_1_1", "mosaic_2"], ["prod_type","image"],["proc_type","raw"],["obs_type","object"],["exposure",30.,36000.]]

adsurl = f'{natroot}/api/adv_search'

outfile='testmap.fits'
nside=256

def Usage():
    print("usage: python mkarchivequery.py -h -o out.fits -s searchlist -f fieldlist -n nside")

def dbquery(searchlist,fieldlist):
    jj = {"outfields" : fieldlist, "search" : searchlist}
    apiurl = f'{adsurl}/find/?limit=10000000'
    df = pd.DataFrame(requests.post(apiurl,json=jj).json()[1:])
    df['OBJECT']=df['OBJECT'].str.replace(","," ")
    df['ra']=df['ra_center'].astype('float')
    df['dec']=df['dec_center'].astype('float')
    df['exposure']=df['exposure'].astype('float')
    df=df[(np.isnan(df['ra']) | np.isnan(df['dec']))==False]
        
    return df

def mkmap(df,nside=256):
    #dra, ddec = np.loadtxt('DECam_profile.txt',unpack=True)
    mosra=np.array([-0.31,-0.31,0.31,0.31])
    mosdec=np.array([-0.31,0.31,0.31,-0.31])
    bokra=np.array([-0.56,-0.56,0.56,0.56])
    bokdec=np.array([-0.56,0.56,0.56,-0.56])
    #dra=np.append(dra,dra[0])
    #ddec=np.append(ddec,ddec[0])

    npix = hp.nside2npix(nside)
    map = np.zeros(npix)


    dfd = df[df['instrument']=='decam']
    dfm = df[df['instrument'].str.contains(r'mosaic')]
    dfb = df[df['instrument'].str.contains(r'90prime')]
    nind=int(1.1**2*np.pi/(hp.nside2resol(nside,arcmin=True)/60)**2*1.1)
    ptab=np.zeros(len(df),dtype=[('np','i8'),('ind','%di8'%nind)])

    radius=2.2/2*np.pi/180
    print(len(dfd),' DECam exposures')
    ptabd=np.zeros(len(dfd),dtype=[('np','i8'),('ind','%di8'%nind)])
    for idx in range(0,len(dfd)):
        row=dfd.iloc[idx]
        if row['ra'] is not None and row['dec'] is not None and row['exposure'] is not None:
            ra0=np.double(row['ra'])
            dec0=np.double(row['dec'])

            #ra1=ra0+dra/np.cos(dec0*np.pi/180)
            #dec1=dec0+ddec
            #dec1=np.clip(dec1,-90,90)

            vec=hp.ang2vec(ra0,dec0,lonlat=True)
            #ipix=hp.query_polygon(nside,vec)
            ipix=hp.query_disc(nside,vec,radius)
            map[ipix]=map[ipix]+float(row['exposure'])
            ptabd['np'][idx]=len(ipix)
            ptabd['ind'][idx][0:len(ipix)]=ipix

    print(len(dfm),' Mosaic exposures')
    ptabm=np.zeros(len(dfm),dtype=[('np','i8'),('ind','%di8'%nind)])
    for idx in range(0,len(dfm)):
        row=dfm.iloc[idx]
        if row['ra'] is not None and row['dec'] is not None and row['exposure'] is not None:
            ra0=np.double(row['ra'])
            dec0=np.double(row['dec'])
            ra1=ra0+mosra/np.cos(dec0*np.pi/180)
            dec1=dec0+mosdec
            dec1=np.clip(dec1,-90,90)

            vec=hp.ang2vec(ra1,dec1,lonlat=True)
            try:
                ipix=hp.query_polygon(nside,vec)
                map[ipix]=map[ipix]+float(row['exposure'])
                ptabm['np'][idx]=len(ipix)
                ptabm['ind'][idx][0:len(ipix)]=ipix
            except RuntimeError:
                print(idx,'Problem with Mosaic polygon')
                
    ptab[df['instrument']=='decam']=ptabd
    ptab[df['instrument'].str.contains(r'mosaic')]=ptabm

    print(len(dfb),' Bok exposures')
    ptabb=np.zeros(len(dfb),dtype=[('np','i8'),('ind','%di8'%nind)])
    for idx in range(0,len(dfb)):
        row=dfb.iloc[idx]
        if row['ra'] is not None and row['dec'] is not None and row['exposure'] is not None:
            ra0=np.double(row['ra'])
            dec0=np.double(row['dec'])
            ra1=ra0+bokra/np.cos(dec0*np.pi/180)
            dec1=dec0+bokdec
            dec1=np.clip(dec1,-90,90)

            vec=hp.ang2vec(ra1,dec1,lonlat=True)
            try:
                ipix=hp.query_polygon(nside,vec)
                map[ipix]=map[ipix]+float(row['exposure'])
                ptabb['np'][idx]=len(ipix)
                ptabb['ind'][idx][0:len(ipix)]=ipix
            except RuntimeError:
                print(idx,'Problem with Bok polygon')
                
    ptab[df['instrument'].str.contains(r'90prime')]=ptabb


    return map,ptab

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"ho:s:f:n:",["help","out=","searchlist=","fieldlist=","nside="])

    except getopt.GetoptError:
        print("Error processing arguments:")

    for opt,arg in opts:
        if opt in ("-h","--help"):
            Usage()
            sys.exit()
        elif opt in ("-o","--out"):
            outfile=arg
        elif opt in ("-s","--searchlist"):
            searchlist=arg
        elif opt in ("-s","--fieldlist"):
            fieldlist=arg
        elif opt in ("-n","--nside"):
            nside=np.int(arg)

    fieldlist = ["archive_filename","ra_center","dec_center","proposal","ifilter","MJD-OBS","OBJECT","survey","dateobs_center","instrument","proc_type","exposure","obs_type","prod_type"]
    searchlist = [["instrument", "decam", "90prime", "mosaic", "mosaic3", "mosaic_1", "mosaic_1_1", "mosaic_2"], ["prod_type","image"],["proc_type","raw"],["obs_type","object"],["exposure",30.,36000.]]

    df=dbquery(searchlist,fieldlist)
    (map,ptab)=mkmap(df,nside=nside)

    hp.write_map(outfile,map,coord='C',overwrite=True)

    if nside < 2048:
        indexfile=os.path.splitext(outfile)[0]+'_index.fits'
        t = Table(ptab)
        t.write(indexfile,overwrite=True)
        datafile=os.path.splitext(outfile)[0]+'_data.fits'
        csv=df.to_csv()
        t2=ascii.read(csv,format='csv')
        t2.write(datafile,overwrite=True)

    pngfile1=os.path.splitext(outfile)[0]+'_300.png'
    pngfile2=os.path.splitext(outfile)[0]+'_600.png'
    cmap = cm.jet
    cmap.set_under('w')
    lon=np.arange(360)
    lat=np.zeros(360)
    
    datestr=time.strftime('%B %d, %Y')
    hp.mollview(np.log10(map+0.1),cmap=cmap,min=0.0,max=6.,rot=-80,
                flip='geo',coord='C',cbar=True,notext=True,
                title=datestr,norm='linear',xsize=1000)
    hp.projplot(lon,lat,'r',lonlat=True,coord='G')
    hp.graticule(dpar=45.,dmer=30.,coord='C')
    plt.savefig(pngfile1,dpi=300)
    plt.savefig(pngfile2,dpi=600)
    #plt.savefig('/Users/olsen/Dropbox/ExposureMap/'+pngfile1,dpi=300)
    #plt.savefig('/Users/olsen/Dropbox/ExposureMap/'+pngfile2,dpi=600)

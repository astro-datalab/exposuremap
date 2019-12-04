import psycopg2
import petl as etl
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import time
import healpy as hp
import getopt,sys,os
from pylab import cm
import matplotlib.pyplot as plt

############
# Globals  #
############

query="""SELECT reference,fits_extension,object,survey,surveyid,
prop_id,start_date,ra,dec,equinox,telescope,exposure,depth,depth_err,
magzero,magerr,seeing,release_date,airmass,date_obs,dtutc,dtpi,dtpropid,
filter,ha,instrument,mjd_obs,observat,obstype,proctype,prodtype,corn1dec,
corn2dec,corn3dec,corn4dec,corn1ra,corn2ra,corn3ra,corn4ra FROM voi.siap 
WHERE (exposure >= 30 AND obstype = 'object') AND (((telescope = 'ct4m' 
AND instrument = 'decam')) OR ((telescope = 'ct4m' 
AND instrument = 'mosaic_2')) OR ((telescope = 'kp4m' 
AND instrument LIKE 'mosaic%')) OR ((telescope = 'bok23m' AND instrument = '90prime'))) AND proctype = 'InstCal' AND prodtype = 'image' ORDER BY dtutc ASC"""
outfile='testmap.fits'
nside=256

def Usage():
    print("usage: python mkquery.py -h -o out.fits -q query_string -n nside")

def dbquery(query):

    try:
        conn = psycopg2.connect("dbname='metadata' user='dbreader' host='db.sdm.noao.edu' password='Mosaic_DHS' port=5432")
    except:
        print("I am unable to connect to the database")

    table = etl.fromdb(conn, query)
    df=etl.todataframe(table)
    df['object']=df['object'].str.replace(","," ")
    df['ra']=df['ra'].astype('float')
    df['dec']=df['dec'].astype('float')
    df['equinox']=df['equinox'].astype('float')
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
            map[ipix]=map[ipix]+np.float(row['exposure'])
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
                map[ipix]=map[ipix]+np.float(row['exposure'])
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
                map[ipix]=map[ipix]+np.float(row['exposure'])
                ptabb['np'][idx]=len(ipix)
                ptabb['ind'][idx][0:len(ipix)]=ipix
            except RuntimeError:
                print(idx,'Problem with Bok polygon')
                
    ptab[df['instrument'].str.contains(r'90prime')]=ptabb


    return map,ptab

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"ho:q:n:",["help","out=","query=","nside="])

    except getopt.GetoptError:
        print("Error processing arguments:")

    for opt,arg in opts:
        if opt in ("-h","--help"):
            Usage()
            sys.exit()
        elif opt in ("-o","--out"):
            outfile=arg
        elif opt in ("-q","--query"):
            query=arg
        elif opt in ("-n","--nside"):
            nside=np.int(arg)

    #query="""SELECT reference,fits_extension,object,survey,surveyid,prop_id,start_date,ra,dec,equinox,telescope,exposure,depth,depth_err,magzero,magerr,seeing,release_date,airmass,date_obs,dtutc,dtpi,dtpropid,filter,ha,instrument,mjd_obs,observat,obstype,proctype,prodtype,corn1dec,corn2dec,corn3dec,corn4dec,corn1ra,corn2ra,corn3ra,corn4ra FROM voi.siap WHERE (exposure >= 30 AND obstype = 'object') AND (((telescope = 'ct4m' AND instrument = 'decam')) OR ((telescope = 'ct4m' AND instrument = 'mosaic_2')) OR ((telescope = 'kp4m' AND instrument LIKE 'mosaic%'))) AND proctype = 'InstCal' ORDER BY dtutc ASC"""
    
    query="""SELECT reference,fits_extension,object,survey,surveyid,
prop_id,start_date,ra,dec,equinox,telescope,exposure,depth,depth_err,
magzero,magerr,seeing,release_date,airmass,date_obs,dtutc,dtpi,dtpropid,
filter,ha,instrument,mjd_obs,observat,obstype,proctype,prodtype,corn1dec,
corn2dec,corn3dec,corn4dec,corn1ra,corn2ra,corn3ra,corn4ra FROM voi.siap 
WHERE (exposure >= 30 AND obstype = 'object') AND (((telescope = 'ct4m' 
AND instrument = 'decam')) OR ((telescope = 'ct4m' 
AND instrument = 'mosaic_2')) OR ((telescope = 'kp4m' 
AND instrument LIKE 'mosaic%')) OR ((telescope = 'bok23m' AND instrument = '90prime'))) AND proctype = 'InstCal' AND prodtype = 'image' ORDER BY dtutc ASC"""

    df=dbquery(query)
    df=df[df['obstype']=='object']
    (map,ptab)=mkmap(df,nside=nside)

    hp.write_map(outfile,map,coord='C')
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

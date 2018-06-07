# plume code
#
import numpy as np
import datetime
import time
from pydap.client import open_url
import scipy.io
import pdb
from datetime import date
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,Point
import scipy.interpolate
from mpl_toolkits.basemap import Basemap
import pyproj
import math
from lonlat2km import lonlat2km


#res=input('Enter the resolution (2 or 6 km): ')
res=2
# 2
#days=input('Enter the number of the days to capture: ')
days=3
year,month,day,hour,minu,secs=input("Enter the start date [yr mn dy hr mn sc] ").split()
#[year,month,day,hour,minu,secs]=input('Enter the start date [yr,mn,dy,hr,mn,sc]): ')
#[2017,5,27,5,0,0]
if int(res)==2:
    url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/2km/hourly/RTV/HFRADAR_US_West_Coast_2km_Resolution_Hourly_RTV_best.ncd'
if int(res)==6:
    url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/6km/hourly/RTV/HFRADAR_US_West_Coast_6km_Resolution_Hourly_RTV_best.ncd'
maxlat,minlat=input('Enter the north and south boundary []: ').split()
maxlat=float(maxlat)
minlat=float(minlat)
# [38.5,37]
minlon,maxlon=input('Enter the west and east boundary []: ').split()
minlon=float(minlon)
maxlon=float(maxlon)
# [-123.5, -122]
#
# need to keep the track of the fractional part.
offset=datetime.datetime(1970,1,1).toordinal()
time1=datetime.datetime(int(year),int(month),int(day),int(hour),int(minu),int(secs)).toordinal()
time2=time1+int(days)
avgtime=(time1+time2)/2
#
# open the url and get lazy load value to the data
#
datapointer=open_url(url)
timevec=datapointer['time']
longitude=datapointer['lon']
latitude=datapointer['lat']
timeoff=datapointer['time_offset']
time=timevec[:]
latitude=latitude[:]
longitude=longitude[:]
#
indlat=[i for i,e in enumerate(latitude) if ((e >= minlat)and(e <= maxlat))]
indlon=[i for i,e in enumerate(longitude) if ((e >= minlon)and(e <= maxlon))]
#
time_off=timevec.units
year=int(time_off[12:16])
month=int(time_off[17:19])
day=int(time_off[20:22])
hour=int(time_off[23:25])
minu=int(time_off[26:28])
secs=int(time_off[29:31])
newoff=datetime.datetime(year,month,day,hour,minu,secs).toordinal()
time=time/24+newoff
indtime=[i for i,e in enumerate(time) if (e > float(time1)and(e <= float(time2)))]
time=time[indtime]
data1=datapointer['u']
newu=data1['u'][indtime[0]:indtime[-1]+1,indlat[0]:indlat[-1]+1,indlon[0]:indlon[-1]+1]
data2=datapointer['v']
newv=data2['v'][indtime[0]:indtime[-1]+1,indlat[0]:indlat[-1]+1,indlon[0]:indlon[-1]+1]
#
newu=newu*100
newv=newv*100
#
longitude=longitude[indlon]
latitude=latitude[indlat]
#
x,y=np.meshgrid(longitude,latitude)
#x=np.transpose(x)
#y=np.transpose(y)
newu=np.transpose(newu,(0,2,1))
newv=np.transpose(newv,(0,2,1))
#
coastline=scipy.io.loadmat('c:\plume_software\SantaClaraPlumeTraj_example\simple_coast.mat')
xb=coastline['xb']
yb=coastline['yb']
xb=np.append(xb,(-121.4955,-124.0604))
yb=np.append(yb,(39.9989,39.9989))
#
# create projection information using a mercator projection for the area of interest.
# will help when creating polygon of the coast to mask data over land
wgs84=pyproj.Proj(init='EPSG:4326')
npmerc=pyproj.Proj("+proj=merc +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# projection of the polygon of the coastline
polyx,polyy=pyproj.transform(wgs84,npmerc,xb,yb)
#
gridx,gridy=pyproj.transform(wgs84,npmerc,x,y)
# actual polygon creation
poly_proj=Polygon(zip(polyx,polyy))
# create blank mask array
mask=np.zeros_like(x)
dist=np.zeros_like(x)
(lx,ly)=gridx.shape
# blank land in mask (value=1)
for i in range(lx):
    for j in range(ly):
        grid_point=Point(gridx[i][j],gridy[i][j])
        dist[i][j]=grid_point.distance(poly_proj)
        if grid_point.within(poly_proj):
            mask[i][j]=1
#
#test=[]
# now sure how we automate this piece took me trial and error to get a good number for upper bound.
ll=(dist > 0)&(dist < 2500)
zboundx=x[ll]
zboundy=y[ll]
#
# reorder values so that we can actually have a coast line boundary and not a mess 
lz=len(zboundx)
testi=[]
for i in range(lz):
    mdis=np.sqrt((xb-zboundx[i])**2+(yb-zboundy[i])**2)
    dmin=np.min(mdis)
    ind=[j for j,e in enumerate(mdis) if e==dmin]
    testi=np.append(testi,ind)
zbound=np.argsort(testi)
zboundx=zboundx[zbound]
zboundy=zboundy[zbound]
indx=[]
indy=[]
for i in range(lz):
    ix=[j for j,e in enumerate(longitude) if e==zboundx[i]]
    iy=[j for j,e in enumerate(latitude) if e==zboundy[i]]
    indx=np.append(indx,ix)
    indy=np.append(indy,iy)
z=complex(0,1)
blist=[]
for i in range(lz):
    ca,cb=lonlat2km(zboundx[i],zboundy[i],xb,yb)
    cc=np.abs(ca+cb*z)
    cmin=np.min(cc)
    cind=[j for j,e in enumerate(cc) if e==cmin]
    blist=np.append(blist,int(cind[0]))
#
blist=blist.astype(int)
#
# these are the nearest grid points to the coast...
# u,v array shape is time,lon,lat
# we need to create a full grid for each time-step.  How best to do that?
#
gt,gx,gy=newu.shape
mask=np.transpose(mask)
land=mask==1
x=np.transpose(x)
y=np.transpose(y)
for k in range(gt):
    tmpx=x[:][:]
    tmpy=y[:][:]
    tmpu=np.squeeze(newu[k][:][:])
    tmpv=np.squeeze(newv[k][:][:])
    lu=np.isnan(tmpu)==0
    tmpx=tmpx[lu]
    tmpy=tmpy[lu]
    tmpu=tmpu[lu]
    tmpv=tmpv[lu]
    gridu=scipy.interpolate.griddata((tmpx,tmpy),tmpu,(x,y))
    gridv=scipy.interpolate.griddata((tmpx,tmpy),tmpv,(x,y))
    gridu[land]=np.nan
    gridv[land]=np.nan
    newu[k][:][:]=gridu
    newv[k][:][:]=gridv
#
# compute along coast distance and angles
#
a,b=lonlat2km(xb[0:-4],yb[0:-4],xb[1:-3],yb[1:-3])
c=np.abs(a+b*z) # magnitude
baxis=np.cumsum(c)
#aa,bb=lonlat2km(xb[0:-5],yb[0:-5],xb[2:-3],yb[2:-3])
#thb=np.arctan2(bb,aa) # angle in radians
thb=np.arctan2(b,a) # angle in radians
ang=thb*180/np.pi
#
#
dcf=0.3
cv=1e-5*3600
npts=50
nlife=24*3
xc=-122.5726
yc=37.7836
#
#
lts=len(time)
indx=indx.astype(int)
indy=indy.astype(int)
for k in range(lts):
    tmpu=np.squeeze(newu[k,indx,indy])
    tmpv=np.squeeze(newv[k,indx,indy])
    # need to shrink tmpu and tmpv to just the boundary points?
    # this needs to be added below
    
    mg=np.abs(tmpu+tmpv*z)*dcf
    ph=np.arctan2(tmpv,tmpu)
    phb=thb[blist]
    # this is not the correct shape?
    mgb=mg*np.cos(ph-phb)
    ub=mgb*np.cos(phb)
    vb=mgb*np.sin(phb)
    
    
#
# compute boundary around a set of points....
#
#from scipy.spactial import ConvexHull
#allPoints=np.colum_stack((lats,longs))
#hullPoints=ConvexHull(allPoints)
pdb.set_trace()
#    fu=scipy.interpolate.interp2d(tmpx,tmpy,tmpu)
#    fv=scipy.interpolate.interp2d(tmpx,tmpy,tmpv)
#    unew=fu(longitude,latitude)
#    vnew=fv(longitude,latitude)
    
# These remove SF bay and tomales bay
#xr=np.arange(263,309)
#xr=np.append(xr,np.arange(411,676))
#xb=np.delete(xb,xr)
#yb=np.delete(yb,xr)
#
#indy=[i for i,e in enumerate(yb) if e > 40]
#xb=np.delete(xb,indy)
#yb=np.delete(yb,indy)
#indx=[i for i,e in enumerate(xb) if e > -122]
#xb=np.delete(xb,indx)
#yb=np.delete(yb,indx)


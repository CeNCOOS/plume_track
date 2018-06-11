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
from shapely.geometry import Polygon,Point,LineString
import scipy.interpolate
from mpl_toolkits.basemap import Basemap
import pyproj
import math
from lonlat2km import lonlat2km,km2lonlat


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
uerr=2
#
#
lts=len(time)
indx=indx.astype(int)
indy=indy.astype(int)
gxx=x
gyy=y
x_=[]
y_=[]
xb=xb[0:-4]
yb=yb[0:-4]
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
    emg=np.interp(baxis,baxis[blist],mgb)
    eub=emg*np.cos(thb)
    evb=emg*np.sin(thb)
    u_=np.squeeze(newu[k,:,:])
    v_=np.squeeze(newv[k,:,:])
    #
    ival=np.isnan(u_)==0
    jval=np.isnan(eub)==0
    #
    # do we only do this for k==1?
    # reseeds for each time step otherwise?
    if k==0:
        x_=np.append(x_,xc*np.ones((npts,1)))
        y_=np.append(y_,yc*np.ones((npts,1)))
    #
    # remove particles beyond lifetime of particle
    #
    if k > nlife:
        x_=np.delete(x_,range(1,npts+1))
        y_=np.delete(y_,range(1,npts+1))
    #
    gtmpx=gxx[ival]
    gtmpy=gyy[ival]
    utmp=u_[ival]
    vtmp=v_[ival]
    testx=np.append(gtmpx, xb[jval])
    testy=np.append(gtmpy, yb[jval])
    testu=np.append(utmp, eub[jval])
    testv=np.append(vtmp, evb[jval])
    u_=scipy.interpolate.griddata((testx,testy),testu,(x_,y_))
    v_=scipy.interpolate.griddata((testx,testy),testv,(x_,y_))
    ntotalpt=len(x_)
    # gerate random directions
    th=2*np.pi*np.random.sample((ntotalpt,))
    # add random error to u and v uerr=2
    un_=u_+uerr*np.cos(th)
    vn_=v_+uerr*np.sin(th)
    mgn=np.abs(un_+z*vn_)
    phn=np.arctan2(vn_,un_)
    # cv is some time-step function?
    dx=un_*cv
    dy=vn_*cv
    xn_,yn_=km2lonlat(x_,y_,dx,dy)
    jj=1
    for j in range(ntotalpt):
        # create shapely object of boundary points as a line
        a1=LineString(zip(xb,yb))
        xx=[x_[j],xn_[j]]
        yy=[y_[j],yn_[j]]
        # create shapely object of 2 velocity starting points
        a2=LineString(zip(xx,yy))
        # find out if they intersect?
        cc=a2.intersection(a1)
        try:
            cx=cc[0]
            cy=cc[1]
            a,b=lonlat2km(cx,cy,xb,yb)
            c=np.abs(a+z*b)
            #pdb.set_trace()
            ii=[i for i,e in enumerate(c) if e==np.min(c)]              
            #ii=c==np.min(c)
            ii=ii[0]
            dx=eub[ii]*cv
            dy=evb[ii]*cv
            xn_[j],yn_[j]=km2lonlat(x_[j],y_[j],dx,dy)
        except:
            continue
    if k==0:
        xf=np.vstack(xn_)
        yf=np.vstack(yn_)
    else:
        xtest=np.expand_dims(xn_,axis=1)
        ytest=np.expand_dims(yn_,axis=1)
        xf=np.hstack((xf,xtest))
        yf=np.hstack((yf,ytest))
    x_=xn_
    y_=yn_
    
#
# compute boundary around a set of points....
#
myc=plt.cm.rainbow(np.linspace(0,1,72))
mp=Basemap(llcrnrlon=360-123.5,llcrnrlat=37,urcrnrlon=360-122,urcrnrlat=38.5,resolution='f',projection='cyl')
mp.fillcontinents()
mp.drawcoastlines()
mp.quiver(360+testx,testy,testu,testv)
for i in range(72):
    plt.plot(xf[:,i]+360,yf[:,i],'.',color=myc[i])
llon=np.arange(360-123.5,360-122,0.25)
llat=np.arange(37,38.5,0.25)
mp.drawmeridians(llon,labels=[1,0,0,1])
mp.drawparallels(llat,labels=[1,0,0,1])
plt.show()


# plume code
#
import numpy as np
import datetime, sys
from pydap.client import open_url
import scipy.io
import pdb
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon,Point
import scipy.interpolate
import xarray as xr
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import pyproj
#from lonlat2km import lonlat2km

class TrackingModel():
    
    def __init__(self):
        self.minlon, self.maxlon, self.minlat, self.maxlat, self.resolution_km, self.url, self.days_to_capture, self.start_time = self.set_model_parameters(default=True)
        self.current_dataset = self._get_HFR_subset()
        self.lat = self.current_dataset['lat'].values
        self.lon = self.current_dataset['lon'].values
        self.x = np.arange(0, len(self.lon), 1) * 1850
        self.y = np.arange(0, len(self.lat), 1) * 1995
        self.x_grid, self.y_grid = np.meshgrid(self.x, self.y)
        self.lon_grid, self.lat_grid = np.meshgrid(self.lon, self.lat)
        self.origin = ccrs.UTM(zone=10).transform_point(self.lon[0], self.lat[0], ccrs.PlateCarree())
        self.time_index = 0
        self.current_time = self.start_time
        self.particles = np.array([])
        
    
    def set_model_parameters(self, default=False):
        ''' Hotwire the default parameters for developement '''
        crs = ccrs.PlateCarree() # Make cartopy projection object
        
        resolution_km = 2
        days_to_capture=3
    
        if int(resolution_km)==2:
            url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/2km/hourly/RTV/HFRADAR_US_West_Coast_2km_Resolution_Hourly_RTV_best.ncd'
        if int(resolution_km)==6:
            url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/6km/hourly/RTV/HFRADAR_US_West_Coast_6km_Resolution_Hourly_RTV_best.ncd'
    
        if default:
            year,month,day,hour,minu,secs = 2017,5,27,5,0,0
            start_time = datetime.datetime(year,month,day,hour,minu,secs)
            maxlat, minlat = 38.5, 37
            maxlat=float(maxlat)
            minlat=float(minlat)
            minlon, maxlon = -123.5, -122
            minlon=float(minlon)
            maxlon=float(maxlon)
    
        else:
            try:
                year,month,day,hour,minu,secs=input("Enter the start date [yr mn dy hr mn sc] ").split()
            except ValueError:
                print('Defaulting to 2017-05-27T05:00:00')
                year,month,day,hour,minu,secs = 2017,5,27,5,0,0
                start_time = datetime.datetime(year,month,day,hour,minu,secs)
            #[year,month,day,hour,minu,secs]=input('Enter the start date [yr,mn,dy,hr,mn,sc]): ')
            #[2017,5,27,5,0,0]
            try:
                maxlat,minlat=input('Enter the north and south boundary []: ').split()
            except ValueError:
                maxlat, minlat = 38.5, 37
            maxlat=float(maxlat)
            minlat=float(minlat)
            # [38.5,37]
            try:
                minlon, maxlon=input('Enter the west and east boundary []: ').split()
            except ValueError:
                minlon, maxlon = -123.5, -122
            minlon=float(minlon)
            maxlon=float(maxlon)
    
        return minlon, maxlon, minlat, maxlat, resolution_km, url, days_to_capture, start_time

    def _get_HFR_subset(self):
    
        """
        Open netcdf from hfrnet as xarray datase, subset the desired data in space and time
        Currently, using global variables for parameters, should reduce as code becomes object-oriented
    
        :return hfr_ds = xarray dataset over the specified spatial grid and time
        """
    
    
        try:
            ds = xr.open_dataset(self.url)
        except Exception as e:
            print("Trouble downloading data, check connections.")
            print(e)
            sys.exit()
    
        return ds.sel(time = slice(self.start_time, self.start_time + datetime.timedelta(days=self.days_to_capture)),
                      lat =  slice(self.minlat,self.maxlat),
                      lon = slice(self.minlon,self.maxlon))


    def add_particle(self, coord):
        utm = ccrs.UTM(zone=10)
        p_dist = utm.transform_point(coord[0], coord[1], ccrs.PlateCarree()) 
        p_dist = [p_dist[0] - self.origin[0], p_dist[1] - self.origin[1]]
        self.particles = np.append(self.particles, Particle(p_dist))
    
    def part_to_coor(self, particle):
        pc_proj = ccrs.PlateCarree()
        p_dist = particle.get_position() 
        p_dist = [p_dist[0] + self.origin[0], p_dist[1] + self.origin[1]]
        p_coor = pc_proj.transform_point(p_dist[0], p_dist[1], ccrs.UTM(zone=10))
        return p_coor

    def draw_map(self, draw_bathy=True):
        ''' Draw a map of the domain '''
    
        land_feature  = cfeature.NaturalEarthFeature('physical','land','50m')
        self.fig = plt.figure()
        self.fig.set_size_inches(8,8)
        self.geo_axes = plt.axes(projection=ccrs.PlateCarree())
        extent = [self.minlon, self.maxlon, self.minlat, self.maxlat]
        self.geo_axes.set_extent(extent, crs=ccrs.PlateCarree())
        self.geo_axes.add_feature(land_feature, edgecolor='k', zorder=40)
        self.geo_axes.gridlines(draw_labels=True, zorder= 20)
        if draw_bathy:
            self._draw_bathymetry_SFBAY()


    def _draw_bathymetry_SFBAY(self):
        ''' 
        Draw and return an axes with bathymetry contours
        
        Bathy grid retrieved from: https://maps.ngdc.noaa.gov/viewers/wcs-client/
        '''
        file = './data/SF_Bay_coastal_relief_model.nc'
        bathy = xr.open_dataset(file)
        blat = bathy['lat'].values
        blon = bathy['lon'].values
        elev = bathy['Band1'].values
        levels = np.arange(10*(-2000//10), -499, 450)
        levels = np.append(levels, np.arange(-475, 25, 25))
        lws = [0.5 if l%100 else 1 for l in levels]
        cs = self.geo_axes.contour(blon, blat, elev, levels, linewidths=lws, linestyles='solid', colors=['black'], alpha=0.4)
        plt.clabel(cs, list(np.arange(-125,1,25)),fmt='%1d', inline=True, fontsize=15, colors='k',inline_spacing=10)        

    def plot_particles(self):
        for p in self.particles[:]:
            pos = self.part_to_coor(p)
            self.geo_axes.scatter(pos[0],pos[1], zorder=50)
    

    def advect_particle(self, particle, dim, h=.25):
        ''' 
        Advect a particle object based on the current
        
        '''
        current_pos = particle.get_position()
        new_x = self._solve_position(current_pos[0], current_pos[1], h, dim=0)
        new_y = self._solve_position(current_pos[0], current_pos[1], h, dim=1)
        particle.update_position([new_x, new_y])

    def _solve_position(self, x, y, h, dim):
            """
            Solves for the next position of a particle in after time, h, in either the x
            or y using a runge-kutta 4th order scheme.
        
            TODO: Update function to get next half timestep if goes into next hour
        
            Arguments
            ---------
            X, Y: mesh grid.
            x, y: coordinates where to begin the evolution.
            f: the current vector f that will be evolved.
            h: the time step
            dim: 0 for x and 1 for y.
        
            Returns
            ---------
            interp_value: interpolated value of f(y,x)
            """
            X = self.x_grid
            Y = self.y_grid
            if dim == 0:
                f = self.current_dataset['u'].isel(time=self.time_index).values * 24* 60 * 60
            if dim == 1:
                f = self.current_dataset['v'].isel(time=self.time_index).values * 24* 60 * 60
                
            k1 = h * self.bilinear_interpolation(X, Y, f, x, y)
            k2 = h * self.bilinear_interpolation(X, Y, f, x + 0.5 * h, y + 0.5 * k1)
            k3 = h * self.bilinear_interpolation(X, Y, f, x + 0.5 * h, y + 0.5 * k2)
            k4 = h * self.bilinear_interpolation(X, Y, f, x + h, y + k3)
            try:
                if dim == 0:
                    return x + 1. / 6 * k1 + 1. / 3 * k2 + 1. / 3 * k3 + 1. / 6 * k4
                elif dim == 1:
                    return y + 1. / 6 * k1 + 1. / 3 * k2 + 1. / 3 * k3 + 1. / 6 * k4
            except Exception as e:
                print(e.with_traceback())
                sys.exit()

    def bilinear_interpolation(self, X, Y, f, x, y):
        """
        Interpolation methods for estimating surface current values in between grid points. Edge cases are outlined in the
        top of the function and may need to be refactored. NaNs is returned for cases where values can not be interpolated.
    
        Arguments
        ---------
        X, Y: Coordinate mesh grid
        f: Grid of velocity values that can be accessed as f(j,i) Remeber row, column
        x, y: coordinates to compute interpolation to f(y,x)
    
        Returns
        ---------
        interp_value: interpolated value of f(y,x)
    
        """
        
        # Grid index shape
        M = np.shape(X[:, 0])[0]
        N = np.shape(X[0, :])[0]
    
        dx, dy = X[0, 1] - X[0, 0], Y[1, 0] - Y[0, 0]
        x_start, y_start = X[0, 0], Y[0, 0]
    
        # Find the index of each value
        i1, i2 = int((x - x_start) / dx), int((x - x_start) / dx) + 1
        j1, j2 = int((y - y_start) / dy), int((y - y_start) / dy) + 1
    
        # Boundary Conditions when interpolating near the edge.
        # 1. Eastern boundary
        if (i1 - N) > 1:
            return np.nan
        if i1 >= N - 1 and j1 <= N - 1 and j1 >= 0: # If on the Eastern edge of the boundary
            return f[j1, N - 1]
        if i1 >= N - 1 and j1 <= 0:
            return f[0, N - 1]
        if i1 >= N - 1 and j1 >= N - 1:
            return f[N - 1, N - 1]
    
        # 2. Western boundary
        if i1 <= 0 and j1 <= N - 1 and j1 >= 0:
            return f[j1, 0]
        if i1 <= 0 and j1 <= 0:
            return f[0, 0]
        if i1 <= 0 and j1 >= N - 1:
            return f[N - 1, 0]
    
        # 3. Northern boundary
        if j1 >= M - 1 and i1 <= M - 1 and i1 >= 0:
            return f[M - 1, i1]
        if j1 >= N - 1 and i1 <= 0:
            return f[M - 1, 0]
    
        # 3. Bottom boundary
        if j1 <= 0 and i1 <= M - 1 and i1 >= 0:
            return f[0, i1]
        if j1 <= 0 and i1 >= M - 1:
            return f[M - 1, 0]
    
        x1, x2 = X[j1, i1], X[j2, i2]
        y1, y2 = Y[j1, i1], Y[j2, i2]
    
        interp_value = (1 / (x2 - x1) * 1 / (y2 - y1) *
                          (f[j1, i1] * (x2 - x) * (y2 - y) + f[j1, i2] * (x - x1) * (y2 - y)
                           + f[j2, i1] * (x2 - x) * (y - y1) + f[j2, i2] * (x - x1) * (y - y1)))
    
        return interp_value




class Particle():
    
    def __init__(self,coord):
        self.coordinates = np.array([coord])
        
    def get_position(self):
        ''' Return latest postion of particle class '''
        return self.coordinates[-1,:]
    
    def update_position(self, pos):
        ''' Append to coord array '''
        pos = np.array(pos)
        self.coordinates = np.vstack((self.coordinates,pos))

def run_old():

    # need to keep the track of the fractional part.
    offset=datetime.datetime(1970,1,1).toordinal()
    time1=datetime.datetime(int(year),int(month),int(day),int(hour),int(minu),int(secs)).toordinal()
    time2=time1+int(days_to_capture)
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
    coastline=scipy.io.loadmat('.')
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



if __name__ == "__main__":
    model = TrackingModel()
    model.draw_map()    
    # Debugging Particles
    model.add_particle([-122.6, 37.73])
    ax = model.geo_axes
    
    ax.quiver(model.lon_grid, model.lat_grid, model.current_dataset['u'].isel(time=0), model.current_dataset['v'].isel(time=0))
    ix = np.where(np.isfinite(model.current_dataset['u'].isel(time=0)))
    ax.scatter(model.lon_grid[ix], model.lat_grid[ix])
    
    for i in range(1):
        model.advect_particle(model.particles[0], h=.25, dim=0)
        model.plot_particles()

#    



    # MASK HACk
#    mask = regionmask.defined_regions.natural_earth.land_110.mask(model.current_dataset, wrap_lon=False) # THis Mask is Onland (110 m resolution)
#    u = t['u'].where(mask != 0)
#    u  = u.isel(time=0)
#    model.geo_axes.pcolormesh(model.lon_grid, model.lat_grid, u)
    
    
    
    
#    ax.quiver(xx,yy,start_ds['u'],start_ds['v'])
#    

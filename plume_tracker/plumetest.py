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
import hfr_util
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
        self.hours_elapsed = 0
        self.current_time = self.start_time
        self.particles = np.array([])
        self.time_step = .25


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

        subset_ds = ds.sel(time = slice(self.start_time, self.start_time + datetime.timedelta(days=self.days_to_capture)),
                      lat =  slice(self.minlat,self.maxlat),
                      lon = slice(self.minlon,self.maxlon))
        clean_u = hfr_util.interp_time_surface_currents(subset_ds['u'].values) # linear interpolate currents through time, given a threshold of availible data
        clean_v = hfr_util.interp_time_surface_currents(subset_ds['v'].values)
        subset_ds['u_clean'] =  xr.DataArray(clean_u, coords={'time': subset_ds['time'].values, 'lon': subset_ds['lon'].values, 'lat': subset_ds['lat'].values}, dims=['time','lat','lon'])
        subset_ds['v_clean'] =  xr.DataArray(clean_v, coords={'time': subset_ds['time'].values, 'lon': subset_ds['lon'].values, 'lat': subset_ds['lat'].values}, dims=['time','lat','lon'])
        return subset_ds

    def add_particle(self, coord):
        utm = ccrs.UTM(zone=10)
        p_dist = utm.transform_point(coord[0], coord[1], ccrs.PlateCarree())
        p_dist = [p_dist[0] - self.origin[0], p_dist[1] - self.origin[1]]
        self.particles = np.append(self.particles, Particle(p_dist))

    def part_to_coor(self, particle, last=False):
        pc_proj = ccrs.PlateCarree()
        p_dist = particle.coordinates
        p_dist = [p_dist[:,0] + self.origin[0], p_dist[:,1] + self.origin[1]]
        if last:
            p_coor = pc_proj.transform_points(ccrs.UTM(zone=10), np.array((p_dist[:,-1])),np.array((p_dist[:,-1])))
        else:
            p_coor = pc_proj.transform_points(ccrs.UTM(zone=10), np.array([p_dist[:,0]]),np.array([p_dist[:,1]]))
        p_coor = p_coor[:,:2]
        return p_coor

    def draw_map(self, draw_bathy=True):
        ''' Draw a map of the domain '''

        land_feature  = cfeature.NaturalEarthFeature('physical','land','50m')
        self.fig = plt.figure()
#        self.fig.set_size_inches(8,8)
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

    def plot_particles(self, last):
        if last:
            #Only if you want to plot the last position (debugging)
            for p in self.particles[:]:
                pos = self.part_to_coor(p)
                self.geo_axes.scatter(pos[0,0],pos[0,1], zorder=50, marker='.', c='grey')
        else:
            for p in self.particles[:]:
                self.geo_axes.plot(p.coordinates[:,0]+self.origin[0], p.coordinates[:,1]+self.origin[1],
                                   zorder=50,
                                   c='blue',
                                   markerfacecolor='.5',
                                   markeredgecolor='None',
                                   marker='.',
                                   markersize=5,
                                   transform=ccrs.UTM(zone=10))

    def advect_particle(self):
        '''
        Advect a particle object based on the current
        
        
        Fixed timestep to 1/4 hour --> If this is to be worked update self.timestep
        '''
        for p in self.particles:
            current_pos = p.get_position()
            if not np.any(np.isnan(current_pos)): # check if advection previously failed if so skip
                new_x = self._solve_position(current_pos[0], current_pos[1], self.time_step, dim=0)
                new_y = self._solve_position(current_pos[0], current_pos[1], self.time_step, dim=1)
                p.update_position([new_x, new_y])
        self.update_time()
    
    def update_time(self):
        ''' Update the time_index and number of elapsed hours with each advection '''
        self.hours_elapsed += self.time_step
        self.time_index = int(self.hours_elapsed)

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
                f = self.current_dataset['u_clean'].isel(time=self.time_index).values * 60 * 60
            if dim == 1:
                f = self.current_dataset['v_clean'].isel(time=self.time_index).values * 60 * 60
            try:
                k1 = h * self.bilinear_interpolation(X, Y, f, x, y)
                k2 = h * self.bilinear_interpolation(X, Y, f, x + 0.5 * h, y + 0.5 * k1)
                k3 = h * self.bilinear_interpolation(X, Y, f, x + 0.5 * h, y + 0.5 * k2)
                k4 = h * self.bilinear_interpolation(X, Y, f, x + h, y + k3)
            except ValueError as e:
                print('Error in Interpolation, trying to interpolate a NAN value')
                print(e)
                return
            
            try:
                if dim == 0:
                    return x + 1. / 6 * k1 + 1. / 3 * k2 + 1. / 3 * k3 + 1. / 6 * k4
                elif dim == 1:
                    return y + 1. / 6 * k1 + 1. / 3 * k2 + 1. / 3 * k3 + 1. / 6 * k4
            except Exception as e:
                print(e)
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

    def seed_particles(self, center_coord, radius):
        ''' Create a of cluster of particles within a radius in (km)'''

        x_pos, y_pos = ccrs.UTM(zone=10).transform_point(center_coord[0],center_coord[1], ccrs.PlateCarree())
        
        col_num = [1,3,5,3,1]
        dx = radius/4 * 1000
        ylevel = (np.arange(max(col_num)) - 3) * dx
        for i, n in enumerate(col_num):
            x = np.arange(n)
            x = x - (n-1)/2
            x = x * dx
            y = np.ones(shape=x.shape) * ylevel[i]
            pos_x = x_pos + x
            pos_y = y_pos + y
            coors = ccrs.PlateCarree().transform_points(ccrs.UTM(zone=10), pos_x, pos_y)
            x_coors = coors[:,0]
            y_coors = coors[:,1]
            for pos in zip(x_coors, y_coors):
                self.add_particle(pos)
        


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



if __name__ == "__main__":
    model = TrackingModel()
    model.draw_map(draw_bathy=True)
    # Debugging Particles
    center_coord = [-122.6, 37.77]
    model.seed_particles(center_coord, radius=3)  
    ax = model.geo_axes
#    ax.quiver(model.lon_grid, model.lat_grid, model.current_dataset['u_clean'].isel(time=0), model.current_dataset['v'].isel(time=0))
    ix = np.where(np.isfinite(model.current_dataset['u_clean'].isel(time=0)))
    ixnan = np.where(np.isnan(model.current_dataset['u_clean'].isel(time=0)))
#    ax.scatter(model.lon_grid[ix], model.lat_grid[ix], marker='.', s=10)
#    ax.scatter(model.lon_grid[ixnan], model.lat_grid[ixnan], marker='x', s=10, c='r')
#
    for i in range(24*4*2):
        try:
            model.advect_particle()
            
        except Exception as e:
            print(e)
            print(round(i/4,2),'hours have passed before breaking')
            break
    model.plot_particles(last=False)


    # MASK HACk
#    mask = regionmask.defined_regions.natural_earth.land_110.mask(model.current_dataset, wrap_lon=False) # THis Mask is Onland (110 m resolution)
#    u = t['u'].where(mask != 0)
#    u  = u.isel(time=0)
#    model.geo_axes.pcolormesh(model.lon_grid, model.lat_grid, u)

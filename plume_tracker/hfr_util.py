# -*- coding: utf-8 -*-
import numpy as np
import pdb

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.
    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def interp_time_surface_currents(current_array_with_nans):
    """ Linear interpolate gaps in timeseries at each grid point for the event. Skips gridpoints where over 90% of the
        data are NANs (these are generally outside of the domain of HFR, but remain in the array to keep the grid rectangular)
    Input:
        - current_array_with_nans, 3d numpy array with possible NaNs, structure (time, latitude, longitude)
    Output:
        - clean_current_array - returns 3d numpy array of same size as current_array_with_nans, with gaps filled
          through linear interpolating at a single point through time.
    Example:
        >>> u_clean = interp_time_surface_currents(u)
    """
    current_array = current_array_with_nans.copy()
    nan_location = ~np.isfinite(current_array)
    fraction_bad_data = np.sum(nan_location, axis=(1, 2)) / (nan_location.shape[1]*nan_location.shape[2])
    keep_surface_current_index = np.where(fraction_bad_data <= .9)

    for i in range(current_array.shape[2]): # interate through x
        for j in range(current_array.shape[1]):
            if np.sum( np.isfinite(current_array[:,j,i])) != 0: # This is out of the radar range, all data is bad
                array_with_nans = current_array[keep_surface_current_index,j,i].flatten()
                nan_index, drop_nan_lambda = nan_helper(array_with_nans)
                array_with_nans[nan_index] = np.interp(drop_nan_lambda(nan_index), drop_nan_lambda(~nan_index), array_with_nans[~nan_index])
                current_array[keep_surface_current_index,j,i] = array_with_nans
    return current_array

def get_griddist(lon,lat):
    """Compute the distance between central grid points in meters from the input grid of latitude and longitude.
    We will assume that the grid is small enought to be approximated by a flat plane.
    Input:
        longitude grid, latitude grid
    Output:
        east_meters, north_meters
    """
    numlons=len(lon)
    numlats=len(lat)
    # remove one since for the zero offset index
    middle_lat_index=int(numlats/2)-1
    middle_lon_index=int(numlons/2)-1
    east1,north1=lonlat2km(lon[middle_lon_index],lat[middle_lat_index],lon[middle_lon_index+1],lat[middle_lat_index])
    east2,north2=lonlat2km(lon[middle_lon_index],lat[middle_lat_index],lon[middle_lon_index],lat[middle_lat_index+1])
    east_meters=int(round(east1*1000))
    north_meters=int(round(north2*1000))
    return east_meters,north_meters


def lonlat2km(lon1,lat1,lon2,lat2):
    """   This function will convert longitude/latitude pairs to distances in 
    kilometers east and west of a reference longitude/latitude point.  The
    equation was obtained from Bowditch's book "The American Practical 
    Navigator, 1995 edition, page 552.
    input: lon1,lat1,lon2,lat2 where lon1,lat1 is the first point. lon2,lat2 is the second point
    output: east,north distance in km between the points
    """
    con=radians(lat1)
    ymeter=111132.92-559.8*np.cos(2*con)+1.175*np.cos(4*con)-0.0023*np.cos(6*con)
    xmeter=111412.84*np.cos(con)-93.5*np.cos(3*con)+0.0118*np.cos(5*con)
    east=(lon2-lon1)*xmeter/1000
    north=(lat2-lat1)*ymeter/1000
    return east,north
def radians(d):
    r=d*np.pi/180
    return r

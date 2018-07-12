# -*- coding: utf-8 -*-
import numpy as np


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
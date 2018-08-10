import numpy as np
def lonlat2km(lon1,lat1,lon2,lat2):
    con=radians(lat1)
    ymeter=111132.92-559.8*np.cos(2*con)+1.175*np.cos(4*con)-0.0023*np.cos(6*con)
    xmeter=111412.84*np.cos(con)-93.5*np.cos(3*con)+0.0118*np.cos(5*con)
    east=(lon2-lon1)*xmeter/1000
    north=(lat2-lat1)*ymeter/1000
    return east,north
def radians(d):
    r=d*np.pi/180
    return r

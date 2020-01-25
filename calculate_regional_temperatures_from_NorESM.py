#!/usr/bin/env python
# Create couple4.regt
# Interpolate and filter NorESM regional temperature data to DIAM grid.

#-----------------------------------------------------------------------

# Modules needed:
from netCDF4 import Dataset
import xarray as xr
import numpy as np
from scipy import interpolate
from scipy.signal import savgol_filter
from optparse import OptionParser
from datetime import timedelta
import pandas as pd

#-----------------------------------------------------------------------

def main():
        # Get parameters
    usage = """usage: %prog --input=file.nc --diam=diam.gin
                            --output=NorESM_spatial_temperature.txt
                            [--output_smooth=NorESM_smoothed_spatial_temperature.txt]"""
    
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input", dest="file",
                      help="netCDF input file ", metavar="file" )
    parser.add_option("-d", "--diam", dest="filename",
                      help="DIAM input file (txt)", metavar="filename" )
    parser.add_option("--output", dest="file_spatial",
                      help="Output txt filename", metavar="file_spatial")
    parser.add_option("--output_smoothed", dest="file_smoothed_spatial",
                      help="Output txt filename with smoothed values (optional)", metavar="file_smoothed_spatial")

    (options, args) = parser.parse_args()

    if not options.file:
        parser.error("First input file must be specified!")
    else:
        file = options.file
    
    if not options.filename:
        parser.error("DIAM input must be specified!")
    else:
        filename = options.filename

    if not options.file_spatial:
        file_spatial = "NorESM_spatial_temperature.txt"
    else:
        file_spatial = options.file_spatial

    if not options.file_smoothed_spatial:
        smooth = False
    else:
        smooth = True
        file_smoothed_spatial = options.file_smoothed_spatial
        
    # Read NorESM data:
    
    dset_atm = xr.open_dataset(file, decode_times=True, use_cftime=True)
    
    # align times (times are apparently wrong e.g. we need to retrieve one day to each date)
    time1 = dset_atm.time.copy()
    for itime in range(time1.sizes['time']):
        bb = dset_atm.time.values[itime].timetuple()
        time1.values[itime] = dset_atm.time.values[itime] - timedelta(days=1)
    dset_atm = dset_atm.assign_coords({'time':time1})

    # Compute yearly average
    avg_temp = dset_atm.groupby('time.year').mean('time')
    
    # Stagger NorESM grid to match DIAM:
    stag_lat = (dset_atm.lat[1] - dset_atm.lat[0]) / 2
    stag_lon = (dset_atm.lon[1] - dset_atm.lon[0]) / 2

    stag_noresm_lat = dset_atm.lat - stag_lat
    stag_noresm_lon = dset_atm.lon - stag_lon 

    # Create lat and lon with spacing that fits with the DIAM data (1x1):
    interp_lon = np.arange(0., 360., 1)  # DIAM has lon -180 to 179
    interp_lat = np.arange(-90., 90., 1)

    #-----------------------------------------------------------------------

    # Define numbers for final grid:
    T = avg_temp.TREFHT.shape[0]    # number of timesteps
    X = interp_lon.shape[0]  # number of longitudes
    Y = interp_lat.shape[0]  # number of latitudes

    interp_temp = np.zeros((T, Y, X))

    # Interpolate data to 1x1:
    for i in range(T):
        f = interpolate.interp2d(stag_noresm_lon,
                             stag_noresm_lat,
                             avg_temp.TREFHT[i,:,:],
                             kind='linear')
        interp_temp[i,:,:] = f(interp_lon, interp_lat) - 273.15

    #-----------------------------------------------------------------------

    # Read coordinates from DIAM:
    data = pd.read_csv(filename, header=None, sep="\s+")
    
    diam_lat = data[2]
    diam_lon = data[3]

    #-----------------------------------------------------------------------

    # Shift coordinates to match DAIM:
    final_interp_temp = np.zeros(interp_temp.shape)
    final_interp_temp[:, :, 0:180] = interp_temp[:, :, 180:]
    final_interp_temp[:, :, 180:] = interp_temp[:, :, 0:180]

    # Write NorESM temperature to file:

    interp_lon = interp_lon - 180

    # Open file for writing:
    file1 = open(file_spatial, 'w')
    if smooth:
        file2 = open(file_smoothed_spatial, 'w')

    for lat, lon in zip(diam_lat, diam_lon):

        index_lat = np.where(interp_lat == lat)[0]
        index_lon = np.where(interp_lon == lon)[0]

        a = final_interp_temp[:, index_lat, index_lon]
        temp = np.reshape(a, T)  # from (100,1) to (100)

        file1.writelines(['%6.i' % item for item in [lat, lon]])
        file1.writelines(['%15.8f' % item for item in temp])
        file1.write('\n')

        if smooth:
            # Smooth the data:
            wl = int(len(temp)/2)*2-1
            smooth_temp = savgol_filter(temp, window_length=wl, 
                                polyorder=1, mode='interp')

            file2.writelines(['%6.i' % item for item in [lat, lon]])
            file2.writelines(['%15.8f' % item for item in smooth_temp])
            file2.write('\n')

    file1.close()
    if smooth:
        file2.close()


#-----------------------------------------------------------------------

if __name__ == "__main__":
    main()
    

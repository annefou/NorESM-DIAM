# Create couple4.regt
# Interpolate and filter NorESM regional temperature data to DIAM grid.

#-----------------------------------------------------------------------

# Modules needed:
from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
from scipy.signal import savgol_filter

#-----------------------------------------------------------------------

# Specify file to calculate from:
file = '/home/jennybj/uio/nird/noresm_cases/post/DIAM_short_test_without_i1/DIAM_short_test_without_i1.TREFHT.nc'
gw_file = '/home/jennybj/Documents/etterbehandling/filer/gw_f19_tn14.nc'

#-----------------------------------------------------------------------

# Read in latitude weights:
ncfile = Dataset(gw_file, 'r')
gw = ncfile.variables['gw'][:]
ncfile.close()

#-----------------------------------------------------------------------

# Read NorESM data:
ncfile = Dataset(file, 'r')

noresm_temp = ncfile.variables['TREFHT'][:]
noresm_lat = ncfile.variables['lat'][:]  # -90 to 90
noresm_lon = ncfile.variables['lon'][:]  # 0 to 360

ncfile.close()

#-----------------------------------------------------------------------

# Calculate yearly average:

years = int(noresm_temp.shape[0]/12)
avg_temp = np.zeros((years, len(noresm_lat), len(noresm_lon)))

for year in range(years):
    avg_temp[year, :, :] = np.average(noresm_temp[year*12:(year + 1)*12,:,:], axis=0)

#-----------------------------------------------------------------------

# Stagger NorESM grid to match DIAM:
stag_lat = (noresm_lat[1] - noresm_lat[0]) / 2
stag_lon = (noresm_lon[1] - noresm_lon[0]) / 2

stag_noresm_lat = noresm_lat - stag_lat
stag_noresm_lon = noresm_lon - stag_lon

# Create lat and lon with spacing that fits with the DIAM data (1x1):
interp_lon = np.arange(0., 360., 1)  # DIAM has lon -180 to 179
interp_lat = np.arange(-90., 90., 1)

#-----------------------------------------------------------------------

# Define numbers for final grid:
T = avg_temp.shape[0]    # number of timesteps
X = interp_lon.shape[0]  # number of longitudes
Y = interp_lat.shape[0]  # number of latitudes

interp_temp = np.zeros((T, Y, X))

# Interpolate data to 1x1:
for i in range(T):
    f = interpolate.interp2d(stag_noresm_lon,
                             stag_noresm_lat,
                             avg_temp[i,:,:],
                             kind='linear')
    interp_temp[i,:,:] = f(interp_lon, interp_lat) - 273.15

#-----------------------------------------------------------------------

# Read coordinates from DIAM:

filename = '/home/jennybj/Documents/etterbehandling/spesifikt/DIAM/data/parse2.gin'

# Read lines from file:
with open(filename, 'r') as myfile:
    lines = myfile.readlines()

rows = len(lines)
columns = len(lines[0].split()) - 2

# Array of all corrdinates from DIAM:
diam_lat = np.zeros(rows)
diam_lon = np.zeros(rows)

# Add coordinates:
for i in range(rows):

    line = lines[i].split()
    index = line.index('"')
    
    lat = float(line[index+1])
    lon = float(line[index+2])

    diam_lat[i] = lat
    diam_lon[i] = lon


#-----------------------------------------------------------------------

# Shift coordinates to match DAIM:
final_interp_temp = np.zeros(interp_temp.shape)
final_interp_temp[:, :, 0:180] = interp_temp[:, :, 180:]
final_interp_temp[:, :, 180:] = interp_temp[:, :, 0:180]

# Write NorESM temperature to file:

interp_lon = interp_lon - 180

# Open file for writing:
file1 = open('NorESM_spatial_temperature.txt', 'w')
file2 = open('NorESM_smoothed_spatial_temperature.txt', 'w')

for lat, lon in zip(diam_lat, diam_lon):

    index_lat = np.where(interp_lat == lat)[0]
    index_lon = np.where(interp_lon == lon)[0]

    a = final_interp_temp[:, index_lat, index_lon]
    temp = np.reshape(a, T)  # from (100,1) to (100)

    file1.writelines(['%6.i' % item for item in [lat, lon]])
    file1.writelines(['%15.8f' % item for item in temp])
    file1.write('\n')

    # Smooth the data:
    smooth_temp = savgol_filter(temp, window_length=21, 
                                polyorder=1, mode='interp')

    file2.writelines(['%6.i' % item for item in [lat, lon]])
    file2.writelines(['%15.8f' % item for item in smooth_temp])
    file2.write('\n')

file1.close()
file2.close()


#-----------------------------------------------------------------------


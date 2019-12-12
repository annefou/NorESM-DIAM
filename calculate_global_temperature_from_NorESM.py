# Create couple4.nor
# Calculate the global average temperature change from pre-industrial
# for the years between 2000 and 2099.
# Write these global averages to the correct file format for DIAM.
#
#---------------------------------------------------------------------

# Modules needed:
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

#-----------------------------------------------------------------------

# Specify NorESM-files:
file_without = '/home/jennybj/uio/nird/noresm_cases/post/DIAM_short_test_without_i1/DIAM_short_test_without_i1.TREFHT.nc'
file_pre_indust = 'N1850OCBDRDDMS_f19_tn14_250119.TREFHT.370_389.nc'
gw_file  = '/home/jennybj/Documents/etterbehandling/filer/gw_f19_tn14.nc'	

#-------------------------------------------------------------------------------

def global_average(atm_file, var, gw_file):

	# Read in latitude weights:
	gw_file = Dataset(gw_file, 'r')
	gw 	    = gw_file.variables['gw'][:]
	gw_file.close()

	# Read in variable to average:
	file  = Dataset(atm_file, 'r')
	variable  = file.variables[var][:]
	lat = file.variables['lat'][:]
	lon = file.variables['lon'][:]

	# Calculate spatial average:
	spatial_avg = np.average(variable, axis=2)
	spatial_avg = np.average(variable, axis=1, weights=gw)

	file.close()

	# Calculate annual average:
	year_avg = []
	for i in range(int(spatial_avg.shape[0]/12)):
		year_avg.append(np.average(spatial_avg[i*12:(i+1)*12]))

	return year_avg

#----------------------------------------------------------------------


# Calculate global, annual averages:
yearly_avg_without    = global_average(atm_file=file_without, var='TREFHT', gw_file=gw_file)
yearly_avg_pre_indust = global_average(atm_file=file_pre_indust, var='TREFHT', gw_file=gw_file)

print(np.asarray(yearly_avg_without) - 273.15)
print(len(yearly_avg_without))

# Calculate the pre-industrial temperature:
pre_indust_avg = np.average(yearly_avg_pre_indust)
print(pre_indust_avg - 273.15)

# Calculate change in globale temperatures:
yearly_avg_without = yearly_avg_without - pre_indust_avg

# Define the years:
years_2000 = (np.arange(2000, 2000 + len(yearly_avg_without) + 1))
years = (np.arange(1, len(yearly_avg_without) + 1))


#----------------------------------------------------------------------

# Read NorESM temperature differences to file:

file = open("NorESM_global_temperature.txt", "w")
for index in range(len(yearly_avg_without)-1):
    file.write('%6i %14.8f \n' % (years[index], yearly_avg_without[index]))
file.close()


#-------------------------------------------------------------------

# APPLY SAVITZKY-GOLAY FILTER TO SMOOTH DATA

smooth_temp = savgol_filter(yearly_avg_without, window_length=21, 
							polyorder=1, mode='interp')

#plt.plot(yearly_avg_without)
#plt.plot(smooth_temp)
#plt.show()


#----------------------------------------------------------------------

# Read NorESM smoothed temperature differences to file:

file = open("NorESM_smoothed_global_temperature.txt", "w")
for index in range(len(smooth_temp)-1):
    file.write('%6i %14.8f \n' % (years[index], smooth_temp[index]))
file.close()

#-------------------------------------------------------------------
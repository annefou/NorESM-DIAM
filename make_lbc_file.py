
# Make a LBC file with constant values for other components than CO2

# 1996 in model corresponds to 1990 in DIAM data...
# Now the model starts in 2000, same as NorESM
# File must contain at least one year before the start year of NorESM

#-----------------------------------------------------------------------

# Modules needed:
import Nio
import numpy as np
from pathlib import Path
import shutil
import sys

#-----------------------------------------------------------------------

filename_original = 'LBC_1750-2015_CMIP6_GlobAnnAvg_c180926.nc'
filename  = 'lbc_short_test_without_tax_from_2000_i1.nc'
years = 269		# number of years in file (must be more than the 266 time dim in the original file)


# Delete if file exists:
file = Path(filename)
if file.exists():
	file.unlink()

shutil.copy(filename_original, filename)

# Create new and open new file for writing:
file_lbc = Nio.open_file(filename, 'w')

#------------------------------------------------------------------------

# Read in original variables:
original_file 		= Nio.open_file(filename_original, 'r')
original_time 		= original_file.variables['time'][246:]			# 1996 = 246
original_date 		= original_file.variables['date'][246:]			# 1996 = 246
original_time_bnds 	= original_file.variables['time_bnds'][246:,:]	# 1996 = 246
original_N2O		= original_file.variables['N2O_LBC'][250,:]		# only read 2000 
original_CH4		= original_file.variables['CH4_LBC'][250,:]		# only read 2000 
original_CFC11eq	= original_file.variables['CFC11eq_LBC'][250,:]	# only read 2000 
original_CF2CL2		= original_file.variables['CF2CL2_LBC'][250,:]	# only read 2000 
original_file.close() 

#------------------------------------------------------------------------

# CHANGE ALL TIME RELATED VARIABLES:
print(original_time)
print(original_date)

# Calculate how time and date changes between time steps:
time_diff 		= np.array([original_time[1] - original_time[0], original_time[2] - original_time[1], \
							original_time[3] - original_time[2], original_time[4] - original_time[3]])
date_diff 		= original_date[1] - original_date[0]
time_bnds_diff 	= np.array([original_time_bnds[1,:] - original_time_bnds[0,:], original_time_bnds[2,:] - original_time_bnds[1,:], \
							original_time_bnds[3,:] - original_time_bnds[2,:], original_time_bnds[4,:] - original_time_bnds[3,:]])


# Make new time list:
new_time = [original_time[0]]

for i in range(int(years/4)):
	new_time.append(new_time[-1] + time_diff[0])
	new_time.append(new_time[-1] + time_diff[1])
	new_time.append(new_time[-1] + time_diff[2])
	new_time.append(new_time[-1] + time_diff[3])

# Make new date list:
new_date = [original_date[0]]

for i in range(years - 1):
	new_date.append(int(new_date[-1] + date_diff))

# Make new time bounds array:
new_time_bnds = np.zeros((years,2), dtype=np.float32)
new_time_bnds[0,:] = original_time_bnds[0,:]

idx = 1
for i in range(int(years/4)):
	new_time_bnds[idx,:] = new_time_bnds[idx-1,:] + time_bnds_diff[0,:]
	idx +=1
	new_time_bnds[idx,:] = new_time_bnds[idx-1,:] + time_bnds_diff[1,:]
	idx +=1
	new_time_bnds[idx,:] = new_time_bnds[idx-1,:] + time_bnds_diff[2,:]
	idx +=1
	new_time_bnds[idx,:] = new_time_bnds[idx-1,:] + time_bnds_diff[3,:]
	idx +=1


# Check that length of new time and time bounds equal number of years:
if len(new_time) != years:
	sys.exit("Length of time list is different from number of years specified.")
elif new_time_bnds.shape[0] != years:
	sys.exit("Length of time bounds array is different from number of years specified.")


# Fill time and date variable with values:
file_lbc.variables['time'][:] 		= new_time	
file_lbc.variables['date'][:] 		= new_date
file_lbc.variables['time_bnds'][:] 	= new_time_bnds	

file_lbc.close()	# closing file so that the time dimension is updated				

#------------------------------------------------------------------------

# REPEAT 2000 VALUES FOR CH4, N2O, CFC11_eq AND CF2CL2 FOR ALL YEARS:

# Create arrays of correct dimension:
CH4 	= np.zeros((years, original_CH4.shape[0]), dtype=np.float32)
N2O 	= np.zeros((years, original_CH4.shape[0]), dtype=np.float32)
CFC11eq = np.zeros((years, original_CH4.shape[0]), dtype=np.float32)
CF2CL2 	= np.zeros((years, original_CH4.shape[0]), dtype=np.float32)

for i in range(original_CH4.shape[0]):
	CH4[:,i] 		= original_CH4[i]
	N2O[:,i] 		= original_N2O[i]
	CFC11eq[:,i] 	= original_CFC11eq[i]
	CF2CL2[:,i] 	= original_CF2CL2[i]

print(original_CH4[0])
print(original_N2O[0])
print(original_CFC11eq[0])
print(original_CF2CL2[0])

# Reopening files:
file_lbc = Nio.open_file(filename, 'w')

# Repeat values of other variables than those for time and CO2:

file_lbc.variables['CH4_LBC'][:] 		= CH4
file_lbc.variables['N2O_LBC'][:] 		= N2O
file_lbc.variables['CFC11eq_LBC'][:] 	= CFC11eq
file_lbc.variables['CF2CL2_LBC'][:] 	= CF2CL2


#------------------------------------------------------------------------

# CO2 VALUES:

# Read in CO2 values from DIAM file:
with open('/home/jennybj/Documents/DIAM/couple4.t') as myfile:
	lines=myfile.readlines()

# Extract the wanted values:
atm_co2 = []

for line in lines:
	atm_co2.append(float(line.split()[3]))	

# Convert from GtC to ppm:
GtC_per_ppm = 2.13
atm_co2 = [line/GtC_per_ppm for line in atm_co2]

# Convert from ppm to mol/mol:
atm_co2 = [line*1e-6 for line in atm_co2]

# Extract the wanted years (index 0 = 1990) (adding 4 extra years):
co2vmr = [atm_co2[0]]*4 + atm_co2[0:years-4]

#------------------------------------------------------------------------

print(co2vmr[0:80])

# Make (time,lat) array for CO2 concentrations:
CO2 = np.zeros((years, original_CH4.shape[0]), dtype=np.float32)

for i in range(years):
	CO2[i,:] = co2vmr[i]

# Check that the number of timesteps is correct:
if not len(co2vmr) == years:
	sys.exit("Number of timesteps in file and CO2 array is not the same!")

# Write CO2 values to file:

file_lbc.variables['CO2_LBC'][:] = CO2

#-------------------------------------------------------------------------

file_lbc.close()

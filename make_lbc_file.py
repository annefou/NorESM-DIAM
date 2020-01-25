#!/usr/bin/env python
# Make a LBC file with constant values for other components than CO2

# 1996 in model corresponds to 1990 in DIAM data...
# Now the model starts in 2000, same as NorESM
# File must contain at least one year before the start year of NorESM

#-----------------------------------------------------------------------

# Modules needed:
# Modules needed
import numpy as np
import pandas as pd
import xarray as xr
import cftime
from optparse import OptionParser

def make_lbc_file(filename_in, couple4file, filename_out, nyears, start_year, select_year):
    """ Create a new netCDF file with nyears starting from start_year.
        All variable values except CO2_LBC are filled with select_year values that are repeated 
        nyears from start_year. CO2_LBC values are filled with the CO2 values from DIAM simulation
        that are contained in couple4file (3rd column).
        """
    dset_in = xr.open_dataset(filename_in, decode_times=True, use_cftime=True)
    # Select a single year (start_year)
    dset_tmp = dset_in.sel(time=slice(cftime.DatetimeGregorian(start_year, 1, 1),
                                      cftime.DatetimeGregorian(start_year+1, 1, 1)))
    # Compute new time and time bounds
    time_start = dset_in.sel(time=cftime.DatetimeGregorian(start_year, 1, 1), method='backfill').time
    new_time = [cftime.DatetimeGregorian(year, time_start.dt.month, time_start.dt.day) 
                for year in range(start_year, start_year+nyears)]
    new_time_bnds = [[cftime.DatetimeGregorian(year, 1, 1), cftime.DatetimeGregorian(year+1, 1, 1)] 
                     for year in range(start_year, start_year+nyears)]
    new_date = [float(cftime.DatetimeGregorian(year, time_start.dt.month, time_start.dt.day).strftime(format="%Y%m%d")) 
                for year in range(start_year, start_year+nyears)]
    # Reindex with new times
    dset_out = dset_tmp.reindex({"time":new_time})
    # Fill time and time_bnds
    dset_out['time'] = dset_out.time
    dset_out['time_bnds'][:,:] = new_time_bnds
    dset_out['N2O_LBC'][:,:] = dset_in.sel(time=slice(cftime.DatetimeGregorian(select_year, 1, 1),
                                      cftime.DatetimeGregorian(select_year+1, 1, 1))).N2O_LBC.values
    dset_out['CH4_LBC'][:,:] = dset_in.sel(time=slice(cftime.DatetimeGregorian(select_year, 1, 1),
                                      cftime.DatetimeGregorian(select_year+1, 1, 1))).CH4_LBC.values
    dset_out['CFC11eq_LBC'][:,:] = dset_in.sel(time=slice(cftime.DatetimeGregorian(select_year, 1, 1),
                                      cftime.DatetimeGregorian(select_year+1,1, 1, 1))).CFC11eq_LBC.values
    dset_out['CF2CL2_LBC'][:,:] = dset_in.sel(time=slice(cftime.DatetimeGregorian(select_year, 1, 1),
                                      cftime.DatetimeGregorian(select_year+1, 1, 1))).CF2CL2_LBC.values
    dset_out['date'][:] = new_date
    
    # Read CO2 values from DIAM simulation
    data = pd.read_csv(couple4file, header=None, sep="\s+")
    # Convert CO2 values (3rd column) from GtC to ppm:
    GtC_per_ppm = 2.13
    atm_co2 = data[3].values / GtC_per_ppm
    # Convert from ppm to mol/mol:
    atm_co2 = atm_co2 * 1e-6
    # Fill CO2_LBC with these new CO2 values
    # - The first 4 values are constant and set to the first CO2 value from DIAM simulation 
    #   (start_year to select_year)
    # - All other years are filled with CO2 values from couple4.t
    dset_out['CO2_LBC'][0:4,:] = atm_co2[0]
    for i in range(len(dset_out.time)-4):
        dset_out['CO2_LBC'][4+i,:] = atm_co2[i]
    
    # Save to output netCDF file
    dset_out.to_netcdf(filename_out)
    

def main():
    usage = """usage: %prog --start_year=YYYY --reference_year=YYYY --nyears=n
                            --input=input_file.nc --output=output_file.nc
                            --couple=diam.txt"""
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--start_year", dest="start_year",
                      help="start year YYYY", metavar="start_year" )
    parser.add_option("-r", "--reference_year", dest="reference_year",
                      help="reference year YYYY", metavar="reference_year" )
    parser.add_option("-n", "--nyears", dest="nyears",
                      help="number of years", metavar="nyears")
    parser.add_option("--input", dest="input_file",
                      help="Input netCDF file with reference data", metavar="input")
    parser.add_option("--output", dest="output_file",
                      help="Output netCDF filename", metavar="output")
    parser.add_option("--couple", dest="couple_diam_file",
                      help="DIAM txt output file", metavar="couple")

    (options, args) = parser.parse_args()

    if not options.start_year:
        parser.error("start year must be specified!")
    else:
        start_year = int(options.start_year)
    
    if not options.reference_year:
        reference_year = start_year
    else:
        reference_year = int(options.reference_year)

    if not options.nyears:
        parser.error("Number of years must be specified!")
    else:
        nyears = int(options.nyears)
    
    if not options.input_file:
        input_file = "LBC_1750-2015_CMIP6_GlobAnnAvg_c180926.nc"
    else:
        input_file = options.input_file
        
    if not options.output_file:
        output_file = "lbc_short_test_without_tax_from_" + str(reference_year) + "_i1.nc"
    else:
        output_file = options.output_file
        
    if not options.couple_diam_file:
        couple_diam_file = "couple4.t"
    else:
        couple_diam_file = options.couple_diam_file
        
        
    make_lbc_file(input_file, couple_diam_file, output_file, nyears, start_year, reference_year)
    
if __name__ == "__main__":
    main()

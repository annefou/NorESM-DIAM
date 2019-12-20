# Create couple4.nor
# Calculate the global average temperature change from pre-industrial
# for the years between 2000 and 2099.
# Write these global averages to the correct file format for DIAM.
#
#---------------------------------------------------------------------

# Modules needed:
import xarray as xr
import numpy as np
import cftime
from optparse import OptionParser
from scipy.signal import savgol_filter
#-----------------------------------------------------------------------
def global_average(atm_file, var, gw_file):

    # Read in latitude weights:
    dset_gw = xr.open_dataset(gw_file)

    # Read in variable to average:
    dset_atm = xr.open_dataset(atm_file, decode_times=True, use_cftime=True)

    # Calculate spatial average:
    spatial_avg = np.average(dset_atm[var], axis=2)
    spatial_avg = np.average(dset_atm[var], axis=1, weights=dset_gw['gw'])

    dset_gw.close()
    dset_atm.close()
    
    # Calculate annual average:
    year_avg = []
    for i in range(int(spatial_avg.shape[0]/12)):
        year_avg.append(np.average(spatial_avg[i*12:(i+1)*12]))

    return year_avg

def save_temperature(tempfile, values, interp):
    # APPLY SAVITZKY-GOLAY FILTER TO SMOOTH DATA
    if interp:    
        smooth_temp = savgol_filter(values, window_length=19, 
                            polyorder=1, mode='interp')

        write_temperature(tempfile, smooth_temp)
    else:
        write_temperature(tempfile, values)
        
#----------------------------------------------------------------------
def write_temperature(tempfile, values):
    # Write NorESM temperature differences to file:
    years = (np.arange(1, len(values) + 1))
    file = open(tempfile, "w")
    for index in range(len(values)-1):
        file.write('%6i %14.8f \n' % (years[index], values[index]))
    file.close()

#----------------------------------------------------------------------
def main():
    # Get parameters
    usage = """usage: %prog --input_without=file_without.nc --input_pre_indust=file_pre_indust.nc
                            --weights=wg.nc --output=global_temperature.txt
                            [--output_smooth=global_smoothed_temperature.txt]"""
    
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--input_without", dest="file_without",
                      help="first input file ", metavar="file_without" )
    parser.add_option("-r", "--input_pre_indust", dest="file_pre_indust",
                      help="second input file (pre-indust)", metavar="file_pre_indust" )
    parser.add_option("--weights", dest="weights",
                      help="netCDF file containing weight values", metavar="weights")
    parser.add_option("--output", dest="output_file",
                      help="Output netCDF filename", metavar="output_file")
    parser.add_option("--output_smoothed", dest="output_smoothed_file",
                      help="Output netCDF filename with smoothed values (optional)", metavar="output_file_smoothed")

    (options, args) = parser.parse_args()

    if not options.file_without:
        parser.error("First input file must be specified!")
    else:
        file_without = options.file_without
    
    if not options.file_pre_indust:
        parser.error("Second input file (pre-industrial) must be specified!")
    else:
        file_pre_indust = options.file_pre_indust
  
    if not options.weights:
        parser.error("Weights must be specified!")
    else:
        weights = options.weights

    if not options.output_file:
        output_file = "NorESM_global_temperature.txt"
    else:
        output_file = options.output_file

    
    if not options.output_smoothed_file:
        interp = False
    else:
        interp = True
        output_smoothed_file = options.output_smoothed_file
        
    # Calculate global, annual averages:
    yearly_avg_without    = global_average(atm_file=file_without, var='TREFHT', gw_file=weights)
    yearly_avg_pre_indust = global_average(atm_file=file_pre_indust, var='TREFHT', gw_file=weights)

    print(np.asarray(yearly_avg_without) - 273.15)

    # Calculate the pre-industrial temperature:
    pre_indust_avg = np.average(yearly_avg_pre_indust)
    print("PRE INDUST AVG")
    print(pre_indust_avg - 273.15)

    # Calculate change in global temperatures:
    yearly_avg_without = yearly_avg_without - pre_indust_avg
    print("change in global temperature: ", yearly_avg_without)

    #----------------------------------------------------------------------

    # Write NorESM temperature differences to file:
    save_temperature(output_file, yearly_avg_without, False)

    #-------------------------------------------------------------------
    if interp:
    # APPLY SAVITZKY-GOLAY FILTER TO SMOOTH DATA
        save_temperature(output_smoothed_file, yearly_avg_without, interp)

    #-------------------------------------------------------------------

if __name__ == "__main__":
    main()
# NorESM-DIAM

Coupling NorESM with DIAM macro-economic model

- Scripts for converting DIAM outputs to NorESM inputs (make_lbc_file.py)
- Scripts for converting NorESM outputs to DIAM inputs (calculate_*)


Inputs for workflow:
- number of years (currently 100)
- NorESM case and resolution


In this configuration:
- DIAM runs until it reaches some kind of equilibrium
    - it computes new Gigatonn Carbon and temperature for each location (saved for temperature in these couple4.* files)
    - Use as input parse2.gin
    - it uses economic information from parse2.gin and temperature/Carbon from couple4.nor (temperature change per year), couple4.rgt (regional temperature per year for each location; t in degrees C), couple4.t (Carbon per year; global values; for norESM we use the 4th column)    
- NorESM runs 100 years

# NorESM to DIAM

### Global temperature

```
python calculate_global_temperature_from_NorESM.py --input_without=DIAM_short_test_without_i1.TREFHT.nc --input_pre_indust=N1850OCBDRDDMS_f19_tn14_250119.TREFHT.370_389.nc --weights=gw_f19_tn14.nc --output=NorESM_global_temperature.txt --output_smoothed=NorESM_smoothed_global_temperature.txt
```

### Regional temperature

```
python calculate_regional_temperatures_from_NorESM.py --input=DIAM_short_test_without_i1.TREFHT.nc --diam=/opt/work/DIAM/parse2.gin --output=NorESM_spatial_temperature.txt --output_smoothed=NorESM_smoothed_spatial_temperature.txt
```

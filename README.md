# NorESM-DIAM

Coupling NorESM with DIAM macro-economic model

- Scripts for converting DIAM outputs to NorESM inputs
- Scripts for converting NorESM outputs to DIAM inputs


Inputs for workflow:
- number of years (currently 100)
- NorESM case and resolution


In this configuration:
- DIAM runs until it reaches some kind of equilibrium
    - it computes new Gigatonn Carbon and temperature for each location (saved for temperature in these couple4.* files)
    - Use as input parse2.gin
    - it uses economic information from parse2.gin and temperature/Carbon from couple4.nor (temperature change per year), couple4.rgt (regional temperature per year for each location; t in degrees C), couple4.t (Carbon per year; global values; for norESM we use the 4th column)    
- NorESM runs 100 years

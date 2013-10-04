Python code for evaluating the 2-D wind multipliers

The Python code was produced for the 2-D wind multipliers based on the existing Octave code written by Chris Thomas (2009). To maintainconsistency, the Python code was structured and functioned closely to the Octave code. 
py_north.py - the main code for calculating 2-D multipliers based on the input DEM and the wind direction (north).  
To adopt Chris Thomas (2009) modifications, the code does not rotate matrices but extracts a line for that direction from the array holding the DEM and operates on that line to produce the hill shape multipliers. The calculated line of multipliers is then written back to the same array that holds the DEM data. The smoothing is done upon output. In the end, outputs including smoothed and unsmoothed data are written into ascii files.
It calls the following functions:
     ascii_read.py          # read input DEM data
     make_path.py        # generate indices of a data line depending on the direction
     multiplier_calc.py  # calculate the multipliers for a data line extracted from the dataset
multiplier_calc.py needs to call two functions:        
     Mh.py                     # calculate Mh using the simpler formula modified by C.Thomas 2009   
     findpeaks.py         # get the indices of the ridges in a data line
Directory structure:

../py-files/: store all the Python codes including main and functions.
../python_output/: generated automatically to store output files.

Notes:
Comparisons show small discrepancies between outputs from the Python and the Octave codes. The difference varies depending on the wind direction. The number of inconsistent data is up to 0.13% of the total data (west wind) with the greatest magnitude of 0.02, which is 2% of the minimum magnitude of output data. 
The discrepancies are believed due to the inexact representation of decimal fractions in the IEEE floating point standard and errors introduced when scaling by powers of ten. Regardless of implementation language, similar problems would arise. 
It also should be aware that some Python functions have slightly different functionality. For example, for values exactly halfway between rounded decimal values, Numpy rounds to the nearest even value.

V 0.2 - parallel implementation
Ported to parallel execution - using PyPar for MPI handling. 
To run in parallel mode:
    mpirun -np 8 python topomult.py -i <input file> -o <output path>


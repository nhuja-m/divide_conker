# Assumes there is histogram.py in the src folder
# Also assumed grids are in grids folder and data folder is not empty
from src.histogram import ProcessHistogram
import sys

# Get the datafile name when running the script include the .fits in the end
args = sys.argv[1]

# If the user didn't provide any arguments, print an error message and exit
if len(args) == 0:
    print("Insert name of data file")
    sys.exit(1)

data_file = args # example: 'fnl0sim1.fits'

###
# CAUTION: the 45 in the argument comes from the total number of angular boxes created during the divide step
# Change if necessary
###
fnlx = ProcessHistogram(data_file, 'randoms.fits', 'ex_par.txt', 45)

# stores output in a new folder called "hist"
fnlx.calc_moment()
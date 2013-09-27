##------------------------------------------------------
## This is the main code for calculating 2-D multipliers - wind direction = north 
## input DEM is stored in ../input
## outputs including smoothed and unsmoothed data are written in ../python_output
## more details in readme.txt

## 2013 Python code - W.Jiang
# --------------------------------------------------------
# import supporting modules
# --------------------------------------------------------
import os
from os.path import join as pjoin, exists
import numpy as np
import math
import logging as log
from ascii_read import ElevationData # read DEM data
import make_path       # generate indices of a data line depending on the direction
import multiplier_calc # calculate the multipliers for a data line extracted from the dataset

from files import flStartLog
import itertools
from scipy import signal

__version__ = '0.2 - parallel implementation'

#direction = 'n'
# --------------------------------------------------------
# check direction specification
# define output directory. If it does not exist, make one.
# --------------------------------------------------------
#direction = direction.strip().lower()
#valid = ['n','s','e','w','ne','nw','se','sw']
#if direction not in valid:
#    raise Exception("Error: invalid direction given, must be one of n,s,e,w,ne,nw,se,sw")

def work(input_dem, mh_data_dir, directions):
    
    for direction in balanced(directions):
        #direction = getDirections(directions)
        topomult(input_dem, mh_data_dir, direction)

def topomult(input_dem, mh_data_dir, direction):
    
    if not exists(mh_data_dir):
        os.makedirs(mh_data_dir)
    
    # --------------------------------------------------------
    # read input data using ascii_read
    # note: data was flattened to a single array
    # --------------------------------------------------------
    log.info('Reading data...')
    
    DEM = ElevationData(input_dem)
    nr = int(DEM.nrows)
    nc = int(DEM.ncols)
    xll = DEM.xllcorner
    yll = DEM.yllcorner
    cellsize = DEM.cellsize
    data =  DEM.data.flatten()
    
    log.info('xll = %f' % xll)
    log.info('yll = %f' %  yll)
    log.info('data_spacing = %f' % cellsize)
    # --------------------------------------------------------
    # Compute the starting positions along the boundaries depending on dir 
    # Together, the direction and the starting position determines a line.
    # Note that the starting positions are defined
    # in terms of the 1-d index of the array.
    # --------------------------------------------------------
    
    ####for direction in balanced(['n','s','e','w','ne','nw','se','sw']):
    if len(direction) == 2:
        data_spacing = DEM.cellsize*math.sqrt(2)
    else:
        data_spacing = DEM.cellsize
        
    Mhdata = np.ones(data.shape)
    strt_idx = []    
    if direction.find('n') >= 0:
        strt_idx = np.append(strt_idx, list(range(0, nr * nc, nr)))
    if direction.find('s') >= 0:
        strt_idx =  np.append(strt_idx, list(range(nr - 1, nr * nc, nr)))
    if direction.find('e') >= 0:
        strt_idx =  np.append(strt_idx, list(range((nc - 1) * nr, nr * nc)))
    if direction.find('w') >= 0:
        strt_idx =  np.append(strt_idx, list(range(0, nr)))
       
    # --------------------------------------------------------
    # for the diagonal directions the corner will have 
    # been counted twice so get rid of the duplicate
    # then loop over the data lines (i.e. over the starting positions)
    # --------------------------------------------------------
    strt_idx = np.unique(strt_idx)
    ctr = 1    # counter in order to report progress 
    
    for idx in strt_idx:
        log.debug( 'processing path %3i' % ctr+' of %3i' % len(strt_idx)+', index %5i.' % idx )
       
        # get a line of the data
        # path is a 1-d vector which gives the indices of the data    
        path = make_path.make_path(nr, nc, idx, direction)
        #print path
        #print len(path)
        line = data[path]
        # compute the multipliers
        M = multiplier_calc.multiplier_calc(line, data_spacing)
          
        # write the line back to the data array
        M = M.conj().transpose()
        Mhdata[path] = M[0,].flatten()
        ctr = ctr + 1
    
    # --------------------------------------------------------
    # reshape the result to matrix like 
    # --------------------------------------------------------
    Mhdata = np.reshape(Mhdata, (nc, nr))
    Mhdata = Mhdata.conj().transpose()
    
    # --------------------------------------------------------
    # output unsmoothed data to an ascii file
    # --------------------------------------------------------
    ofn = pjoin(mh_data_dir, 'mh_'+ direction + '.asc')
    log.info( 'outputting unsmoothed data to: %s' % ofn )
    
    fid = open(ofn,'w')
    
    fid.write('ncols         '+str(nc)+'\n')
    fid.write('nrows         '+str(nr)+'\n')
    fid.write('xllcorner     '+str(xll)+'\n')
    fid.write('yllcorner     '+str(yll)+'\n')
    fid.write('cellsize       '+str(cellsize)+'\n')
    fid.write('NOdata_struct_value  -9999\n')
    
    np.savetxt(fid, Mhdata, fmt ='%4.2f', delimiter = ' ', newline = '\n') 
    
    # --------------------------------------------------------
    # output smoothed data to an ascii file
    # --------------------------------------------------------
    
    ofn = pjoin(mh_data_dir, 'mh_'+ direction + '_smooth.asc')
    log.info( 'outputting smoothed data to: %s' % ofn )
     
    fid = open(ofn,'w')
    fid.write('ncols         '+str(nc)+'\n')
    fid.write('nrows         '+str(nr)+'\n')
    fid.write('xllcorner     '+str(xll)+'\n')
    fid.write('yllcorner     '+str(yll)+'\n')
    fid.write('cellsize       '+str(cellsize)+'\n')
    fid.write('NOdata_struct_value  -9999\n')
    
    g = np.ones((3, 3))/9.
    
    mhsmooth = signal.convolve2d(Mhdata, g, mode='same', boundary='fill', 
                                 fillvalue=1)
   
    np.savetxt(fid, mhsmooth, fmt ='%4.2f', delimiter = ' ', newline = '\n') 
    
    fid.close()
    
    log.info('Finished direction %s' % direction)
    
def getDirections(directions):
    for d in balanced(directions):
        msg = 'Calculating multiplier for direction %s' % d
        log.info(msg)
        yield d

def balanced(iterable):
    """
    Balance an iterator across processors.

    This partitions the work evenly across processors. However, it requires
    the iterator to have been generated on all processors before hand. This is
    only some magical slicing of the iterator, i.e., a poor man version of
    scattering.
    """
    P, p = pp.size(), pp.rank()
    return itertools.islice(iterable, p, None, P)

def attemptParallel():
    """
    Attempt to load Pypar globally as `pp`.  If pypar cannot be loaded then a
    dummy `pp` is created.
    """
    global pp

    try:
        # load pypar for everyone

        import pypar as pp

    except ImportError:

        # no pypar, create a dummy one

        class DummyPypar(object):

            def size(self):
                return 1

            def rank(self):
                return 0

            def barrier(self):
                pass

        pp = DummyPypar()


def run():
    
    logfile = 'topomult.log'
    loglevel = 'INFO'
    verbose = True
    attemptParallel()
    if pp.size() > 1 and pp.rank() > 0:
        logfile += '-' + str(pp.rank())
        verbose = False  # to stop output to console

    flStartLog(logfile, loglevel, verbose)
    
    pp.barrier()
    work('../input/dem.asc', '../python_output/',
             ['n','s','e','w','ne','nw','se','sw'])
    pp.barrier()
    
    
if __name__ == '__main__':
    run()
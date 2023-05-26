import astropy.units as u
from astropy import coordinates
from astropy.table import Table
import sys, argparse
import numpy as np
import os
import time
import pandas as pd
import pdb


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=\
        '''
        Creates the finder charts for the whole night, provided a target list file.
        
        Usage: prepare_obs_run.py target_list_filename telescope
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    # parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help="filename of the target list")

    args = parser.parse_args()

    filename = args.filename
    converters = {4: str} #This is so that -00 gets preserved as -00

    targets = pd.read_csv(filename, delimiter= '\s+',header=None,
        usecols=range(0,11),skiprows=1,converters=converters)
    targets=targets.drop_duplicates()

    names = np.array(targets[0])

    RA = np.array([str(x[1]) + ':' +str(x[2]) +':' + str(x[3]) for i,x in targets.iterrows()])
    Dec = np.array([str(x[4]) + ':' +str(x[5]) +':' + str(x[6]) for i,x in targets.iterrows()])

    coords = coordinates.SkyCoord(RA, Dec, unit = (u.hourangle, u.deg))

    mags = np.array(targets[10])

    #Add RA in deg to the dataframe
    targets["RA_deg"] = coords.ra.deg
    
    #sort
    foo = Table(targets.sort_values("RA_deg").drop("RA_deg", 
        axis = 1).to_numpy())

    foo.write('ra_sorted_'+filename,format = 'ascii.fixed_width_no_header',
        bookend=False, delimiter=None, overwrite=True) 

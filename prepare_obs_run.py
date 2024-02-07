"""
This script takes the downloaded target list from YSE PZ and loops through the targets,
producing finder chart and finding offset stars for them. If host is indicated in the input
target list in a line below the target, the host and the PA are included in the finder chart. 

In the end, you get a finder chart for each target, and an updated target list with the PA 
(if host is provided) and the offset stars. 
"""
import astropy.io.ascii as asci
import astropy.units as u
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astropy.table import Table
import sys, argparse
import numpy as np
import os

from finder_chart import get_finder, get_host_PA_and_sep


def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in ra and ':' in dec):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def do_obsrun(filename, telescope, rotate=False, debug=False, 
    max_offset_star_radius=None, max_mag=None):


    table = asci.read(filename, names=('name','ra','dec','mag'))

    coords = np.array([parse_coord(r['ra'], r['dec']) for r in table])
    names = np.array(table['name'].data)
    mags = np.array(table['mag'])

    all_starlist = ''

    skip_host = False
    ######Now, call finder_chart.py
    outimages = []
    for i in range(len(table)):
        if skip_host == False: #this entry is not host, treat as target
            name = names[i]
            ra_deg = coords[i].ra.deg
            dec_deg = coords[i].dec.deg
            mag = mags[i]
            print("Preparing a finder chart for %s."%name)

            if i < len(table)-1: #not the last item
                next_entry = names[i+1]
                if len(next_entry.split('_'))>1 and next_entry.split('_')[0] == name and next_entry.split('_')[1] == 'host':
                        host_ra = coords[i+1].ra.deg
                        host_dec = coords[i+1].dec.deg
                        print("Host galaxy coordinates provided. RA = %.5f Dec = %.5f"%(host_ra, host_dec))

                        skip_host = True #Because of this flag, the next item is skipped.
                else:
                    host_ra = None
                    host_dec = None
            else: #For the last item, no host is provided.
                host_ra = None
                host_dec = None

            # print(host_ra, host_dec)
            if telescope == 'Keck':
                finder_size = 4.0/60 #3 arcmin, LRIS
                max_separation = 3*60 #3 arcmin, LRIS
                min_mag = 13
                max_mag = 18.5

            elif telescope == 'Lick':
                finder_size = 4/60 #4 arcmin, Kast
                max_separation = 3*60 #3 arcmin, Kast
                min_mag = 5
                max_mag = 18
            elif telescope == 'SOAR':
                finder_size = 5/60 #5 arcmin, SOAR
                max_separation = 3*60 #3 arcmin, SOAR
                min_mag = 8
                max_mag = 19.5
            elif telescope == 'Gemini':
                finder_size = 5/60 # 5 arcmin, Gemini
                max_separation = 1*60 # 2 arcmin, Gemini
                min_mag = 11.0
                max_mag = 19.5
            else:
                print("Telescope should be Keck, Lick, SOAR, or Gemini; default=Lick.")
                finder_size = 4/60 #4 arcmin, Kast
                max_separation = 3*60 #3 arcmin, Kast
                min_mag = 11
                max_mag = 17

            if max_mag: max_mag = max_mag

            #Obtain PA and separation from target ra/dec and host ra/dec
            #To do: make it not duplicate for the "get_finder" function.
            if (host_ra is not None) and (host_dec is not None):
                host_pa, host_sep = get_host_PA_and_sep(ra_deg, dec_deg, 
                    host_ra, host_dec)
                pa_offset = 30 #30 degree offset in slit viewing camera PA

            starlist_entry, outimagename = get_finder( ra_deg, dec_deg, name,  finder_size,
                            mag=mag, minmag=min_mag, maxmag=max_mag,
                            num_offset_stars = 3, min_separation = 15,
                            max_separation = None,\
                            host_ra = host_ra, host_dec = host_dec,
                            directory = 'finders',
                            starlist=None, print_starlist=False, 
                            return_starlist = True, debug=debug, 
                            output_format='png')
            if os.path.exists(outimagename):
                outimages.append(outimagename)
            print(starlist_entry)
            all_starlist += starlist_entry
            ### rotated finder chart
            if rotate and ((host_ra is not None) and 
                (host_dec is not None)):
                print("Make rotated finder chart using PA + 30 deg for LRIS.")
                try:
                    finder_name = name+'_finder.png'
                    rotated_finder = name+'_finder_rot.png'
                    off='%.2f'%(host_pa+pa_offset)
                    m=f'convert {finder_name} -virtual-pixel white +distort '
                    m+=f'SRT {off} {rotated_finder}'
                    os.system(m)
                except:
                    print("Check if imagemagick is installed.")
        elif skip_host == True: #In this case, the next item should be a target
            skip_host = False 
    print(all_starlist)
    ###Write to file
    out_file = open(filename.split('.')[0]+'_final.txt', "w")
    out_file.write(all_starlist)
    out_file.close()

    return(outimages)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=\
    '''
    Creates the finder charts for the whole night, provided a target list file.

    Usage: prepare_obs_run.py target_list_filename telescope

    ''', formatter_class=argparse.RawTextHelpFormatter)

    # parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help="filename of the target list")

    parser.add_argument("telescope", type = str,
                        help="Keck, Lick, SOAR, or Gemini")

    parser.add_argument("-r", "--rotate", action="store_true",
                        help="produce rotated finder chart")

    parser.add_argument("-d", "--debug", action="store_true",
                    help="debug mode")

    parser.add_argument("-m", "--max-offset-star-radius", type=float, 
        default=None, help="Maximum offset star radius in arcsec")

    parser.add_argument("--max-mag", type=float, default=None,
        help="Maximum offset star magnitude")

    args = parser.parse_args()
    filename = args.filename
    telescope = args.telescope

    do_obsrun(filename, telescope, rotate=args.rotate,
        debug=args.debug, max_mag=args.max_mag,
        max_offset_star_radius=args.max_offset_star_radius)

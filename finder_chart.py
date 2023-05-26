##This is based on the finder_chart code by Nadia Blagorodnova
##Modified by Kaew Tinyanont
from __future__ import print_function

# For Python 3.0 and later
from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib import request
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import pylab 
from matplotlib.ticker import MaxNLocator
import argparse
import os
import astropy.wcs
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.ndimage.filters import gaussian_filter
import warnings
from astroquery.mast import Catalogs
from astroquery.vizier import Vizier
import requests
warnings.filterwarnings("ignore")

Vizier.ROW_LIMIT = -1

def deg2hour(ra, dec, sep=":"):
    '''
    Transforms the coordinates in degrees into HH:MM:SS DD:MM:SS with the requested separator.
    '''

    if ( type(ra) is str and type(dec) is str ):
        return ra, dec

    c = SkyCoord(ra, dec, frame='icrs', unit='deg')

    ra = c.ra.to_string(unit=u.hourangle, sep=sep, precision=2, pad=True)
    dec = c.dec.to_string(sep=sep, precision=2, alwayssign=True, pad=True)

    return str(ra), str(dec)


def hour2deg(ra, dec):
    '''
    Transforms string HH:MM:SS DD:MM:SS coordinates into degrees (floats).
    '''
    try:
        ra = float(ra)
        dec = float(dec)

    except:
        c = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))

        ra = c.ra.deg
        dec = c.dec.deg

    return ra, dec


def get_offset(ra1, dec1, ra2, dec2):
    '''
    Computes the offset in arcsec between two coordinates.
    The offset is from (ra1, dec1) - generally an offset star to (ra2, dec2) - the fainter target.
    '''

    from astropy.coordinates import SkyCoord
    bright_star = SkyCoord(ra1, dec1, frame='icrs', unit=(u.deg, u.deg))
    target = SkyCoord(ra2, dec2, frame='icrs', unit=(u.deg, u.deg))
    dra, ddec = bright_star.spherical_offsets_to(target)

    return dra.to(u.arcsec).value, ddec.to(u.arcsec).value

def query_sky_mapper_catalogue(ra, dec, radius_deg, minmag=14, maxmag=18.5):
    '''
    Sends a VO query to the SkyMapper catalogue.
    '''

    coord = SkyCoord(ra, dec, unit='deg')

    filters = {'Ngood':'>5',
        'ClassStar':'>0.7'}

    vquery = Vizier(columns=['RAICRS','DEICRS','ClassStar',
        'Ngood','rPSF','e_rPSF'], column_filters=filters)

    radius = 2*radius_deg*u.deg

    resp = vquery.query_region(coord, width=radius, catalog=['II/358/smss'])

    if len(resp)>0:
        catalog = resp[0]
    else:
        return(None)

    newcat = np.zeros(len(catalog), dtype=[("ra", np.double), ("dec", np.double), ("mag", np.float64)])
    newcat["ra"] = catalog["RAICRS"]
    newcat["dec"] = catalog["DEICRS"]
    newcat["mag"] = catalog["rPSF"]

    return newcat

def query_ps1_catalogue(ra, dec, radius_deg, minmag=13, maxmag=18.0, debug = False):
    '''
    Sends a VO query to the PS1 catalogue.
    Filters the result by mangitude (between 15 and 18.5)
    and by the PSF-like shape of the sources.
    '''

    for i in np.arange(5):
        try:
            result = Vizier.query_region(SkyCoord(ra, dec, unit='deg'),
                width=str(radius_deg)+'d', catalog='II/349')
            break
        except requests.exceptions.ConnectionError:
            pass

    if len(result)>0:
        catalog=result[0]
    else:
        return(None)

    if debug:
        nobj = len(catalog)
        mask1 = catalog['rmag']>minmag
        mask2 = catalog['rmag']<maxmag
        mask3 = (np.abs(catalog['imag'].data-catalog['iKmag'].data)<2.0)

        n1 = len(catalog[mask1])
        n2 = len(catalog[mask2])
        n3 = len(catalog[mask3])

        print('r-mag minimum cuts',nobj-n1)
        print('r-mag maximum cuts',nobj-n2)
        print('i-mag Kron cuts',nobj-n3)

    mask = (catalog['Nd']>1) & (catalog['rmag']>minmag) &\
        (catalog['rmag']<maxmag) &\
        (np.abs(catalog['imag'].data-catalog['iKmag'].data)<2.0)

    if debug: print(catalog["imag"] - catalog["iKmag"])

    catalog = catalog[mask]

    newcat = np.zeros(len(catalog), dtype=[("ra", np.double), ("dec", np.double), ("mag", np.float64)])
    newcat["ra"] = catalog["RAJ2000"]
    newcat["dec"] = catalog["DEJ2000"]
    newcat["mag"] = catalog["rmag"]

    print('There are',len(newcat),'sources')

    return newcat

######################MAST PS1 FUNCTIONS##########################
#https://ps1images.stsci.edu/ps1_dr2_api.html
def ps1cone(ra,dec,radius,table="mean",release="dr1",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a cone search of the PS1 catalog
    
    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2)
    """
    
    data = kw.copy()
    data['ra'] = ra
    data['dec'] = dec
    data['radius'] = radius
    return ps1search(table=table,release=release,format=format,columns=columns,
                    baseurl=baseurl, verbose=verbose, **data)


def ps1search(table="mean",release="dr1",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a general search of the PS1 catalog (possibly without ra/dec/radius)
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
    """
    
    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    checklegal(table,release)
    if format not in ("csv","votable","json"):
        raise ValueError("Bad value for format")
    url = "{baseurl}/{release}/{table}.{format}".format(**locals())
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in ps1metadata(table,release)['name']:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data['columns'] = '[{}]'.format(','.join(columns))

# either get or post works
#    r = requests.post(url, data=data)
    r = requests.get(url, params=data)

    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text


def checklegal(table,release):
    """Checks if this combination of table and release is acceptable
    
    Raises a VelueError exception if there is problem
    """
    
    releaselist = ("dr1", "dr2")
    if release not in ("dr1","dr2"):
        raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
    if release=="dr1":
        tablelist = ("mean", "stack")
    else:
        tablelist = ("mean", "stack", "detection")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))


def ps1metadata(table="mean",release="dr1",
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
    """Return metadata for the specified catalog and table
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    baseurl: base URL for the request
    
    Returns an astropy table with columns name, type, description
    """
    
    checklegal(table,release)
    url = "{baseurl}/{release}/{table}/metadata".format(**locals())
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
               names=('name','type','description'))
    return tab


def mastQuery(request):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object

    Returns head,content where head is the response HTTP headers, and content is the returned data
    """
    
    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content


def resolve(name):
    """Get the RA and Dec for an object using the MAST name resolver
    
    Parameters
    ----------
    name (str): Name of object

    Returns RA, Dec tuple with position"""

    resolverRequest = {'service':'Mast.Name.Lookup',
                       'params':{'input':name,
                                 'format':'json'
                                },
                      }
    headers,resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject['resolvedCoordinate'][0]['ra']
        objDec = resolvedObject['resolvedCoordinate'][0]['decl']
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return (objRa, objDec)

##################################################################

def query_ps1_new_mast(ra, dec, radius_deg, minmag=10, maxmag=18.5,debug = False):
    """
    This function has the same functionalities as query_ps1_catalogue, but uses the more modern astroquery 
    """ 

    # Read RA, Dec and magnitude from CSV 
    # catalog = Table.read("/tmp/ps1_cat.csv", format="ascii.csv", header_start=1)
    results = ps1cone(ra,dec,radius_deg,release='dr2',verbose=debug)
    catalog = asci.read(results)
    if debug:
        print(catalog)

    mask = (catalog["nDetections"]>4) * (catalog["rMeanPSFMag"] > minmag) * (catalog["rMeanPSFMag"] < maxmag) *\
    (catalog["iMeanPSFMag"] - catalog["iMeanKronMag"] < 0.05) #This last one to select stars.
    if debug: print(catalog["iMeanPSFMag"] - catalog["iMeanKronMag"])
    # (catalog["iMeanPSFMag"] - catalog["iMeanKronMag"] < 0.05) #This last one to select stars.

    
    #*(catalog["rMeanPSFMag"] > minmag) * (catalog["rMeanPSFMag"] < maxmag) 
    catalog = catalog[mask]

    # print(catalog.colnames)
    newcat = np.zeros(len(catalog), dtype=[("ra", np.double), ("dec", np.double), ("mag", np.float)])
    newcat["ra"] = catalog["raMean"]
    newcat["dec"] = catalog["decMean"]
    newcat["mag"] = catalog["rMeanPSFMag"]
    
    return newcat
            

def get_host_PA_and_sep(ra, dec, host_ra, host_dec):
    """
    Takes the coordinates of the transient and the host and compute the PA and offset

    Parameters
    ----------
    ra : float
        RA of our target in degrees.
    dec : float
        DEC of our target in degrees.
    host_ra : float
        RA of the host galaxy in degrees.
    host_dec : float
        DEC of the host galaxy in degrees.
    slit_length : float (optional) ---> Remove this for now. Observer should be able to judge
        The length of the spectrograph's slit in arcsec. 
        If offset > slit_length, warning is provided.
    
    Outputs:
    --------
    Return PA in degrees and separation between target and host in arcsec
    """
    
    #Define SkyCoord object for target and host
    transient_coord = SkyCoord(ra, dec, frame='icrs', unit='deg')
    host_coord = SkyCoord(host_ra, host_dec, frame='icrs', unit='deg')
    
    #Compute separation and PA
    pa = transient_coord.position_angle(host_coord)
    sep = transient_coord.separation(host_coord)
    
    return pa.deg, sep.arcsec

def get_skymapper_image(ra, dec, rad):

    url="http://skymappersiap.asvo.nci.org.au/dr1_cutout/query?POS=%.6f,%.6f&SIZE=%.3f&FORMAT=image/fits&INTERSECT=center&RESPONSEFORMAT=CSV"%(ra, dec, rad)
    page = urlopen(url)
    content = page.read()
    f = open("/tmp/skymapper_image_index.csv", "wb")
    f.write(content)
    f.close()

    ix = Table.read('/tmp/skymapper_image_index.csv', format="ascii.csv")

    mask = ((ix['band']=='r')|(ix['band']=='g'))

    ix = ix[mask]
    ix.sort(keys='exptime')

    image_url = ix['get_image'][-1]
    urlretrieve(image_url, '/tmp/tmp.fits')

def get_fits_image(ra, dec, rad, debug=False):
    '''
    Connects to the PS1 or SkyMapper image service to retrieve the fits file to be used as a bse for the finder chart.
    '''
    #If dec> -30, we have Pan-STARRS
    if dec > -35:
        # Construct URL to download Pan-STARRS image cutout, and save to tmp.fits

        # First find the index of images and retrieve the file of the image that we want to use.
        image_index_url = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=r'.format(ra, dec)
        urlretrieve(image_index_url, '/tmp/ps1_image_index.txt')
        ix = Table.read('/tmp/ps1_image_index.txt', format="ascii")
        if len(ix)==0:
            get_skymapper_image(ra, dec, rad)
        else:
            f = ix['filename'].data[0]

            image_url = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red={0}&format=fits&size={1}&ra={2}&dec={3}".format(f, int(np.round(rad*3600*4, 0)), ra, dec)
            if (debug):
                print("URL:", image_url)
                print("Downloading PS1 r-band image...")

            #Store the object to a fits file.
            urlretrieve(image_url, '/tmp/tmp.fits')


    #Otherwise, we have SkyMapper
    else:
        get_skymapper_image(ra, dec, rad)


    #Finally, once we have Pan-STARRS or SkyMapper images, we try to open them.
    #If there has been any problem with that, we will just go to the DSS image service.
    try:
        image = fits.open("/tmp/tmp.fits")
        #If everything went well, it shall be a fits image and opening it shall cause no issue.
        return '/tmp/tmp.fits'

        #If there was an error with the fits, we shall go for the DSS image
    except IOError:
        #One of the services may fail, so we need to account for that and provide a backup DSS image service.
        try:
            image_url = 'http://archive.eso.org/dss/dss/image?ra=%.5f&dec=%.5f&x=%.2f&y=%.2f&Sky-Survey=DSS1&mime-type=download-fits' % \
                ((ra), (dec), (rad*60), (rad*60))
            if debug: print("Downloading DSS image...")
            urlretrieve(image_url, '/tmp/tmp.fits')
        except:
            image_url = 'http://archive.stsci.edu/cgi-bin/dss_search?ra=%.6f&dec=%.6f&generation=DSS2r&equinox=J2000&height=%.4f&width=%.4f&format=FITS' % \
                (ra, dec, rad*60, rad*60)
            urlretrieve(image_url, '/tmp/tmp.fits')

        #We try one more time to open it. If not successful, we return None as the image filename.
        try:
            fits.open("/tmp/tmp.fits")
        except IOError:
            print("Your fits image could not be retrieved.")
            return None



def get_cutout(ra, dec, name, rad, debug=True):
    '''
    Obtains the color composite cutout from the PS1 images.
    '''
    try:
        ra=float(ra)
        dec=float(dec)
    except:
        ra, dec = hour2deg(ra, dec)

    catalog = query_ps1_catalogue(ra, dec, rad, debug = debug)

    if (debug):
        print(catalog)


    # Construct URL to download DSS image cutout, and save to tmp.fits
    image_index_url_red = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=y'.format(ra, dec)
    image_index_url_green = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=i'.format(ra, dec)
    image_index_url_blue = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=g'.format(ra, dec)

    urlretrieve(image_index_url_red, '/tmp/image_index_red.txt')
    urlretrieve(image_index_url_green, '/tmp/image_index_green.txt')
    urlretrieve(image_index_url_blue, '/tmp/image_index_blue.txt')

    ix_red = np.genfromtxt('/tmp/image_index_red.txt', names=True, dtype=None)
    ix_green = np.genfromtxt('/tmp/image_index_green.txt', names=True, dtype=None)
    ix_blue = np.genfromtxt('/tmp/image_index_blue.txt', names=True, dtype=None)

    image_url = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=%s&green=%s&blue=%s&filetypes=stack&auxiliary=data&size=%d&ra=%.6f&dec=%.6f&output_size=256"%\
    (ix_red["filename"], ix_green["filename"], ix_blue["filename"], rad*3600*4, ra, dec)

    if (debug):
        print(image_url)
        print("Downloading PS1 r-band image...")
    urlretrieve(image_url, '/tmp/tmp_%s.jpg'%name)


def get_finder(ra, dec, name, rad, debug=False, starlist=None,
                print_starlist=True,
                return_starlist = False,
                host_ra = None, host_dec = None,
                telescope="Keck", directory=".",
                minmag=13, maxmag=19.5, num_offset_stars = 3,
                min_separation = 1, max_separation = None,
                mag=np.nan, \
                marker = 'circle', source_comments = None,
                output_format = 'png'):
    '''
    Creates a PDF with the finder chart for the object with the specified name and coordinates.
    It queries the PS1 catalogue to obtain nearby offset stars and get an R-band image as background.

    Parameters
    ----------
    ra : float
        RA of our target in degrees.
    dec : float
        DEC of our target in degrees.
    name : str
        The name of your target
    rad : float
        Search radius for the finder in degrees.
    debug : bool (optional)
        Option to activate/ deactivate additional output.
    starlist : str (optional)
        Path/name of the file where the coordinates for the object are going to be saved.
        If the file exists, the content will be appended at the end.
        If no value is provided, the output just writes to the standard output (in case print_starlist is True).
    print_starlist : boolean (optional)
        Indicates if the starlist shall be printed in the standard output.
    return_starlist : boolean (optional, default is False)
        If true, this function return the strings to be written to a starlist file.
    host_ra, host_dec : float (optional)
        Coordinates of the host galaxy. If provided, the PA and separation between host and target are computed
        and written to the starlist.
    telescope : str (optional)
        The current two values accepted are "Keck" and "P200".
    directory : str (optional)
        The directory where the PDF with the finder chart shall be stored.
        If no value given, the file will be store in the current directory where the script is run.
    minmag : float (optional)
        The minimum magnitud (brightest in this case) star that we would like to use as an offset star.
    maxmag : float (optional)
        The maximum magnitude (faintest) star that we would like to use as an offset star.
    num_offset_stars : int (optional)
        The number of offset stars to include. The default is 3.
    min_separation : float (optional)
        The minimum separation in arcsec between the offset star and the SN. We don't want something too close.
    max_separation : float (optional)
        The maximum separation in arcsec between the offset star and the SN. Default is None.
    mag : float or `None` (optional)
        The magnitude of our target.
    marker : text
        'circle', or 'cross' to mark the star. Target and host galaxy (if provided) are always cross.
    source_comments : text
        If not None, append this comment into the starlist file
    output_format : text
        pdf or png. Only png is supported for rotation right now
    '''

    try:
        ra=float(ra)
        dec=float(dec)
    except:
        ra, dec = hour2deg(ra, dec)

    if dec < -35:
        catalog = query_sky_mapper_catalogue(ra, dec, (rad/2.)*0.95, minmag=minmag, maxmag=maxmag)
    else:
        catalog = query_ps1_catalogue(ra, dec, (rad/2.)*0.95, minmag=minmag, maxmag=maxmag, debug = debug)

    if (debug):
        print("Catalog returned by the function.")
        print(catalog)


    if (len(catalog)==0):
        if debug: print("Looking for a bit fainter stars up to mag: %.2f"%(maxmag+0.5))
        catalog = query_ps1_catalogue(ra, dec, (rad/2.)*0.95, minmag=minmag, maxmag=maxmag+0.5)

    if (not catalog is None and len(catalog)>0):
        np.random.shuffle(catalog)

    ###########Reject stars that are too close
    offsets_stars = SkyCoord(catalog["ra"], catalog["dec"], unit = (u.deg, u.deg))

    #Define SkyCoord object for target and host
    transient_coord = SkyCoord(ra, dec, frame='icrs', unit='deg')

    offsets = transient_coord.separation(offsets_stars)

    no_self_object = offsets.arcsec > min_separation

    catalog = catalog[no_self_object]
    if (debug):
        print('After removing stars too close.')
        print(catalog)
    
    ###########Reject stars that are too far
    if max_separation is not None:
        # not_too_far = (np.abs(catalog["ra"]-ra)*np.cos(np.deg2rad(dec))<max_separation/3600)*(np.abs(catalog["dec"]-dec)<max_separation/3600)
        not_too_far = offsets.arcsec < max_separation
        if len(catalog[not_too_far]) <= num_offset_stars+1:
            print("Not enough stars, ignore maximum separation and use the whole FoV.")
        else:
            catalog = catalog[not_too_far]

        if (debug):
            print('After removing stars too far.')
            print(catalog)

    ###########printwarning if too few stars, set the number of offset stars to length of catalog
    if num_offset_stars >= 4:
        print("Do you really need this many stars? If n < 7 you should be fine...")
    if len(catalog) <= num_offset_stars:
        print("Only %d offset star(s) available"%len(catalog))
        num_offset_stars = len(catalog)
        # print('get to here')


    if len(catalog) >1:
        catalog.sort(order='mag')

    if (debug): print(catalog)

    ###########Get FITS image of the FoV from DSS
    image_file = get_fits_image(ra, dec, rad, debug=debug)

    if image_file is None:
        print("FATAL ERROR! Your FITS image could not be retrieved.")
        return

    image = fits.open(image_file)

    # Get pixel coordinates of SN, reference stars in DSS image
    wcs = astropy.wcs.WCS(image[0].header)

    ###Get pixel coordinates of offset stars
    ###Now support any provided number of offset stars
    ref_pix = []
    for i in range(num_offset_stars):
        ref_pix += [wcs.wcs_world2pix(np.array([[catalog["ra"][i], catalog["dec"][i]]], np.float64), 1)]

    ###Get pixel coordinates of the target
    target_pix = wcs.wcs_world2pix([(np.array([ra,dec], np.float64))], 1)

    ###Get pixel coordinates of the host (if not None)
    if (host_ra is not None) and (host_dec is not None):
        host_pix = wcs.wcs_world2pix([(np.array([host_ra,host_dec], np.float64))], 1)

    # Plot finder chart

    #Adjust some of the counts to make easier the plotting.
    image[0].data[image[0].data>30000] = 30000
    image[0].data[np.isnan(image[0].data)] = 0

    fig, ax = plt.subplots(1,1,figsize = (8,6), subplot_kw={'projection':wcs})

    plt.set_cmap('gray_r')
    smoothedimage = gaussian_filter(image[0].data, 1.3)
    ax.imshow(smoothedimage, origin='lower',vmin=np.percentile(image[0].data.flatten(), 10), \
                vmax=np.percentile(image[0].data.flatten(), 99.0))

    # Mark target
    plt.plot([target_pix[0,0]+35,(target_pix[0,0]+10)],[target_pix[0,1],(target_pix[0,1])], 'g-', lw=2)
    plt.plot([target_pix[0,0],(target_pix[0,0])],[target_pix[0,1]+10,(target_pix[0,1])+35], 'g-', lw=2)
    plt.annotate(name, xy=(target_pix[0,0], target_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points',
        color='green')

    # Mark host, if applicable
    if (host_ra is not None) and (host_dec is not None):
        #Get pa and sep
        host_pa, host_sep = get_host_PA_and_sep(ra,dec, host_ra,host_dec)
        if debug:
            print(host_pa, host_sep )
        plt.plot([host_pix[0,0]-35,(host_pix[0,0]-10)],[host_pix[0,1],(host_pix[0,1])], 'c-', lw=1)
        plt.plot([host_pix[0,0],(host_pix[0,0])],[host_pix[0,1]-10,(host_pix[0,1])-35], 'c-', lw=1)
        plt.annotate('Host', xy=(host_pix[0,0], host_pix[0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points', color = 'c')

    cols = ["b","r","orange","m", "brown", "k", "c"]
    if num_offset_stars > 7:
        print("Why do you need so many offset stars!?")
    for i in range(num_offset_stars):
        if marker == 'cross':
            plt.plot([ref_pix[i][0,0]+35,(ref_pix[i][0,0]+10)],[ref_pix[i][0,1],(ref_pix[i][0,1])], ls = '-',color = cols[i], lw=2)
            plt.plot([ref_pix[i][0,0], (ref_pix[i][0,0]) ] , [ref_pix[i][0,1]+10,(ref_pix[i][0,1])+35], ls = '-',color = cols[i], lw=2)
        if marker == 'circle':
            # print((ref_pix[i][0,0], ref_pix[i][0,1]))
            circ = plt.Circle((ref_pix[i][0,0], ref_pix[i][0,1]), 35, color=cols[i], fill=False, lw = 3)
            ax = plt.gca()
            ax.add_artist(circ)
        plt.annotate("S%d"%(i+1), xy=(ref_pix[i][0,0], ref_pix[i][0,1]),  xycoords='data',xytext=(22,-3), textcoords='offset points', color=cols[i])

    # Set limits to size of DSS image
    ax.set_xlim([0,(image[0].data.shape[0])])
    ax.set_ylim([0,(image[0].data.shape[1])])

    # Plot compass
    plt.plot([(image[0].data.shape[0])-10,(image[0].data.shape[0]-40)],[10,10], 'k-', lw=2)
    plt.plot([(image[0].data.shape[0])-10,(image[0].data.shape[0])-10],[10,40], 'k-', lw=2)
    plt.annotate("N", xy=((image[0].data.shape[0])-20, 40),  xycoords='data',xytext=(-4,5), textcoords='offset points')
    plt.annotate("E", xy=((image[0].data.shape[0])-40, 20),  xycoords='data',xytext=(-12,-5), textcoords='offset points')

    # Set axis tics (not implemented correctly yet)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.set_xlabel('%.1f\''%(rad*60))
    ax.set_ylabel('%.1f\''%(rad*60))

    # Set size of window (leaving space to right for ref star coords)
    plt.subplots_adjust(right=0.65,left=0.05, top=0.99, bottom=0.05)

    # List name, coords, mag of references etc
    plt.text(1.02, 0.85, name+" mag=%.1f"%mag, transform=ax.transAxes, fontweight='bold')
    plt.text(1.02, 0.80, "%.5f %.5f"%(ra, dec),transform=ax.transAxes)
    
    rah, dech = deg2hour(ra, dec)
    plt.text(1.02, 0.75,rah+"  "+dech, transform=ax.transAxes)
    
    #If host info is provided
    if (host_ra is not None) and (host_dec is not None):
        plt.text(1.02, 0.7, "Host PA = %.2f deg, Sep = %.2f\""%(host_pa, host_sep), transform=ax.transAxes)

    for i in range(num_offset_stars):
        ofR = get_offset(catalog["ra"][i], catalog["dec"][i], ra, dec)
        S = deg2hour(catalog["ra"][i], catalog["dec"][i], sep=":")

        plt.text(1.02, 0.60-0.2*i,"S%d, mag=%.2f"%(i+1,catalog["mag"][i]), transform=ax.transAxes, color=cols[i])
        plt.text(1.02, 0.55-0.2*i,'%s %s'%(S[0], S[1]), transform=ax.transAxes, color=cols[i])
        plt.text(1.02, 0.5 -0.2*i,"E: %.2f N: %.2f"%(ofR[0], ofR[1]),transform=ax.transAxes, color=cols[i]) 



    

    if not os.path.exists(directory):
        os.makedirs(directory)

    # # Save to pdf
    outimagename = ''
    if output_format == 'pdf':
        outimagename = os.path.join(directory, str(name+'_finder.pdf'))
        pylab.savefig(os.path.join(directory, str(name+'_finder.pdf')), bbox_inches = 'tight')
        if debug: print("Saved to %s"%os.path.join(directory, str(name+'_finder.pdf')))

    # Save to png
    if output_format == 'png':
        outimagename = os.path.join(directory, str(name+'_finder.png'))
        pylab.savefig(os.path.join(directory, str(name+'_finder.png')), bbox_inches = 'tight', dpi = 150)
        if debug: print("Saved to %s"%os.path.join(directory, str(name+'_finder.png')))

    pylab.close("all")
    
    if telescope == "Keck":
        commentchar = "#"
        separator = ""
    else:
        commentchar = "!"
        separator = "!"
        
    ######printoffset stars to a starlist file
    #If no magnitude was supplied, just do not put it on the chart.
    if not np.isnan(mag):
        rmag = "r=%.2f"%mag
    else:
        rmag = ""

    # if print_starlist or not starlist is None:
    if print_starlist:
        #Target if no host
        if (host_ra is None) and (host_dec is None):
            print( "{:s}{:s} {:s}  2000.0  {:s} {:s}".format(name.ljust(16), *deg2hour(ra, dec, sep=" "), commentchar, rmag ) )
        else:
            print( "{:s}{:s} {:s}  2000.0  rotmode=pa rotdest={:.2f} {:s} {:s".format(name.ljust(16), *deg2hour(ra, dec, sep=" "), host_pa, commentchar, rmag ) )
        #offset stars
        for i in range(num_offset_stars):
            S = deg2hour(catalog["ra"][i], catalog["dec"][i], sep=" ")
            ofR = get_offset(catalog["ra"][i], catalog["dec"][i], ra, dec)
            if (host_ra is None) and (host_dec is None):
                print( "{:s}{:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} {:s} r={:.1f} ".format( (name+"_S%d"%(i+1)).ljust(16), S[0], S[1], separator, ofR[0], ofR[1], commentchar, catalog["mag"][i]))
            else:
                print( "{:s}{:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} rotmode=pa rotdest={:.2f} {:s} r={:.1f} ".format( (name+"_S%d"%(i+1)).ljust(16), S[0], S[1], separator, ofR[0], ofR[1], host_pa, commentchar, catalog["mag"][i]))
    if return_starlist == True:
        if (host_ra is None) and (host_dec is None):
            starlist_str = "{:s}{:s} {:s}  2000.0  {:s} {:s}\n".format(name.ljust(16), *deg2hour(ra, dec, sep=" "), commentchar, rmag )
        else:
            starlist_str = "{:s}{:s} {:s}  2000.0  rotmode=pa rotdest={:.2f} {:s} {:s}\n".format(name.ljust(16), *deg2hour(ra, dec, sep=" "), host_pa, commentchar, rmag)
        for i in range(num_offset_stars):
            S = deg2hour(catalog["ra"][i], catalog["dec"][i], sep=" ")
            ofR = get_offset(catalog["ra"][i], catalog["dec"][i], ra, dec)
            if (host_ra is None) and (host_dec is None):
                starlist_str += "{:s}{:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} {:s} r={:.1f} \n".format( (name+"_S%d"%(i+1)).ljust(16), S[0], S[1], separator, ofR[0], ofR[1], commentchar, catalog["mag"][i])
            else:
                starlist_str += "{:s}{:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} rotmode=pa rotdest={:.2f} {:s} r={:.1f} \n".format( (name+"_S%d"%(i+1)).ljust(16), S[0], S[1], separator, ofR[0], ofR[1], host_pa, commentchar, catalog["mag"][i])
 
    #Write to the starlist if the name of the starlist was provided.
    if (not starlist is None) and (telescope =="Keck"):
        with open(starlist, "a") as f:
            if (host_ra is None) and (host_dec is None):
                starlist_str2 = "{:s}{:s} {:s}  2000.0  {:s} {:s}\n".format(name.ljust(16), *deg2hour(ra, dec, sep=" "), commentchar, rmag )
            else:
                starlist_str2 = "{:s}{:s} {:s}  2000.0  rotmode=pa rotdest={:.2f} {:s} {:s}\n".format(name.ljust(16), *deg2hour(ra, dec, sep=" "), host_pa, commentchar, rmag)
            for i in range(num_offset_stars):
                S = deg2hour(catalog["ra"][i], catalog["dec"][i], sep=" ")
                ofR = get_offset(catalog["ra"][i], catalog["dec"][i], ra, dec)
                if (host_ra is None) and (host_dec is None):
                    starlist_str2 += "{:s}{:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} {:s} r={:.1f} \n".format( (name+"_S%d"%(i+1)).ljust(16), S[0], S[1], separator, ofR[0], ofR[1], commentchar, catalog["mag"][i])
                else:
                    starlist_str2 += "{:s}{:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} rotmode=pa rotdest={:.2f} {:s} r={:.1f} \n".format( (name+"_S%d"%(i+1)).ljust(16), S[0], S[1], separator, ofR[0], ofR[1], host_pa, commentchar, catalog["mag"][i])
            f.write(starlist_str2)

    if (not starlist is None) and (telescope =="P200"):
        with open(starlist, "a") as f:
            f.write( "{0} {1} {2}  2000.0 ! \n".format(name.ljust(19), *deg2hour(ra, dec, sep=" ")) )
            # if (len(catalog)>0):
            #     f.write ( "{:s} {:s} {:s}  2000.0 ! raoffset={:.2f} decoffset={:.2f} r={:.1f}  \n".format( (name+"_S1").ljust(19), S1[0], S1[1], ofR1[0], ofR1[1], catalog["mag"][0]))
            # if (len(catalog)>1):
            #     f.write ( "{:s} {:s} {:s}  2000.0 ! raoffset={:.2f} decoffset={:.2f} r={:.1f}  \n".format( (name+"_S2").ljust(19), S2[0], S2[1], ofR2[0], ofR2[1], catalog["mag"][1]))
            for i in range(num_offset_stars):
                S = deg2hour(catalog["ra"][i], catalog["dec"][i], sep=" ")
                ofR = get_offset(catalog["ra"][i], catalog["dec"][i], ra, dec)
                f.write ( "{:s} {:s} {:s}  2000.0 ! raoffset={:.2f} decoffset={:.2f} r={:.1f}  \n".format( (name+"_S%d"%(i+1)).ljust(19), S[0], S[1], ofR[0], ofR[1], catalog["mag"][i]))
            f.write('\n')     

    if return_starlist == True:
        return (starlist_str, outimagename)



if __name__ == '__main__':
    

    parser = argparse.ArgumentParser(description=\
        '''
        Creates the finder chart for the given RA, DEC and NAME.
        
        Usage: finder_chart.py <RA> <Dec> <Name> <rad [deg]> <telescope [P200|Keck]>
            
        ''', formatter_class=argparse.RawTextHelpFormatter)
        
        
    #Check if correct number of arguments are given
    if len(sys.argv) < 4:
        print("Usage: finder_chart.py <RA> <Dec> <Name>  <rad [deg]> <telescope [P200|Keck]>")
        sys.exit()
     
    ra=sys.argv[1]
    dec=sys.argv[2]
    name=str(sys.argv[3])
    if (len(sys.argv)>=5):
        rad = float(sys.argv[4])
        if (rad > 30./60):
            print('Search radius larger than 30 arcmin. Not sure why you need such a large finder chart... reducing to 10 armin for smoother operations...')
            rad = 10./60
    else:
        rad = 2./60
        
    print('Using search radius of %.1f arcsec.'%(rad*3600))

        
    if (len(sys.argv)>5):
        telescope = sys.argv[5]
    else:
        telescope = "Keck"
        print('Assuming that the telescope you observe will be Keck. If it is "P200", please specify otherwise.')

    if (len(sys.argv)>6):
        num_offset_stars = int(sys.argv[6])
    else:
        num_offset_stars = 3

    
    get_finder(ra, dec, name, rad, telescope=telescope, debug=False, minmag=13, 
        maxmag=19.5, num_offset_stars=num_offset_stars)

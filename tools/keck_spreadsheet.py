import astropy.units as u
from astropy import coordinates
import sys, argparse
import numpy as np
import os
import time
import pandas as pd
from astropy.time import Time,TimeDelta
from datetime import datetime, timedelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astroplan import Observer
import pdb
import math

# The figure below shows the telescope limits for the Keck telescopes. 
# In the regions where an elevation limit of ~35° is shown, 
# the telescope hits the Nasmyth deck. Note that this is a problem in the South and West for Keck II
# while for Keck I it is more of a problem in the North and East.

def keck1_rising_time(sky_coord, date_UT):
    date_UT = Time(date_UT)
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-7, 7, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)

    obj_altaz = sky_coord.transform_to(observer_frame)
    rising = obj_altaz.az < 180*u.deg #select only times where the object is setting. 
    #now check for Az when Alt is either 

    nasmyth_platform = obj_altaz[rising].alt > 33.3*u.deg #time steps where object could be below the nasmyth platform, K1
    general_shutter  = obj_altaz[rising].alt > 18*u.deg
    # print(np.sum(nasmyth_platform), np.sum(general_shutter))
    #Object never rises
    #Need to be more check cases here. 
    if np.sum(general_shutter) == 0:
        rise_time = None
    else: #object do rise
        # print('foo ', obj_altaz[rising][nasmyth_platform][0])
        if (len(obj_altaz[rising][nasmyth_platform]) > 0) and (obj_altaz[rising][nasmyth_platform][0].az >5.3*u.deg and obj_altaz[rising][nasmyth_platform][0].az < 146.2*u.deg):
            rise_time = rtime(obj_altaz[rising][nasmyth_platform][0].obstime.datetime)
        elif len(obj_altaz[rising][general_shutter]) > 0:
            rise_time = rtime(obj_altaz[rising][general_shutter][0].obstime.datetime)
        else:
            rise_time = None
    return rise_time

def keck2_setting_time(sky_coord, date_UT):
    date_UT = Time(date_UT)
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-7, 7, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)

    obj_altaz = sky_coord.transform_to(observer_frame)
    setting = obj_altaz.az > 180*u.deg #select only times where the object is setting. 
    #now check for Az when Alt is either 

    nasmyth_platform = obj_altaz[setting].alt < 36.8*u.deg #time steps where object could be below the nasmyth platform, K2
    general_shutter  = obj_altaz[setting].alt < 18*u.deg

    #object never sets
    if np.sum(general_shutter) == 0:
        set_time = None
    else: #object sets
        if obj_altaz[setting][nasmyth_platform][0].az >185.3*u.deg and obj_altaz[setting][nasmyth_platform][0].az < 332.8*u.deg:
            set_time = rtime(obj_altaz[setting][nasmyth_platform][0].obstime.datetime)
        else:
            set_time = rtime(obj_altaz[setting][general_shutter][0].obstime.datetime)
    return set_time
    

def rtime(astroplan_datetime):
    # Rounds astroplan datetime to nearest minute by adding a timedelta minute if second >= 30

    t=astroplan_datetime
    risetime=t.replace(second=0, microsecond=0, minute=t.minute, hour=t.hour)\
               +timedelta(minutes=t.second//30)
    return (risetime.strftime("%H:%M"))


def keck_wrap_time(sky_coord, date_UT):

    # Technically this should be split into Keck I and Keck II
    # Here I am ignoring the Nazymth platform 
    # Just making sure than objects are >18 degree alt 
    #
    
    date_UT = Time(date_UT)
    utc_offset = -10*u.hour #Hawaii standard time
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-7, 7, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)
    
    obj_altaz = sky_coord.transform_to(observer_frame)
    
    # Get Azimuth ranges for targets
    # 0 is North, 90 is East
    # North wrap is: 325-90
    # South Wrap is : 90-235
    # North or South is 235-325
    nwrap=np.logical_or((obj_altaz.az >= 325*u.deg) & (obj_altaz.az <= 360*u.deg), (obj_altaz.az >= 0*u.deg) & (obj_altaz.az <= 90*u.deg))
    swrap=(obj_altaz.az >= 90*u.deg) & (obj_altaz.az <= 235*u.deg)
    nswrap=(obj_altaz.az > 235*u.deg) & (obj_altaz.az < 325*u.deg)
    
    # Make sure alt is greater than 18 degrees
    notset = obj_altaz.alt > 18*u.deg #select only times where the object is setting. 
    
    nwrap=np.logical_and(nwrap,notset)
    swrap=np.logical_and(swrap,notset)
    nswrap=np.logical_and(nswrap,notset)
    
    if len(obj_altaz[nwrap])>0:
     n_start=rtime(obj_altaz[nwrap].obstime[0].datetime)
     n_end=rtime(obj_altaz[nwrap].obstime[-1].datetime)
    if len(obj_altaz[nwrap])==0:
      n_start,n_end=None,None
      #print('Not in North Wrap')
    
    if len(obj_altaz[swrap])>0:
     s_start=rtime(obj_altaz[swrap].obstime[0].datetime)
     s_end=rtime(obj_altaz[swrap].obstime[-1].datetime)
    if len(obj_altaz[swrap])==0:
     s_start,s_end=None,None
     #print('Not in South Wrap')

    if len(obj_altaz[nswrap])>0:
     ns_start=rtime(obj_altaz[nswrap].obstime[0].datetime)
     ns_end=rtime(obj_altaz[nswrap].obstime[-1].datetime)
    if len(obj_altaz[nswrap])==0:
     #print('Not in North/South Wrap')
     ns_start,ns_end=None,None
    
    if s_start==None and s_end==None and ns_start==None and ns_end==None:
     string='N (%s-%s)'%(n_start,n_end)
        
    if n_start==None and n_end==None and ns_start==None and ns_end==None:
     string='S (%s-%s)'%(s_start,s_end)
   
    if n_start==None and n_end==None and ns_start!=None and ns_end!=None:
     string='S,NS (%s-%s,%s-%s)'%(s_start,s_end,ns_start,ns_end)
        
    if s_start==None and s_end==None and ns_start!=None and ns_end!=None:
     string='N,NS (%s-%s,%s-%s)'%(n_start,n_end,ns_start,ns_end)

    if ns_start==None and ns_end==None and n_start!=None and n_end!=None and s_start!=None and s_end!=None:
     string='N,S (%s-%s,%s-%s)'%(n_start,n_end,s_start,s_end)
        
    if ns_start!=None and ns_end!=None and n_start!=None and n_end!=None and s_start!=None and s_end!=None:
     string='N,S,N/S (%s-%s,%s-%s,%s-%s)'%(n_start,n_end,s_start,s_end,ns_start,ns_end)

  
    #return('N',n_start,n_end,'S',s_start,s_end,'N/S',ns_start,ns_end,string)
    #print(ns_start,ns_end,n_start,n_end,s_start,s_end)
    return(string)

def moon_distance(sky_coord, date_UT,degrees=15):
    
    date_UT = Time(date_UT)
    utc_offset = -10*u.hour #Hawaii standard time
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-7, 7, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)
    
    obj_altaz = sky_coord.transform_to(observer_frame)
    keck = Observer(location=keck_location, name="keck", timezone="US/Hawaii")
    moon_altaz = keck.moon_altaz(midnight+delta_midnight) 
    
    #moon=moon_altaz.alt> 0*u.degree
    moon_sep = obj_altaz.separation(moon_altaz)
    
    moon_warning= moon_altaz[moon_sep <degrees*u.deg]

    if len(moon_warning)>0:
        return("Moon distance: %.1f deg"%(np.median(moon_sep).deg))
    if len(moon_warning)==0:
        return('')


def keck_alt_warning(sky_coord, date_UT):
    
    date_UT = Time(date_UT)
    utc_offset = -10*u.hour #Hawaii standard time
    midnight = date_UT - utc_offset
    # start_time = date_UT - utc_offset - 6*u.hour #6pm
    # end_time =    date_UT - utc_offset + 6*u.hour #6am
    delta_midnight = np.linspace(-7, 7, 720*2)*u.hour #one per 0.5 minute
    observer_frame = AltAz(obstime=midnight+delta_midnight,
                              location=keck_location)
    obj_altaz = sky_coord.transform_to(observer_frame)
    alt_high = obj_altaz.alt > 85*u.deg
    
    if len(obj_altaz[alt_high])>0:
     alt_start=rtime(obj_altaz[alt_high].obstime[0].datetime)
     alt_end=rtime(obj_altaz[alt_high].obstime[-1].datetime)
     return('Unstable guiding (alt > 85 deg): %s-%s'%(alt_start,alt_end))
    else:
     alt_start,alt_end,str=None, None,None
     return('')




def NIRES_telluric_exp_time(optical_mag):
    optical_mag=float(optical_mag)
    if optical_mag<7:
        exp_time_str='4x2s (bright, keep off slit centre)'
        exp_time_ind=2
    elif 7.0 <= optical_mag< 7.5:
        exp_time_str='2 / 4x2'
        exp_time_ind=2
    elif 7.5 <= optical_mag< 8:
        exp_time_str='2 / 4x2'
        exp_time_ind=2
    elif 8 <= optical_mag< 8.5:
        exp_time_str='2 / 4x3'
        exp_time_ind=3
    elif 8.5 <= optical_mag< 9:
        exp_time_str='2 / 4x5'
        exp_time_ind=5
    elif 9 <= optical_mag< 9.5:
        exp_time_str='3 / 4x10'
        exp_time_ind=10
    elif 9.5 <= optical_mag< 10:
        exp_time_str='3 / 4x15'
        exp_time_ind=15
    elif 10 <= optical_mag< 11:
        exp_time_str='3 / 4x20'
        exp_time_ind=20
    elif 11 <= optical_mag< 12:
        exp_time_str='3 / 4x30'
        exp_time_ind=30
    elif optical_mag>= 12:
        exp_time_str='faint telluric, pick brighter one'
        exp_time_ind=60
    else:
        exp_time_str='telluric does not have V mag. Check.'
        exp_time_ind=60        
    # print(optical_mag, exp_time_str, exp_time_ind)
    return exp_time_str,exp_time_ind

def NIRES_target_exp_time(optical_mag):
    optical_mag=float(optical_mag)
    if optical_mag<13:
        exp_time_str='20 / 4x60s'
        exp_time_ind=60
        dither="ABBA"
    if 13 <= optical_mag< 14:
        exp_time_str='20 / 4x80'
        exp_time_ind=80
        dither="ABBA"
    if 14 <= optical_mag< 15:
        exp_time_str='30 / 4x100'
        exp_time_ind=100
        dither="ABBA"
    if 15 <= optical_mag< 16:
        exp_time_str='30 / 4x120'
        exp_time_ind=120
        dither="ABBA"
    if 16 <= optical_mag< 16.5:
        exp_time_str='30 / 4x150'
        exp_time_ind=150
        dither="ABBA"
    if 16.5 <= optical_mag< 17:
        exp_time_str='30 / 4x200'
        exp_time_ind=200
        dither="ABBA"
    if 17 <= optical_mag< 17.5:
        exp_time_str='30 / 4x250'
        exp_time_ind=250
        dither="ABBA"
    if 17.5 <= optical_mag< 18:
        exp_time_str='30 / 4x300'
        exp_time_ind=300
        dither="ABBA"
    if 18.0 <= optical_mag< 18.5:
        exp_time_str='40 / 6x300'
        exp_time_ind=300
        dither="ABBAAB"
    if 18.5 <= optical_mag< 19.0:
        exp_time_str='40 / 8x300'
        exp_time_ind=300
        dither="ABBAABBA"
    if optical_mag>= 19:
        exp_time_str='faint, need ~1 hour on it, do you really want to do this?'
        exp_time_ind=300
        dither="ABBAABBA"

    return exp_time_str,exp_time_ind,dither


def LRIS_exp_time(optical_mag):
    optical_mag=float(optical_mag)
    if optical_mag<13:
        exp_time_str='110, 100'
        exp_time_b=110
        exp_time_r=100
        n_r=1
        n_b=1
    if 13 <= optical_mag< 15:
        exp_time_str='120, 110'
        exp_time_b=120
        exp_time_r=110
        n_r=1
        n_b=1
    if 15 <= optical_mag< 17:
        exp_time_str='250, 240'
        exp_time_b=250
        exp_time_r=240
        n_r=1
        n_b=1
    if 17 <= optical_mag< 18:
        exp_time_str='310, 300'
        exp_time_b=310
        exp_time_r=300
        n_r=1
        n_b=1
    if 18 <= optical_mag< 19:
        exp_time_str='610, 600'
        exp_time_b=610
        exp_time_r=600
        n_r=1
        n_b=1
    if 19 <= optical_mag< 19.5:
        exp_time_str='900, 2x417'
        exp_time_b=900
        exp_time_r=417
        n_r=2
        n_b=1
    if 19.5 <= optical_mag< 20:
        exp_time_str='1000, 2x467'
        exp_time_b=1000
        exp_time_r=467
        n_r=2
        n_b=1
    if 20 <= optical_mag< 20.5:
        exp_time_str='1100, 2x517'
        exp_time_b=1100
        exp_time_r=517
        n_r=2
        n_b=1
    if 20.5 <= optical_mag< 21.5:
        exp_time_str='1200, 2x567'
        exp_time_b=1200
        exp_time_r=567
        n_r=2
        n_b=1
    if 21.5 <= optical_mag< 22.5:
        exp_time_str='1800, 3x556'
        exp_time_b=1800
        exp_time_r=556
        n_r=3
        n_b=1
    if 22.5 <= optical_mag< 23.0:
        exp_time_str='2100, 3x656'
        exp_time_b=2100
        exp_time_r=656
        n_r=3
        n_b=1
    if 23.0 <= optical_mag< 23.5:
        exp_time_str='2x1200, 4x567'
        exp_time_b=1200
        exp_time_r=567
        n_r=4
        n_b=2
    if 23.5 <= optical_mag< 24.0:
        exp_time_str='2x1400, 4x667'
        exp_time_b=1400
        exp_time_r=667
        n_r=4
        n_b=2
    if 24.0 <= optical_mag< 24.5:
        exp_time_str='2x1800, 6x556'
        exp_time_b=1800
        exp_time_r=556
        n_r=6
        n_b=2
    if 24.5 <= optical_mag< 25.0:
        exp_time_str='2x2000, 6x623'
        exp_time_b=2000
        exp_time_r=623
        n_r=6
        n_b=2
    if optical_mag>= 25:
        exp_time_str=''
        exp_time_b=''
        exp_time_r=''
        n_r=''
        n_b=''
    return exp_time_str,exp_time_b,exp_time_r,n_b,n_r

"""import sys
if len(sys.argv) < 3:
    print('Usage: python LRIS_exp_time.py individual_blue_exp_time #red_desired #blue_desired')
    print('If #blue_desired not provided, defaults to 1.\n\n')
    print('Example: python LRIS_exp_time.py 1200 2 \n--> will give you the individual exposure time for the red side if the blue side is 1x1200 sec, \nand you want to take 2 red exposures.')
else:
    blue_exp = int(sys.argv[1])
    num_red = int(sys.argv[2])
    if len(sys.argv) == 4:
        num_blue = int(sys.argv[3])
    else:
        num_blue = 1

    #Readout times from https://www2.keck.hawaii.edu/inst/lris/win_bin.html
    blue_readout = 54 #s, 1x1 binning, full frame 4 amps
    red_readout = 65 #x, Spec 2x1 2 amp

    ### We want to make sure that the last frame of blue and red finish exposing at the same time so that we can slew. 
    ###  (blue_exp + blue_readout) * num_blue - blue_readout = (red_exp + red_readout) * num_red - red_readout
    red_exp = (((blue_exp + blue_readout)*num_blue - blue_readout) + red_readout)/num_red -red_readout

    print("Appropriate exposure time for each red frame is:\n%d"%red_exp)"""



# Keck I   Azimuth range  5.3° to 146.2°   , Elevations accessible 33.3° to 88.9°
# Keck II Azimuth range  185.3° to 332.8° , Elevations accessible 36.8° to 89.5°°
# Stable guiding cannot be guaranteed for targets transiting at elevations higher than 85°.

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=\
        '''
        Calculates observing contraints for a given night
        
        Usage: observing_contraints.py target_list 
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('filename',help = 'Final LRIS or NIRES target list.')
    parser.add_argument('instrument',help = 'Instrument. Options: LRIS, NIRES')
    parser.add_argument('--utdate', dest='utdate',type = Time,help = 'fmt: yyyy-mm-dd. Date of observations. Superseeds taking the time from the filename.')
    parser.add_argument("--fhalf", action="store_true",help="first half night")
    parser.add_argument("--shalf", action="store_true",help="second half night")
    parser.add_argument("--comments", "-c", help="Original YSE PZ target list with comments", default=None)

    args = parser.parse_args()

    lat = 19.8263
    lon = -155.47441
    height= 4145 #4159.581216 # metres
    
    cmt='#'

    #If comments are given, open that file and store comments and SN name
    cmt_dict = {}
    if args.comments is not None:
        comments = open(args.comments, 'r')
        for line in comments:
            if line[0] != '!':
                sn_name = line.split(' ')[0]
                comment = (line.split('comment = ')[1]).replace("\n","")
                # print(sn_name, comment)
                if sn_name in cmt_dict.keys():
                    cmt_dict[sn_name] = cmt_dict[sn_name]+' '+comment
                else:
                    cmt_dict[sn_name] = comment
        print(cmt_dict)
    #keck = Observer.at_site("Keck")
    keck_location = EarthLocation(lat=lat*u.deg,lon=lon*u.deg,height=height*u.m)
    keck = Observer(location=keck_location, name="keck", timezone="US/Hawaii")
    utc_offset = -10*u.hour #Hawaii standard time
    filename = args.filename

    if args.utdate:
        utdate=args.utdate
        date = Time(utdate)-TimeDelta(1)
        astropy_utdate = Time(utdate)
        print("User specified UT date: ",utdate)

    else:
        if args.instrument=='NIRES':
            date = filename.split(".")[0].split("Keck_II_")[1].split("_")[0]
            utdate = Time(date)+TimeDelta(1)
        if args.instrument=='LRIS':
            date = filename.split(".")[0].split("Keck_I_")[1].split("_")[0]
            utdate = Time(date)+TimeDelta(1)
        print("Using UT date from the filename: ",utdate)


    # USNO, defintion of sunset time is when the solar disk center is at -0.8333 degrees altitude
    # to account for the solar radius and atmospheric refraction, I found -0.7 matches well with keck, but this can be changed!
    times = {'utdate':utdate,
            'caldate':date,
            'sunset': keck.sun_set_time(utdate, which='next',horizon=0*u.deg,n_grid_points=200).datetime,
            'etwi05': keck.sun_set_time(utdate, which='next',horizon=-0.5*u.deg,n_grid_points=200).datetime,
            'etwi06': keck.sun_set_time(utdate, which='next',horizon=-0.6*u.deg,n_grid_points=200).datetime,
            'etwi07': keck.sun_set_time(utdate, which='next',horizon=-0.7*u.deg,n_grid_points=200).datetime,
            'etwi4': keck.sun_set_time(utdate, which='next',horizon=-4.0*u.deg,n_grid_points=200).datetime,
            'etwi5': keck.sun_set_time(utdate, which='next',horizon=-5.0*u.deg,n_grid_points=200).datetime,
            'etwi6': keck.sun_set_time(utdate, which='next',horizon=-6.0*u.deg,n_grid_points=200).datetime,
            'etwi8': keck.sun_set_time(utdate, which='next',horizon=-8.0*u.deg,n_grid_points=200).datetime,
            'etwi10': keck.sun_set_time(utdate, which='next',horizon=-10.0*u.deg,n_grid_points=200).datetime,
            'etwi12': keck.sun_set_time(utdate, which='next',horizon=-12.0*u.deg,n_grid_points=200).datetime,
            'etwi18': keck.sun_set_time(utdate, which='next',horizon=-18.0*u.deg,n_grid_points=200).datetime,
            'midnight': keck.midnight(utdate, which='next').datetime,
            'mtwi18': keck.sun_rise_time(utdate, which='next',horizon=-18*u.deg,n_grid_points=200).datetime,
            'mtwi12': keck.sun_rise_time(utdate, which='next',horizon=-12*u.deg,n_grid_points=200).datetime,
            'mtwi10': keck.sun_rise_time(utdate, which='next',horizon=-10*u.deg,n_grid_points=200).datetime,
            'mtwi8': keck.sun_rise_time(utdate, which='next',horizon=-8*u.deg,n_grid_points=200).datetime,
            'mtwi6': keck.sun_rise_time(utdate, which='next',horizon=-6*u.deg,n_grid_points=200).datetime,
            'mtwi5': keck.sun_rise_time(utdate, which='next',horizon=-5*u.deg,n_grid_points=200).datetime,
            'mtwi4': keck.sun_rise_time(utdate, which='next',horizon=-4*u.deg,n_grid_points=200).datetime,
            'mtwi07': keck.sun_rise_time(utdate, which='next',horizon=-0.7*u.deg,n_grid_points=200).datetime,
            'mtwi06': keck.sun_rise_time(utdate, which='next',horizon=-0.6*u.deg,n_grid_points=200).datetime,
            'mtwi05': keck.sun_rise_time(utdate, which='next',horizon=-0.5*u.deg,n_grid_points=200).datetime,
            'sunrise': keck.sun_rise_time(utdate, which='next',horizon=0*u.deg,n_grid_points=200).datetime,

        }
    print('caldate  %s'%times['caldate'])
    print('utdate   %s'%times['utdate'])
    print('sunset   %s'%rtime(times['sunset']))
    print('etwi05   %s'%rtime(times['etwi05']))
    print('etwi06   %s'%rtime(times['etwi06']))
    print('etwi07   %s'%rtime(times['etwi07']))
    print('etwi4    %s'%rtime(times['etwi4']))
    print('etwi5    %s'%rtime(times['etwi5']))
    print('etwi6    %s'%rtime(times['etwi6']))
    print('etwi8    %s'%rtime(times['etwi8']))
    print('etwi10   %s'%rtime(times['etwi10']))
    print('etwi12   %s'%rtime(times['etwi12']))
    print('etwi18   %s'%rtime(times['etwi18']))
    print('midnight %s'%rtime(times['midnight']))
    print('mtwi18   %s'%rtime(times['mtwi18']))
    print('mtwi12   %s'%rtime(times['mtwi12']))
    print('mtwi10   %s'%rtime(times['mtwi10']))
    print('mtwi8    %s'%rtime(times['mtwi8']))
    print('mtwi6    %s'%rtime(times['mtwi6']))
    print('mtwi5    %s'%rtime(times['mtwi5']))
    print('mtwi4    %s'%rtime(times['mtwi4']))
    print('mtwi07   %s'%rtime(times['mtwi07']))
    print('mtwi06   %s'%rtime(times['mtwi06']))
    print('mtwi05   %s'%rtime(times['mtwi05']))
    print('sunrise  %s'%rtime(times['sunrise']))

    #print('mtwi06   %s'%rtime(times['mtwi06']))

    

    # NIRES should start/end ~20 min after/before sunset/sunrise, or 5 degree
    if args.instrument=='NIRES':
        if args.fhalf:
            times.update({'startime':times['etwi5']})
            times.update({'endtime':times['midnight']})
        if args.shalf:
            times.update({'startime':times['midnight']})
            times.update({'endtime':times['mtwi5']})
        if not (args.fhalf or args.shalf):
            print('NIRES FULL NIGHT')
            times.update({'startime':times['etwi5']})
            times.update({'endtime':times['mtwi5']})

    # LRIS can start at like 8 degree ish
    if args.instrument=='LRIS':
        if args.fhalf:
            times.update({'startime':times['etwi8']})
            times.update({'endtime':times['midnight']})
        if args.shalf:
            times.update({'startime':times['midnight']})
            times.update({'endtime':times['mtwi8']})
        if not (args.fhalf or args.shalf):
            print('NIRES FULL NIGHT')
            times.update({'startime':times['etwi8']})
            times.update({'endtime':times['mtwi8']})


    file = open(filename, 'r')
    lines = file.read().split('\n')
    out_file = open('spreadsheet_'+filename, 'w')
    print("Writing to: "+'spreadsheet_'+filename)

    if args.instrument=='NIRES':
        main_nires_header='\tName\tRa\tDec\tMag (B,V)\t Exposure seq (svc/spec)\t Dither \t Rotdest \tTime(m) w overhead\tIndividual Time(s)\tTotal Time(s)\tTotal Time (min)\tSetting time(UT)\tWraps \t Warnings\n'
        shalf_nires_header='Start(UT)\tEnd(18,12,5,0 deg) (UT)'
        fhalf_nires_header='Start(0,5,12,18 deg)(UT)\tEnd (UT)'
        full_nires_header='Start(0,5,12,18 deg)(UT)\tEnd(18,12,5,0 deg) (UT)'

        if args.shalf:
            header=shalf_nires_header+main_nires_header
            out_file.write(header)
            out_file.write(rtime(times['midnight'])+'\t'+rtime(times['mtwi18'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi5'])+','+rtime(times['mtwi07'])+'\n')
            out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
 
        if args.fhalf:
            header=fhalf_nires_header+main_nires_header
            out_file.write(header)
            out_file.write(rtime(times['etwi07'])+','+rtime(times['etwi5'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi18'])+'\t'+rtime(times['midnight'])+'\n')
            out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
 
        if not (args.fhalf or args.shalf):
            header=full_nires_header+main_nires_header
            out_file.write(header)
            out_file.write(rtime(times['etwi07'])+','+rtime(times['etwi5'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi18'])+'\t'+rtime(times['mtwi18'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi5'])+','+rtime(times['mtwi07'])+'\n')
            out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')


    if args.instrument=='LRIS':
     # need to copy lris log...
        main_lris_header='\tName\tRa\tDec\tMag (r)\t Exposure seq (b / r)\t Rotdest \tTime(m) w overhead\tTotal Time(s)\tTotal Time(min)\tRising time(UT)\tWraps\t Warnings\n'
        shalf_lris_header='Start(UT)\tEnd(18,12,8,0 deg) (UT)'
        fhalf_lris_header='Start(0,8,12,18 deg)(UT)\tEnd (UT)'
        full_lris_header='Start(0,8,12,18 deg)(UT)\tEnd(18,12,8,0 deg) (UT)'

        if args.shalf:
            header=shalf_lris_header+main_lris_header
            out_file.write(header)
            out_file.write(rtime(times['midnight'])+'\t'+rtime(times['mtwi18'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi8'])+','+rtime(times['mtwi07'])+'\n')
            out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
 
        if args.fhalf:
            header=fhalf_lris_header+main_lris_header
            out_file.write(header)
            out_file.write(rtime(times['etwi07'])+','+rtime(times['etwi8'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi18'])+'\t'+rtime(times['midnight'])+'\n')
            out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
 
        if not (args.fhalf or args.shalf):
            header=full_lris_header+main_lris_header
            out_file.write(header)
            out_file.write(rtime(times['etwi07'])+','+rtime(times['etwi8'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi18'])+'\t'+rtime(times['mtwi18'])+','+rtime(times['mtwi12'])+','+rtime(times['mtwi8'])+','+rtime(times['mtwi07'])+'\n')
            out_file.write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')


    for ind, i in enumerate(lines):

        if len(i.strip()) > 1:
            if i.strip()[0] != cmt : #The entire line is not commented
                if "_S" not in i.split()[0]: #this is a target, not an offset star
                    split=i.split()
                    #print(len(split))                    
                    name=split[0]
                    ra = split[1]+' '+split[2]+' '+split[3]
                    dec= split[4]+' '+split[5]+' '+split[6]

                    ra_sheet = "'"+split[1]+':'+split[2]+':'+split[3]
                    dec_sheet= "'"+split[4]+':'+split[5]+':'+split[6]

                    source_coords = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg))
                    
                    wraps=keck_wrap_time(source_coords,times['utdate'])
                    warning_alt=keck_alt_warning(source_coords,times['utdate'])
                    warning_moon=moon_distance(source_coords,times['utdate'],15)
                    #print(name,ra_sheet,dec_sheet,warning_alt,warning_moon)

                    #t_altaz_midnight = keck.altaz(times['midnight'], source_coords)
                    #t_altaz_start = keck.altaz(times['startime'], source_coords)
                    #t_altaz_end = keck.altaz(times['endtime'], source_coords)

                    #print('azimuth range %.0f-%.0f, dec %s'%(t_altaz_start.az.deg,t_altaz_end.az.deg,dec))
                    
                    #print(name+" K1 rising time is ", k1_rise)
                    #print(type(k2_set.datetime))
                    #print(k2_set.datetime))
                    
                    if args.instrument=='NIRES':
                        k2_set=keck2_setting_time(source_coords, times['utdate'])
                        #print("K2 setting time is ",k2_set)
                        if k2_set==None:
                            k2_set='-'
                        #print(rtime(k2_set))

                        if "HIP" in name:
                            vmag=split[17]
                            bmag=split[14]
                            dither='ABBA'
                            exp_time_str=NIRES_telluric_exp_time(vmag)[0]
                            exp_time=NIRES_telluric_exp_time(vmag)[1]
                            telluric_tot='00:10:00'
                            rotdest=float(split[9][8:])
                            #total_exp_time_s=exp_time*len(dither)
                            #total_exp_time_m=total_exp_time_s/60.
                            out_file.write('\t\t %s \t %s \t %s \t %s,%s\t%s\t%s\t%.0f\t%s \t %s \t \t \t %s\t %s \t %s \n'%\
                                (name,ra_sheet,dec_sheet,bmag,vmag,exp_time_str,dither,rotdest,telluric_tot,exp_time,k2_set,wraps,warning_alt+' '+warning_moon))
                      
                        if not "HIP" in name:
                            sn_mag=split[11].split('r=')[1]
                            exp_time_str=NIRES_target_exp_time(sn_mag)[0]
                            exp_time=NIRES_target_exp_time(sn_mag)[1]
                            dither=NIRES_target_exp_time(sn_mag)[2]
                            rotdest=float(split[9][8:])
                            total_exp_time_s=exp_time*len(dither)
                            total_exp_time_m=math.ceil(total_exp_time_s/60.)
                            time_w_overheads=total_exp_time_m+7 # increased overheads to 7 minutes per target 
                            if name in cmt_dict.keys():
                                out_file.write('\t\t%s\t%s\t%s\t%s \t %s \t %s \t %.0f \t 00:%s:00 \t %s \t %s \t %.1f \t %s \t %s \t %s \t %s \n'%\
                                    (name,ra_sheet,dec_sheet,sn_mag,exp_time_str,dither,rotdest,time_w_overheads,exp_time,total_exp_time_s,total_exp_time_m,k2_set,wraps,warning_alt+' '+warning_moon, cmt_dict[name]))
                            else:
                                out_file.write('\t\t%s\t%s\t%s\t%s \t %s \t %s \t %.0f \t 00:%s:00 \t %s \t %s \t %.1f \t %s \t %s \t %s \n'%\
                                    (name,ra_sheet,dec_sheet,sn_mag,exp_time_str,dither,rotdest,time_w_overheads,exp_time,total_exp_time_s,total_exp_time_m,k2_set,wraps,warning_alt+' '+warning_moon))


                    if args.instrument=='LRIS':
                        print(name) 
                        k1_rise=keck1_rising_time(source_coords, times['utdate'])
                        #print('K1 rising time is:',k1_rise)
                        if k1_rise==None:
                            k2_rise='-'

                        rmag=split[11]
                        rotdest=float(split[9][8:])
                        lris_et=LRIS_exp_time(rmag[2:])
                        exp_time_str=lris_et[0]
                        exp_time_b=lris_et[1]
                        exp_time_r=lris_et[2]
                        n_b=lris_et[3]
                        n_r=lris_et[4]

                        total_exp_time_s=exp_time_b*n_b
                        total_exp_time_m=math.ceil(total_exp_time_s/60.)
                        time_w_overheads=(total_exp_time_m)+5 

                        #print('\t\t%s \t%s \t%s \t%s \t%s \t%.0f\t 00:%s:00 \t%s \t%s \t%s \t%s \t%s \n'%\
                        #     (name,ra_sheet,dec_sheet,rmag[2:],exp_time_str,rotdest,time_w_overheads,total_exp_time_s,total_exp_time_m,k1_rise,wraps,warning_alt+' '+warning_moon))
                        if name in cmt_dict.keys():
                            out_file.write('\t\t%s \t%s \t%s \t%s \t%s \t %.0f \t 00:%s:00 \t%s \t%s \t%s \t%s \t%s \t %s \n'%\
                                (name,ra_sheet,dec_sheet,rmag[2:],exp_time_str,rotdest,time_w_overheads,total_exp_time_s,total_exp_time_m,k1_rise,wraps,warning_alt+' '+warning_moon, cmt_dict[name]))
                        else:
                            out_file.write('\t\t%s \t%s \t%s \t%s \t%s \t %.0f \t 00:%s:00 \t%s \t%s \t%s \t%s \t%s \n'%\
                                (name,ra_sheet,dec_sheet,rmag[2:],exp_time_str,rotdest,time_w_overheads,total_exp_time_s,total_exp_time_m,k1_rise,wraps,warning_alt+' '+warning_moon))




import numpy as np
from astropy.io import ascii as asci
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

def getTelluric(coords, max_distance = 20, spec_type = 'A0V', Vmin = 6.5, Vmax = 11):
    """
    Given a source coordinates in SkyCoord format, query the Hipparcos catalog on Vizier to get the closest A0V star. 
    Return a result table with RA, Dec, Spectral type, Vmag, distance from source, and RA difference    
    """
    #Define what columns to query and what the conditions are
    hip_viz = Vizier(columns=['HIP', 'RAhms', 'DEdms', 'Vmag', 'SpType'], \
               column_filters={"Vmag":"%s..%s"%(str(Vmin), str(Vmax)), "SpType":spec_type}, row_limit = -1)

    #query. Use twice the max distance and add flag if the closest one is too far
    result = hip_viz.query_region(coords, radius=2*max_distance*u.deg, catalog='I/239')
    # try:
    print(result)
    try:
        a0v_res = result[0] #[0] is Hipparcos; [1] is Tycho
    except:
        result = hip_viz.query_region(coords, radius=4*max_distance*u.deg, catalog='I/239')
        a0v_res = result[0]

    #Get B mag from Simbad, also check spectral type
    Simbad.add_votable_fields("flux(B)")
    Simbad.add_votable_fields("flux(V)")
    Simbad.add_votable_fields("sptype")

    add_mag = Simbad.query_objects(["HIP"+str(x) for x in a0v_res['HIP']])

    res_coords = SkyCoord(ra = a0v_res['RAhms'], dec = a0v_res['DEdms'], unit = (u.hourangle, u.deg))
    a0v_res['distance'] = coords.separation(res_coords)
    a0v_res['dRA'] = (coords.ra.deg - res_coords.ra.deg)
    a0v_res['abs_dRA'] = np.abs(coords.ra.deg - res_coords.ra.deg)
    a0v_res['Simbad_spt'] = add_mag['SP_TYPE']
    a0v_res['Simbad_B'] = add_mag['FLUX_B']
    a0v_res['Simbad_V'] = add_mag['FLUX_V']

    actually_A0V = add_mag['SP_TYPE'] == 'A0V'

    #print(a0v_res)

    # if sortby == 'dra':
    #     a0v_res.sort('dRA_(source-cal)')
    # else:
    #     a0v_res.sort('distance')
    # except:
    #     print("No nearby A0V telluric")
    #     a0v_res = None

    return a0v_res[actually_A0V]

#Now the part we execute

if __name__ == "__main__":
    format = 'keck' #or 'iobserve', or 'keck'
    #growth is name hh mm ss dd mm ss !comments
    #iobserve is name hh:mm:ss dd:mm:ss #comments
    outformat = 'keck'

    #TODO update this to argparse
    if len(sys.argv) == 1:
        file_name = input("Type file_name of coordinates ")
    elif len(sys.argv) == 2:
        file_name = sys.argv[1]
    elif len(sys.argv) == 3:
        file_name = sys.argv[1]
        format = sys.argv[2]
    elif len(sys.argv) == 4:
        file_name = sys.argv[1]
        format = sys.argv[2]
        outformat = sys.argv[3]

    if format == 'growth':
        cmt = '!'
        spc = ' '
        print("format is %s"%format)
    elif format == 'iobserve':
        cmt = '#'
        spc = ':'
        print("format is %s"%format)
    elif format == 'keck':
        cmt = '#'
        spc = ' '
    #Open the file
    file = open(file_name, 'r')
    out_file = open('with_telluric_'+file_name, 'w')

    lines = file.read().split('\n')
    #Go through the input line by line
    if outformat == 'irtf':
        irtf_counter = 0
    for ind, i in enumerate(lines):
        if len(i.strip()) > 1:
            if i.strip()[0] != cmt : #The entire line is not commented
                if "Telluric" not in (i.split(cmt)[1]): #Not already a Telluric source
                    splitted = i.split()
                    if "_S" not in i.split()[0]: #this is a target, not an offset star
                        if format == 'iobserve': # hh:mm:ss dd:mm:ss
                            source_coords  = SkyCoord(ra =splitted[1], dec = splitted[2], unit = (u.hourangle, u.deg))
                            ra_h = splitted[1].split(':')[0]
                            ra_m = splitted[1].split(':')[1]
                            ra_s = splitted[1].split(':')[2]
                            dec_d = splitted[2].split(':')[0]
                            dec_m = splitted[2].split(':')[1]
                            dec_s = splitted[2].split(':')[2]
                            # dec = splitted[4]+' '+splitted[5]+' '+splitted[6]
                            # source_coords = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg))
                            # print(source_coords)
                        elif format == 'growth' or format == 'keck': # hh mm ss  dd mm ss
                            ra_h =  splitted[1]
                            ra_m =  splitted[2]
                            ra_s =  splitted[3]
                            dec_d = splitted[4]
                            dec_m = splitted[5]
                            dec_s = splitted[6]                            
                            ra = splitted[1]+' '+splitted[2]+' '+splitted[3]
                            dec = splitted[4]+' '+splitted[5]+' '+splitted[6]
                            source_coords = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg))
                        print("Getting telluric for %s"%splitted[0])
                        tellurics = getTelluric(source_coords)
                        #Write source to output file
                        if outformat == 'irtf':
                            irtf_out_string = '{:s}     {:s}{:s}:{:s}:{:s}  {:s}:{:s}:{:s} 2000.0 0.0 0.0'.format(str(irtf_counter).zfill(2), splitted[0].ljust(14), \
                                                                                                          ra_h, ra_m, ra_s, \
                                                                                                        dec_d,dec_m,dec_s)
                            out_file.write(irtf_out_string+'\n')
                            irtf_counter+=1
                        else:
                            out_file.write(i+'\n')
                        #Check if rotdest is set. If so, make sure telluric has the same rotdest
                        if "rotdest" in i:
                            #loop through to find the component with rotdest
                            for j in splitted:
                                if "rotdest" in j:
                                    rotdest = j.split('=')[1]
                        else:
                            rotdest = None
                        #Write best telluric to output file
                        if tellurics is not None:
                            #Two Closest telluric in distance
                            ### FROM EXPERIENCE, closest telluric in RA is useless. Get 2 closest in distance. 
                            # best_tel_dist = tellurics[tellurics['distance'] == np.min(tellurics['distance'])][0]
                            tellurics.sort('distance')
                            min_RA_sep = 0.01 #arcmin
                            good_RA_diff = np.logical_and(tellurics["abs_dRA"] > min_RA_sep/60*15 , \
                                                            tellurics["abs_dRA"] < 1.2*15)
                            #ra difference bigger than 10 mins, smaller than 1.2 hour
                            # print(tellurics[good_RA_diff])
                            ###########################################TO DO, MAKE SURE YOU ALWAYS GET ONE TO EAST AND ONE TO WEST##################
                            if np.sum(good_RA_diff.astype(int)) > 2:
                                best_tel_dist = tellurics[good_RA_diff][0]
                                best_tel_dist2 = tellurics[good_RA_diff][1]
                            else:
                                print("All telluric are closer than %d arcmin in RA or further than 1 hour. Check results"%min_RA_sep)
                                best_tel_dist = tellurics[0]
                                best_tel_dist2 = tellurics[1]      
                            # tellurics.sort('abs_dRA')
                            # best_tel_ra = tellurics[0]
                            # print(best_tel_ra)
                            if outformat in ['iobserve', 'growth', 'keck']:
                                if rotdest is None:
                                    outstring1 = ('HIP'+str(best_tel_dist['HIP'])).ljust(16)+best_tel_dist['RAhms'].replace(' ',spc)+' '+best_tel_dist['DEdms'].replace(' ',spc).ljust(12)+\
                                                    '  2000.0  '+cmt+' Telluric B = %.2f V = %.2f distance = %.2f deg, dRA = %.2f deg SpT = %s \n'%(best_tel_dist['Simbad_B'], best_tel_dist['Simbad_V'],
                                                     best_tel_dist['distance'], best_tel_dist['dRA'], best_tel_dist['Simbad_spt'])
                                    # outstring2 = ('HIP'+str(best_tel_ra['HIP'])).ljust(20)+best_tel_ra['RAhms'].replace(' ',spc)+'  '+best_tel_ra['DEdms'].replace(' ',spc)+\
                                    #                 '  2000.0  '+cmt+' Telluric V = %.2f distance = %.2f deg, dRA = %.2f deg \n'%(best_tel_ra['Vmag'], best_tel_ra['distance'], best_tel_ra['dRA'])   
                                    outstring2 = ('HIP'+str(best_tel_dist2['HIP'])).ljust(16)+best_tel_dist2['RAhms'].replace(' ',spc)+' '+best_tel_dist2['DEdms'].replace(' ',spc).ljust(12)+\
                                                    '  2000.0  '+cmt+' Telluric B = %.2f V = %.2f distance = %.2f deg, dRA = %.2f deg SpT = %s \n'%(best_tel_dist2['Simbad_B'], best_tel_dist2['Simbad_V'], 
                                                    best_tel_dist2['distance'], best_tel_dist2['dRA'], best_tel_dist2['Simbad_spt'])  
                                else:
                                    outstring1 = ('HIP'+str(best_tel_dist['HIP'])).ljust(16)+best_tel_dist['RAhms'].replace(' ',spc)+' '+best_tel_dist['DEdms'].replace(' ',spc).ljust(12)+\
                                                    '  2000.0  rotmode=pa rotdest=%s '%(rotdest)+cmt+' Telluric B = %.2f V = %.2f distance = %.2f deg, dRA = %.2f deg SpT = %s \n'\
                                                    %(best_tel_dist['Simbad_B'], best_tel_dist['Simbad_V'], best_tel_dist['distance'], best_tel_dist['dRA'], best_tel_dist['Simbad_spt'])   
                                    outstring2 = ('HIP'+str(best_tel_dist2['HIP'])).ljust(16)+best_tel_dist2['RAhms'].replace(' ',spc)+' '+best_tel_dist2['DEdms'].replace(' ',spc).ljust(12)+\
                                                    '  2000.0  rotmode=pa rotdest=%s '%(rotdest)+cmt+' Telluric B = %.2f V = %.2f distance = %.2f deg, dRA = %.2f deg SpT = %s \n'\
                                                    %(best_tel_dist2['Simbad_B'], best_tel_dist2['Simbad_V'], best_tel_dist2['distance'], best_tel_dist2['dRA'], best_tel_dist2['Simbad_spt'])                                      
                            elif outformat == 'irtf':
                                outstring1 = str(irtf_counter).zfill(2)+'     '+('HIP'+str(best_tel_dist['HIP'])).ljust(14)+best_tel_dist['RAhms'].replace(' ',':')+'  '+best_tel_dist['DEdms'].replace(' ',':').ljust(12)+\
                                                ' 2000.0 0.0 0.0\n'
                                irtf_counter+=1
                                outstring2 = str(irtf_counter).zfill(2)+'     '+('HIP'+str(best_tel_dist2['HIP'])).ljust(14)+best_tel_dist2['RAhms'].replace(' ',':')+'  '+best_tel_dist2['DEdms'].replace(' ',':').ljust(12)+\
                                                ' 2000.0 0.0 0.0\n'     
                                irtf_counter+=1
                            out_file.write(outstring1)   
                            out_file.write(outstring2)                                 
                        else:
                            outstring = cmt+' No telluric found for the source above'
                            out_file.write(outstring+'\n')
                    else: #An offset star, just add it back
                        if outformat == 'irtf':
                            irtf_out_string = '{:s}     {:s}{:s}:{:s}:{:s}  {:s}:{:s}:{:s} 2000.0 0.0 0.0'.format(str(irtf_counter).zfill(2), splitted[0].ljust(14), \
                                                                                                          splitted[1],splitted[2],splitted[3], \
                                                                                                        splitted[4],splitted[5],splitted[6])
                            irtf_counter += 1
                            out_file.write(irtf_out_string+'\n')
                        else:
                            out_file.write(i+'\n')
            else: #for commented lines, just add it back
                out_file.write(i+'\n')
    out_file.close()
    file.close()
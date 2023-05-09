# Create SOAR Observing Night Webpages
#
# To run:
# 
#   1. Create/edit user_inputs.ini file...
#   2. cd to GitHub/SOAR_Goodman_ObsWebPage/python.
#   3. python SOARObsWebpage.py <path to user_inputs.ini file>
#      e.g.:
#      python SOARObsWebpage.py ~/Desktop/user_inputs.ini


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Initial Setup
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import pandas as pd
from scipy import interpolate
import glob
import math
import os
import sys
from shutil import copy2, copyfile
import datetime

from matplotlib import path
from matplotlib.patches import Arrow
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.table import Column

from astropy.visualization import SqrtStretch
from astropy.visualization import LinearStretch
from astropy.visualization import LogStretch

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import ZScaleInterval
from astropy.visualization import MinMaxInterval

from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs.utils import skycoord_to_pixel

from astropy.nddata import Cutout2D

import prepare_obs_run

def add_options():
    import argparse
    parser = argparse.ArgumentParser(description=\
    '''
    Creates webpage and finders for SOAR observing run.

    Usage: SOARObsWebpage.py

    ''', formatter_class=argparse.RawTextHelpFormatter)

    # parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help="filename of the target list")

    parser.add_argument("--date", type=str, default='',
                        help="Date to use for creating webpage.")
    parser.add_argument("--outdir", type=str, default='.',
                        help="Output date to serialize webpage.")
    parser.add_argument("--priTargetName", type=str, default='GW',
                        help="Primary target name.")

    args = parser.parse_args()

    return(args)


args = add_options()
if args.date:
    t = Time(args.date)
else:
    t = Time(datetime.datetime.utcnow())

nite=t.datetime.strftime('%Y-%m-%d')
YYYY=t.datetime.strftime('%Y')
MM=t.datetime.strftime('%m')
DD=t.datetime.strftime('%d')

telCalYear=YYYY
telCalSemester = 'B'
if str(MM) in ['02','03','04','05','06','07']:
    telCalSemester = 'A'
if str(MM) in ['01']:
    telCalYear = str(int(telCalYear)-1)

telCalCode=telCalYear+telCalSemester

priTargetName = args.priTargetName

# Create an outputDateLine for the footer of each new webpage created today...
today = datetime.date.today()
todayString = today.strftime("%d %B %Y")
outputDateLine = """<EM>Last updated: %s </EM><BR>\n""" % todayString

# Create directories, if not already existing...
mainDir = args.outdir
niteDir=os.path.join(mainDir, nite)
print(mainDir)
print(f'Night directory is {niteDir}')
if not os.path.exists(niteDir):
    os.makedirs(niteDir)

targetsDir=os.path.join(niteDir,'Targets')
if not os.path.exists(targetsDir):
    os.makedirs(targetsDir)
    
stdStarsDir=os.path.join(targetsDir,'StdStars')
if not os.path.exists(stdStarsDir):
    os.makedirs(stdStarsDir)

priTargetsDir=os.path.join(targetsDir,priTargetName)
if not os.path.exists(priTargetsDir):
    os.makedirs(priTargetsDir)

priTargetsDataDir=os.path.join(priTargetsDir,'data')
if not os.path.exists(priTargetsDataDir):
    os.makedirs(priTargetsDataDir)

findersDir=os.path.join(targetsDir,args.priTargetName,'finders')
if not os.path.exists(findersDir):
    os.makedirs(findersDir)


# Copy some static (non-changing) files to new nite directory...

copyList = ['data/SOAR/AcqCam_Cheat_Sheet_Dec2015.pdf',    
            'data/SOAR/AcqCam_User_Guide_Dec2015.pdf',
            'data/SOAR/Goodman_Tutorial_2017.pdf',
            'data/SOAR/OldSchoolTargetAcqForm.pdf', 
            'data/SOAR/soar_quick_guide.pdf',
            'data/SOAR/goodman_hgar.dat', 
            'data/SOAR/goodman_hgarne.dat',
            'data/SOAR/Observing_Plan.html',
            'data/SOAR/Preferred_Instrument_Setup.html',
            'data/SOAR/round13-poly.csv'
           ]

for fname in copyList:
    copy2(fname, niteDir)

# Copy the static (non-changing) StdStar.html file to the new
#  standard star targets directory...

fname = 'data/SOAR/StdStars.html'
copy2(fname, stdStarsDir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 4. Create the main webpage for this nite
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

mainPage=os.path.join(niteDir,'index.html')

fout = open(mainPage,'w')

fout.write("""<HTML>\n""")
fout.write("""<HEAD>\n""")
fout.write("""<TITLE>SOAR-4m GW Target of Opportunity (ToO) Nite:  """+nite+"""</TITLE>\n""")
fout.write("""</HEAD\n""")
fout.write("""\n""")
fout.write("""<BODY>\n""")
fout.write("""\n""")
fout.write("""<H1>SOAR-4m GW Target of Opportunity (ToO) Nite:  """+nite+"""</H1>\n""")
fout.write("""\n""")
fout.write("""<b><i>SOAR ToO Policy Overview ( <A HREF=http://www.ctio.noirlab.edu/soar/content/targets-opportunity-overview>webpage</A> )""")
fout.write("""     and Specific ToO Protocols for SOAR 2018B ( <A HREF=./SOAR_TOO_Protocols.pdf>PDF</A> ) \n""")
fout.write("""</i></b>\n""")
fout.write("""\n""")
fout.write("""<P>\n""")
fout.write("""<H2> Communication</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""  <LI> <A HREF=https://darkenergysurvey.slack.com/messages/C1C5JCNLE/>DES GW Slack Channel (URL)</A>\n""")
fout.write("""  <LI> <A HREF=https://wiki.gw-astronomy.org/OpenLVEM/>OpenLVEM Wiki</A>\n""")
fout.write(f"""  <LI> <A HREF=https://time-allocation.noirlab.edu/#/telescopes/classic-schedules/?year={YYYY}&month={MM}&semesterCode={telCalCode}&telescopeId=2>SOAR TelOps Calendar</A>\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""<P>\n""")
fout.write("""<b><i></i></b>\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""<H2>Observing Strategy</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""<LI> <A HREF=./Preferred_Instrument_Setup.html> Preferred Instrument Setup</A>\n""")
fout.write("""<LI> <A HREF=./Observing_Plan.html>  Nightly Observing Plan</A>\n""")
fout.write("""<LI> <A HREF=http://www.briancasey.org/artifacts/astro/almanac.cgi?command=display&year="""+YYYY+"""&month="""+MM+"""&day="""+DD+"""&observatory=t&obslong=&obslat=&actual=&effective=&sitename=&zone=&zonename=&zoneabbrev=&daylight=>(CTIO) Sky Almanac for """+nite+"""</A>\n""")
fout.write("""<LI> <A HREF=./Targets/Targets.html>Targets</A>:\n""")
fout.write("""	<UL>\n""")
fout.write("""	  <LI> <A HREF=./Targets/"""+priTargetName+"""/"""+priTargetName+""".html>"""+priTargetName+"""</A>\n""")
fout.write("""	  <LI> <A HREF=./Targets/StdStars/StdStars.html> Standard Stars</A>\n""")
fout.write("""	</UL>\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""<H2>Guides</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""<LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/observing-soar>Observing with SOAR</A>\n""")
fout.write("""<LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/soar-astronomersobservers>SOAR Instrumentation</A>\n""")
fout.write("""<LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/goodman-high-throughput-spectrograph>Goodman High Throughput Spectrograph</A> \n""")
fout.write("""  <UL>\n""")
fout.write("""  <LI> <A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/GOODMAN/Images/New_Goodman_Cheat_Sheet.pdf>Goodman Spectrograph Cheat Sheet</A>  \n""")
fout.write("""  <LI> Step-by-Step User's Guide (2017) \n""")
fout.write("""    (<A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/GOODMAN/Goodman_Tutorial_2017.pdf>web version</A>, \n""")
fout.write("""     <A HREF=./Goodman_Tutorial_2017.pdf>local version</A>) \n""")
fout.write("""  <LI> <A HREF=./soar_quick_guide.pdf>Melissa Butner's Quick Guide</A>\n""")
fout.write("""  <LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/goodman-long-slits> Available long slits</A>\n""")
fout.write("""  </UL>\n""")
fout.write("""<LI> Goodman Acquistion Camera (GACAM) Users's Manual\n""")
fout.write("""     (<A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/GOODMAN/AcqCam_User_Guide_Dec2015.pdf>web version</A>, \n""")
fout.write("""     <A HREF=./AcqCam_User_Guide_Dec2015.pdf>local version</A>)\n""")
fout.write("""   <UL>\n""")
fout.write("""     <LI> Goodman Acquistion Camera (GACAM) Cheat Sheet\n""")
fout.write("""          (<A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/GOODMAN/AcqCam_Cheat_Sheet_Dec2015.pdf>web version</A>, \n""")
fout.write("""           <A HREF=./AcqCam_Cheat_Sheet_Dec2015.pdf>local version</A>) \n""")
fout.write("""   </UL>\n""")
fout.write("""<LI> <A HREF=https://docs.google.com/document/d/1GT4qqnBWrycDa3hsiRCXYAEHA0RdDQjJ-oFvDxPOJec/edit#heading=h.cfhdub4lvrzk>Remote Observing from FNAL (DECam-specific)</A> \n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""<H2>Sky Conditions</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""  <LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/weather-sky-monitoring-tools>Pachon and Tololo Weather, Sky Conditions, and Monitoring Tools</A>\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""<H2>Utilities</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""<LI> <A HREF=http://skyview.gsfc.nasa.gov/>SkyView</A>\n""")
fout.write("""<LI> <A HREF=http://aps.umn.edu/catalog/finder_chart/> MAPS finder chart engine</A>\n""")
fout.write("""<LI> <A HREF=http://deslabs.ncsa.illinois.edu/>DES Labs</A>\n""")
fout.write("""	<UL> \n""")
fout.write("""	  <LI> <A HREF=http://descut.cosmology.illinois.edu/>Cutout Server</A>\n""")
fout.write("""	  <LI> <A HREF=http://descut.cosmology.illinois.edu/footprint>Footprint Server</A>\n""")
fout.write("""	</UL>\n""")
fout.write("""\n""")
fout.write("""<LI> <A HREF=http://www.briancasey.org/artifacts/astro/> Brian Casey's Astronomy Tools:</A>\n""")
fout.write("""	<UL> \n""")
fout.write("""	  <LI> <A HREF=http://www.briancasey.org/artifacts/astro/airmass.cgi> Hourly Airmass Tables</A>\n""")
fout.write("""	  <LI> <A HREF=http://www.briancasey.org/artifacts/astro/skycalendar.cgi> Monthly Sky Calendar</A>\n""")
fout.write("""	  <LI> <A HREF=http://www.briancasey.org/artifacts/astro/observability.cgi> Object Seasonal Observability</A>\n""")
fout.write("""	</UL>\n""")
fout.write("""<LI> <A HREF=http://www.fourmilab.ch/cgi-bin/Yoursky?z=1&lat=30.165&ns=South&lon=70.815&ew=West> Sky Chart for CTIO</A>\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""<H2>Data Processing Material</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""\n""")
fout.write("""  <LI> <A HREF=https://github.com/DouglasLeeTucker/SOAR_Goodman_QuickReduce>SOAR Goodman Quick Reduce</A> <b><i>(Prototype for quick on-the-fly reductions.)</i></b>\n""")
fout.write("""  \n""")
fout.write("""  <LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/goodman-data-reduction-pipeline>The Goodman Data Reduction Pipeline</A> <b><i>(NEW -- but not currently useable for our 400_WD configuration.)</i></b>\n""")
fout.write("""  \n""")
fout.write("""  <LI> Comparison Arc Spectra\n""")
fout.write("""    <UL>\n""")
fout.write("""      <LI> HgAr: <A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/Instrument_Plots/400m1_HgAr_3000-7000.pdf>plot (3000A-7000A)</A>, <A HREF=./goodman_hgar.dat>linelists</A>\n""")
fout.write("""      <LI> HgAr+Ne: <A HREF=./goodman_wd400_hgarne.png>plot (3750A-7150A)</A>, <A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/GOODMAN/HgArNe_400M2_GG455_full.pdf>plot (5000A-9000A)</A>, <A HREF=./goodman_hgarne.dat>linelists</A>\n""")
fout.write("""    </UL>\n""")
fout.write("""    \n""")
fout.write("""  \n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""<H2>Forms</H2>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""  <LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/instrument-setup-forms-and-observing-reports>\n""")
fout.write("""	  General SOAR Instrument Setup Forms and Observing Reports</A>\n""")
fout.write("""     <UL>\n""")
fout.write("""     <LI> <A HREF=http://www.ctio.noirlab.edu/SOAR/Forms/INST/setup.php>Instrument Setup Form</A>\n""")
fout.write("""     <LI> <A HREF=http://www.ctio.noirlab.edu/SOAR/Forms/EON/Form.php?telescope=SOAR> SOAR End-of-night form</A>\n""")
fout.write("""     <LI> <A HREF=http://www.ctio.noirlab.edu/new/Tools/Forms/EOR/Form.php?telescope=SOAR> SOAR End-of-run form</A>\n""")
fout.write("""     </UL>\n""")
fout.write("""<LI> Goodman Paper Log\n""")
fout.write("""     (<A HREF=http://www.ctio.noirlab.edu/soar/sites/default/files/documents/GOODMAN-log.pdf>web version</A>, \n""")
fout.write("""      <A HREF=./GOODMAN-log.pdf>local version</A>)\n""")
fout.write("""<LI> <A HREF=./OldSchoolTargetAcqForm.pdf>Faint-Object Target Acquisition Checklist</A> (print out a copy for each target fainter than i=19)\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""<HR>\n""")
fout.write(outputDateLine)
fout.write("""\n""")
fout.write("""</BODY>\n""")
fout.write("""</HTML>\n""")
fout.write("""\n""")

fout.close()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 5. Create the main targets directory for this nite
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

targetsPage=os.path.join(targetsDir,'Targets.html')

fout = open(targetsPage,'w')

fout.write("""<HTML>\n""")
fout.write("""<HEAD>\n""")
fout.write("""<TITLE>SOAR Observing Run:  Targets</TITLE>\n""")
fout.write("""</HEAD>\n""")
fout.write("""\n""")
fout.write("""<BODY>\n""")
fout.write("""\n""")
fout.write("""<H2> Targets</H2>\n""")
fout.write("""\n""")
fout.write("""<HR>  \n""")
fout.write("""\n""")
fout.write("""<H3>"""+priTargetName+""" Targets</H3>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""  <LI> <A HREF=./"""+priTargetName+"""/"""+priTargetName+""".html>Main Page</A>\n""")
fout.write("""  <LI> <A HREF=./"""+priTargetName+"""/"""+priTargetName+"""targets.obs_strategy.html>Finding Charts, Estimated Exposure Times, Airmass Tables, etc.)</A>\n""")
fout.write("""  <LI> <A HREF=./"""+priTargetName+"""/"""+priTargetName+"""targets.goodman_targetlist.txt>SOAR TCS-formatted Target List</A>\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""<HR>\n""")
fout.write("""\n""")
fout.write("""<H3>Standard Stars</H3>\n""")
fout.write("""\n""")
fout.write("""<UL>\n""")
fout.write("""  <LI> <A HREF=./StdStars/StdStars.html>Main Page</A>\n""")
fout.write("""  <LI> <A HREF=http://www.ctio.noirlab.edu/soar/content/hamuy-spectrophotometric-standards>Hamuy Spectrophotometric Standards</A>\n""")
fout.write("""</UL>\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""<HR>\n""")
fout.write("""\n""")
fout.write("""<A HREF=../index.html>Return to Main Page</A>\n""")
fout.write("""\n""")
fout.write("""<HR>\n""")
fout.write(outputDateLine)
fout.write("""<BR>\n""")
fout.write("""\n""")
fout.write("""</BODY>\n""")
fout.write("""</HTML>\n""")
fout.write("""\n""")

fout.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 6. Create the various figures, finding charts, and webpages for the primary
#    targets for this nite.
#
#    Originally based on createGWTargetListHTML.py from the Feb 2019 SOAR-4m
#    observing run.
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# DES footprint coordinates...
#  Read in round13-poly.csv file...
polygonDataDF = pd.read_csv('data/SOAR/round13-poly.csv')
ra_des = coord.Angle(polygonDataDF['radeg'], unit=u.degree)
dec_des = coord.Angle(polygonDataDF['decdeg'], unit=u.degree)

# Read in GW candidate list file...
table = ascii.read(args.filename)

## For generic SOAR targets for which we want S/N~30:
##  total exptime = 10.**(0.400*i-4.133)
## (See Magellan Jan 2019 Run python script).
# For the GW program, we want a single exposure to be roughly:
#   900sec for i<=18.5      (this is probably a bit too much)
#  1200sec for 18.5<i=<19.5 (this is probably a bit too much)
#  1800sec for 19.5<i<20.5  (this is probably not enough, especially on the fainter end)
# This should be improved, especially for fainter objects not much above sky background.
# It may be hard to use SOAR for objects fainter than i=20.0, although we got not unreasonable
#  S/N in January 2019 for a practice target at G=19.9 in 1800 seconds. 
col = Column(np.power(10.,0.400*table['mag']-4.133)/3., name='exptime')
table.add_column(col)

# For particularly bright objects, no need to go below 30 sec for an individual exposure.
mask = table['mag'] < 15.2
if len(table[mask])>0:
    table[mask]['exptime_ind'] = 30.

# Plot sky distribution of each set of candidates...
# GW coordinates...
coords = SkyCoord(table['ra'], table['dec'], unit=(u.hour, u.deg))
ra_gw = coords.ra.degree
dec_gw = coords.dec.degree

# Create ds9 region file containing regions for all the candidates...
regFile = 'GWtargets.'+nite+'.reg'
freg_out = open(os.path.join(priTargetsDataDir,regFile),'w')
# Write header infor into ds9 region file...
freg_out.write("""# Region file format: DS9 version 4.1\n""")
freg_out.write("""global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n""")
freg_out.write("""fk5\n""")


# Open priTargetsPage (fout)...
priTargetsPage=os.path.join(priTargetsDir, priTargetName+'.html')
fout = open(priTargetsPage,'w')

#Write HTML header/first part of main body into priTargetsPage (fout)
fout.write("""<HTML>\n""")
fout.write("""<HEAD>\n""")
fout.write("""<TITLE>SOAR Observing Run:  """+priTargetName+""" Targets</TITLE>\n""")
fout.write("""</HEAD>\n""")
fout.write("""\n""")
fout.write("""<BODY>\n""")
fout.write("""\n""")
fout.write("""<H2>"""+priTargetName+""" Targets</H2>\n""")
fout.write("""\n""")
fout.write("""<H3> Target Lists:</H3>\n""")
fout.write("""<UL>\n""")
fout.write("""<LI> <A HREF=./"""+priTargetName+"""targets.obs_strategy.html> Formatted for Observing Strategy</A> (includes Estimated Exposure Times, Finding Charts, Airmass Tables, etc.)\n""")
fout.write("""<LI> <A HREF=./"""+priTargetName+"""targets.goodman_targetlist.txt> Formatted for the SOAR-4m TCS</A>\n""")
fout.write("""</UL>\n""")
fout.write("""<H3> <H3>Mosaic of Target Images (from DECam or SkyView): </H3>\n""")
fout.write("""\n""")
#fout.write("""<table border cellspacing=0 cellpadding=0>\n""")
#fout.write("""<table border\n""")
fout.write("""<table border cellspacing=0 cellpadding=20>\n""")

# Open GWtargets.obs_strategy.html (fout2)...
outputFile2 = os.path.join(priTargetsDir,'GWtargets.obs_strategy.html')
fout2 = open(outputFile2,'w')

#Write HTML Header into GWtargets.obs_strategy.html
fout2.write('<HTML>\n')
fout2.write('<HEAD>\n')
fout2.write('<TITLE>GW Candidates Target List</TITLE>\n')
fout2.write('</HEAD>\n')
fout2.write('\n')
fout2.write('<BODY>\n')
fout2.write('\n')
fout2.write('<H2>GW Candidates Target List</H2>\n')
fout2.write('\n')
fout2.write('<P>\n')
fout2.write('<i>Note 1:  SkyView finding charts are 10arcmin on a side (in imaging mode, the Goodmann Spectrograph has a FOV of 7.2arcmin)</i>\n')
fout2.write('<P>\n')
fout2.write('<i>Note 2:  Airmass table is for the night of '+nite+'.</i>\n')
fout2.write('<P>\n')
fout2.write('<i>Note 3:  xephem data file can be found <A HREF=./GWtargets.xephem.edb>here</A>.  (Download and copy into your ~/.xephem directory.)</i>\n')
fout2.write('<P>\n')
#fout2.write('<i>Note 4:  Finding chart FITS files, a global ds9 region file for these FITS files, and corresponding JPG images can be found in des30:/opt/des-cal/calib/DESCalibObsRuns/SOAR-Jan2019/Targets/GW/data/FITS on the FNAL DES cluster.  (Links to the FITS files and the global ds9 region file can be found in the Comments column of the target list below.)</i>\n')
#fout2.write('<P>\n')
fout2.write('\n')
fout2.write('\n')
fout2.write('<PRE>\n')

outputFile3 = os.path.join(priTargetsDir,'GWtargets.goodman_targetlist.txt')
fout3 = open(outputFile3,'w')
fout3.write('#OBJECT_ID RAHMS DECDMS EPOCH\n')

outputFile3a = os.path.join(priTargetsDir,'GWtargets.des_targetviewer_list.csv')
fout3a = open(outputFile3a,'w')

outputFile4 = os.path.join(priTargetsDir,'GWtargets.xephem.edb')
fout4 = open(outputFile4,'w')


# Iterate through GW "df_gw", outputting the target list (fout2,fout3,fout3a,fout4) 
#  and image mosaic (fout)...
#fout2.write('\n')
fout2.write('<font color="blue"><b>========================================</b></font>\n')
fout2.write('<font color="blue"><b>GW Candidates Target List</b></font>\n')
fout2.write('<font color="blue"><b>========================================</b></font>\n')
#fout2.write('\n')
fout2.write('\n')
fout2.write('ID#    Name                         RA(2000)       Dec(2000)         mag     exptime(sec)                          Finding Chart Files                                      Airmass Table    Comments\n')
fout2.write('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
fout2.write('\n')

nrow = len(table)

outimages = prepare_obs_run.do_obsrun(args.filename, 'SOAR')
finderimages = []
for file in outimages:
    basename = os.path.basename(file)
    newfile = os.path.join(findersDir, basename)
    copyfile(file, newfile)
    finderimages.append(basename)

for index, row in enumerate(table):
    
    name_gw = row['name']

    coord_hmsdms = coords[index].to_string(style='hmsdms',sep=':',precision=2)

    ra_gw_hms = coord_hmsdms.split()[0]
    dec_gw_dms = coord_hmsdms.split()[1]

    ra_gw_hms_list = ra_gw_hms.split(':')
    ra_gw_hh = str(ra_gw_hms_list[0])
    ra_gw_mm = str(ra_gw_hms_list[1])
    ra_gw_ss = str(ra_gw_hms_list[2])

    dec_gw_dms_list = dec_gw_dms.split(':')
    dec_gw_dd = str(dec_gw_dms_list[0])
    dec_gw_mm = str(dec_gw_dms_list[1])
    dec_gw_ss = str(dec_gw_dms_list[2])
    
    ra_gw = coords[index].ra.degree
    dec_gw = coords[index].dec.degree

    mag_gw = row['mag']
    
    totExpTime = row['exptime']       


    # If this is an odd-numbered row, start a new table row for the mosaic of target images...
    if (index%2) == 1:
        fout.write("""<tr>\n""")
    
    # Produced entry for this target in the mosaic of target images...
    fout.write("""<td> <img src=finders/%s>\n""" % (finderimages[index]))
    fout.write("""<P>\n""")
    fout.write("""<center>\n""")
    fout.write("""<table border>\n""")
    fout.write("""  <tr>\n""")
    fout.write("""  <th>RA</th>\n""")
    fout.write("""  <th>DEC</th>\n""")
    fout.write("""  <th>mag</th>\n""")
    fout.write("""  <th>ExpTime (sec)</th>\n""")
    fout.write("""  <th>Airmass Table</th>\n""")
    fout.write("""  </tr>\n""")
    fout.write("""  <tr>\n""")
    fout.write("""  <td>%s</td>\n""" % (ra_gw_hms))
    fout.write("""  <td>%s</td>\n""" % (dec_gw_dms))
    fout.write("""  <td>%5.2f </td>\n""" % (mag_gw))
    fout.write("""  <td>%d </td>\n""" % (int(totExpTime)))
    fout.write("""  <td><a href=http://www.briancasey.org/artifacts/astro/airmass.cgi?observatory=t&obslong=&obslat=&actual=&effective=&site=&zone=&zonename=&zoneabbrev=&daylight=1&year=%s&month=%s&day=%s&starname=%s&rightasc=%s+%s+%s&declin=%s+%s+%s&epoch=2000&altitude=-18&command=display>%s-%s-%s</a>\n""" % (YYYY, MM, DD, name_gw, ra_gw_hh, ra_gw_mm, ra_gw_ss, dec_gw_dd, dec_gw_mm, dec_gw_ss, YYYY, MM, DD))
    fout.write("""  </tr>\n""")
    fout.write("""</table>\n""")
    fout.write("""</center>\n""")
    fout.write("""<P>\n""")
    fout.write("""<center>\n""")
    fout.write("""<table border>\n""")
    fout.write("""  <tr>\n""")
    fout.write("""  <th>5'x5' DECam cutout</th>\n""")
    fout.write("""  <th>Full DECam CCD image</th>\n""")
    fout.write("""  <th>ds9 region file</th>\n""")
    fout.write("""  <th>10'x10' DSS2 red image</th>\n""")
    fout.write("""  </tr>\n""")
    fout.write("""  <tr>\n""")
    fout.write("""  <td><a href=finders/%s>local png</a></td>\n""" % (finderimages[index]))
    fout.write("""  <td><a href=http://skyview.gsfc.nasa.gov/cgi-bin/runquery.pl?Survey=DSS2+Red&GRID=J2000&gridlabels=1&position=%s+%s+%s,%s+%s+%s&Size=0.166667&Pixels=600&Return>SkyView</a>\n""" % (ra_gw_hh, ra_gw_mm, ra_gw_ss, dec_gw_dd, dec_gw_mm, dec_gw_ss) )
    fout.write("""  <P>\n""")
    fout.write("""  </tr>\n""")
    fout.write("""</table>\n""")
    fout.write("""</center>\n""")
    fout.write("""<P>\n""")
    fout.write("""\n""")
    fout.write("""</td>\n""")
    
    # If this is an even-numbered row or is the last row, end the table row for the mosaic of target images...
    if (index%2 == 0) or (index == nrow):
        fout.write("""</tr>\n""")
             

#outputLine2 = '       ----------------------------------------------------------------------------------------------------------------------\n'
outputLine2 = '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'
outputLine2 = outputLine2+'\n'
fout2.write(outputLine2)


########################################################################

# Write out final part of main body of fout...
fout.write("""</table>\n""")
fout.write("""\n""")
fout.write("""<H3>Sky Distribution of Targets</H3>\n""")
fout.write("""\n""")
fout.write("""<P>\n""")
fout.write("""<IMG SRC=./"""+priTargetName+"""targets.radec.png></CENTER>\n""")
fout.write("""<P>\n""")
fout.write("""\n""")
fout.write("""\n""")

# Info for footer of fout and fout2...
today = datetime.date.today()
todayString = today.strftime("%d %B %Y")
outputDateLine = """<EM>Last updated: %s </EM><BR>""" % todayString

#Write HTML Footer (fout)
fout.write("""<HR>\n""")
fout.write("""\n""")
fout.write("""<A HREF=../../index.html>Return to Main Page</A>\n""")
fout.write("""\n""")
fout.write("""\n""")
fout.write("""<HR>\n""")
fout.write(outputDateLine)
fout.write("""\n""")
fout.write("""</BODY>\n""")
fout.write("""</HTML>\n""")

#Write HTML Footer (fout2)
fout2.write('')
fout2.write('</PRE>')
fout2.write('<HR>')
fout2.write('<A HREF=../../index.html>Return to Main Page</A>')
fout2.write('<HR>')
fout2.write(outputDateLine)
fout2.write('')
fout2.write('</BODY>')
fout2.write('</HTML>')

#Clean up...
fout.close()
fout2.close()
fout3.close()
fout3.close()
fout3a.close()
fout4.close()
freg_out.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


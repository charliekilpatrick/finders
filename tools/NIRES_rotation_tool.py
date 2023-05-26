import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import ascii as asci

from matplotlib.patches import Rectangle, Circle
from matplotlib import transforms
from matplotlib.widgets import Slider, Button

import sys, glob, argparse

from PIL import Image

#############Define the matplotlib slider tool
#Vertical Slider bar from https://stackoverflow.com/questions/25934279/add-a-vertical-slider-with-matplotlib
#Update when this is supported by matplotlib

from matplotlib.widgets import AxesWidget
import six

class VertSlider(AxesWidget):
    """
    A slider representing a floating point range.

    For the slider to remain responsive you must maintain a
    reference to it.

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *hline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, valstep = None, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*.

        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...).

        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in

        label : str
            Slider label

        valmin : float
            The minimum value of the slider

        valmax : float
            The maximum value of the slider

        valinit : float
            The slider initial position

        label : str
            The slider label

        valfmt : str
            Used to format the slider value, fprint format string

        closedmin : bool
            Indicate whether the slider interval is closed on the bottom

        closedmax : bool
            Indicate whether the slider interval is closed on the top

        slidermin : Slider or None
            Do not allow the current slider to have a value less than
            `slidermin`

        slidermax : Slider or None
            Do not allow the current slider to have a value greater than
            `slidermax`


        dragging : bool
            if the slider can be dragged by the mouse
        valstep : float, optional, default: None
            If given, the slider will snap to multiples of `valstep`.
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.valstep = valstep
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.hline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center', 
                             fontsize = 18)

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center',
                               fontsize = 18)

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if self.valstep:
            val = np.round((val - self.valmin)/self.valstep)*self.valstep
            val += self.valmin
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
#         print('foo')
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)

###########Function to download DSS image
def download_DSS(coord, size = 800, image_server = 'DSS2 Red'):
    """
    Given a SkyCoord object, download a DSS image from the given size using
    a given server. The default is DSS2 Red to have all-sky coverage in the 
    band closest to the guider. 
    """
    img = SkyView.get_images(position=coord,survey=[image_server]\
                         ,pixels='%s,%s'%(str(size), str(size)),
                         coordinates='J2000',grid=True,gridlabels=True)
    return img

###########Main function to bring up the plot, given the coordinates
def plot_NIRES_fov(coords, coords_offset, target_name):
    """
    This function takes the coordinates of the SN and offset star, download the DSS image
    and make a plot with the slit locations and the off-axis guider FoV.
    The plot will have a slider for the user to adjust the PA. Once a PA with a guide star in
    the guider FoV at both the SN and the offset position is found, click 'Save' to save the 
    finder and return the PA. 
    Input:
        - coords:        SkyCoord of the SN
        - coords_offset: SkyCoord of the offset star
        - target_name:   Name of the target. Used to label plot and name output finder.
    """
    #Download the DSS image
    img = download_DSS(coords)
    #Get the wcs object from the image
    wcs = WCS(img[0][0].header)
    #Define the plot
    plt.figure(figsize = (10,10))
    ax = plt.subplot(projection = wcs)
    ax.imshow(img[0][0].data, origin = 'lower')

    #picture coordinates of source and offset star
    x, y = wcs.all_world2pix([coords.ra.deg], [coords.dec.deg], 1)
    if coords_offset is not None:
        xo, yo = wcs.all_world2pix([coords_offset.ra.deg], [coords_offset.dec.deg], 1)

    #dimension of slit
    dx = 0.55/2 #0.55"
    dy = 18/2 #18"

    #Define slit for source and offset star
    slit = Rectangle((x[0]-dx, y[0]-dy), 2*dx, 2*dy, color = 'w', lw = 0.5)
    if coords_offset is not None:
        slit_o = Rectangle((xo[0]-dx, yo[0]-dy), 2*dx, 2*dy, color = 'r', lw = 0.5)

    #Define guider field for source and offset star
    guider_dx = 3.44*60*1.1 #...shouldn't this number be known? This is the guider offset of around 3.8 arcmin
    guider_dim = 1.8*60 #1.8 arcmin. Yep, it's small
    #Rectangle objects representing the guider FoV for the SN and the offset star. 
    Guider = Rectangle((x[0]+guider_dx-guider_dim/2, y[0]-guider_dim/2), guider_dim,guider_dim, \
                    color = 'w', fill = False, lw = 0.5)#, \
                    #transform = transforms.Affine2D.rotate_around(x = x[0],y =  y[0],theta =  np.radians(PA)))
    if coords_offset is not None:
        Guider_o = Rectangle((xo[0]+guider_dx-guider_dim/2, yo[0]-guider_dim/2), guider_dim,guider_dim, \
                    color = 'r', fill = False, lw = 0.5)

    #Define guider circle for source and offset star
    Guider_path_in  = Circle((x[0], y[0]), guider_dx-guider_dim/2, fill = False, color = 'w', lw = 0.5,ls = '--')
    Guider_path_out = Circle((x[0], y[0]), guider_dx+guider_dim/2, fill = False, color = 'w', lw = 0.5,ls = '--')

    if coords_offset is not None:
        Guider_o_path_in  = Circle((xo[0], yo[0]), guider_dx-guider_dim/2, fill = False, color = 'r', lw = 0.5, ls = '--')
        Guider_o_path_out = Circle((xo[0], yo[0]), guider_dx+guider_dim/2, fill = False, color = 'r', lw = 0.5, ls = '--')


    #ROTATE

    # for i in range(0,370,30):
    PA = 0

    rot = transforms.Affine2D().rotate_around(x[0], y[0],np.radians(PA))+ ax.transData
    if coords_offset is not None:
       rot_o = transforms.Affine2D().rotate_around(xo[0], yo[0],np.radians(PA))+ ax.transData

    Guider.set_transform(rot)
    slit.set_transform(rot)
    if coords_offset is not None:
        Guider_o.set_transform(rot_o)
        slit_o.set_transform(rot_o)

    #Add all these overlays to the DSS image
    ax.add_patch(slit)
    ax.add_patch(Guider)
    if coords_offset is not None:
        ax.add_patch(slit_o)
        ax.add_patch(Guider_o)

    ax.add_patch(Guider_path_in )
    ax.add_patch(Guider_path_out)

    if coords_offset is not None:
        ax.add_patch(Guider_o_path_in )
        ax.add_patch(Guider_o_path_out)

    ax.tick_params(labelsize = 14)
    ax.set_xlabel('RA', fontsize = 16)
    ax.set_ylabel('Dec', fontsize = 16)
    ax.set_title('%s, PA = %d deg'%(target_name,PA), fontsize = 18)

    #Slider to adjust the PA
    plt.subplots_adjust(right=0.8, bottom=0., top = 1.1)

    axPA= plt.axes([0.85, 0.21, 0.03,0.67], facecolor='lightgoldenrodyellow')
    sPA = VertSlider(axPA, 'PA', 0, 360, valinit=0, valstep=1, valfmt = '%d', dragging = True)

    #Add Button to save and quit
    def quit_save(foo):
        plt.savefig('%s_NIRES.pdf'%target_name, bbox_inches = 'tight')
        plt.close()
        # return PA
    axsave = plt.axes([0.75, 0.1, 0.2, 0.05], facecolor = 'blue')
    bsave = Button(axsave, 'Save and Quit')
    bsave.label.set_fontsize(14)
    bsave.on_clicked(quit_save)

    def skip(foo):
        plt.close()
    # axskip = plt.axes([0.6, 0.1, 0.1, 0.05])
    # bskip = Button(axskip, 'Skip')
    # bskip.label.set_fontsize(14)
    # bskip.on_clicked(quit_skip)

    # sPA.set_val(10)
    def update(val):
        global finalPA
        PA = sPA.val
        rot = transforms.Affine2D().rotate_around(x[0], y[0],np.radians(PA))+ ax.transData
        rot_o = transforms.Affine2D().rotate_around(xo[0], yo[0],np.radians(PA))+ ax.transData
        #rerotate and draw boxes
        Guider.set_transform(rot)
        slit.set_transform(rot)
        Guider_o.set_transform(rot_o)
        slit_o.set_transform(rot_o)

        ax.add_patch(slit)
        ax.add_patch(Guider)
        ax.add_patch(slit_o)
        ax.add_patch(Guider_o)
        ax.set_title('%s, PA = %d deg'%(target_name,PA), fontsize = 16)

    #when the PA slider is changed, run update
    cid =sPA.on_changed(update)
    plt.show()
    # print(sPA.val)
    return sPA.val

if __name__ == '__main__':
    #When running this script, take a target list already with ONE offset star per SN
    #For each object, 
    #Read target and offset

    #Call plot_NIRES_fov(coords, coords_offset, target_name)

    #Update the rotdest line

    #Save the updated target list 

    print("###############NIRES PA Selection Tool###############")
    print("Your supplied target list should be in the Keck fixed-width format.")
    print("Each target could have ONE offset star in the next line.")
    print("The format is:")
    print("SN2034abc       xxxxxxx")
    print("SN2034abc_S1    xxxxxxx")

    parser = argparse.ArgumentParser(description=\
        '''
        Creates the finder charts for the whole night, provided a target list file.
        
        Usage: prepare_obs_run.py target_list_filename telescope
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    # parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help="filename of the target list")

    parser.add_argument("-d", "--debug", action="store_true",
                    help="debug mode")

    parser.add_argument("-r", "--rotate", action="store_true",
                        help="produce rotated finder chart")

    args = parser.parse_args()

    #get file name
    filename = args.filename
    # converters = {'col5': [asci.convert_numpy(np.str)]} #This is so that -00 gets preserved as -00
    # targets = asci.read(filename, comment = '!', format = 'no_header', converters = converters)
    try:
        targets = asci.read( filename, format="fixed_width_no_header",
                    col_starts=(0, 16, 28, 42, 50),
                    col_ends=(15, 27, 41, 49, 140),
                    names = ('names', 'RA', 'Dec', 'epoch', 'comments'))
        print(targets)
    except:
        print('Make sure your supplied target list is in the Keck fixed width format.')

    names = np.array(targets['names'])
    #Make SkyCoord
    coords = SkyCoord(targets['RA'], targets['Dec'], unit = (u.hourangle, u.deg))

    all_starlist = ''

    is_offset_star = False

    if args.debug:
        print(targets)

    ######Loop through the target list 
    for i in range(len(targets)):
        if is_offset_star == False: #this entry is not offset star, treat as target
            name = names[i]
            # ra_deg = coords[i].ra.deg
            # dec_deg = coords[i].dec.deg
            target_coord = coords[i]
            # mag = mags[i]
            print("Making NIRES PA Selection GUI for %s."%name)
            if i < len(targets)-1: #not the last item
                next_entry = names[i+1]
                if next_entry.split('_')[0] == name and "S" in next_entry.split('_')[1]: #next entry is offset star
                    # offset_ra = coords[i+1].ra.deg
                    # offset_dec = coords[i+1].dec.deg
                    offset_coord = coords[i+1]
                    print("Offset star provided. RA = %.5f Dec = %.5f"%(coords[i+1].ra.deg, coords[i+1].dec.deg))
                    is_offset_star = True #Because of this flag, the next item is skipped. 
                else:
                    offset_coord = None
            else: #For the last item, no host is provided. 
                offset_coord = None

            # Run the PA selection tool
            finalPA = plot_NIRES_fov(target_coord, offset_coord, name)

            starlist_entry = "{:s}{:s} {:s}  {:s}  rotmode=pa rotdest={:.2f} {:s} \n".\
                format( targets[i]['names'].ljust(16), targets[i]['RA'], targets[i]['Dec'], str(targets[i]['epoch']), finalPA,  targets[i]['comments'])
            if is_offset_star:
                starlist_entry += "{:s}{:s} {:s}  {:s}  rotmode=pa rotdest={:.2f} {:s} \n".\
                format( targets[i+1]['names'].ljust(16), targets[i+1]['RA'], targets[i+1]['Dec'], str(targets[i+1]['epoch']), finalPA,  targets[i+1]['comments'])
            if args.debug:
                print(starlist_entry)
            all_starlist += starlist_entry

        elif is_offset_star == True: #In this case, the next item should be a target
            is_offset_star = False 
        ###Keep writing to file so we don't lose progress.
        out_file = open(filename.split('.')[0]+'_final.txt', "w")
        out_file.write(all_starlist)        
        out_file.close()

        if args.rotate:
                print("Making rotated finder chart for NIRES")
                try:
                    finder_name = name+'_finder.png'
                    rotated_finder = name+'_finder_rot.png'
                    # NIRES SVC is north down, so first we flip finders 180
                    # then we rotate by rotdest clockwise (which is negative in PIL rotate)
                    finder_rot=180-finalPA
                    im = Image.open(finder_name)
                    im_rotate=im.rotate(finder_rot, resample=Image.BICUBIC, expand = True,fillcolor=(255,255,255))
                    im_rotate.save(rotated_finder,dpi=(400,400))

                except:
                    print("Check if PIL library is installed.")

    print(all_starlist)


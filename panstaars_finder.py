###############Code by Georgios Dimitriadis

def panstamps_lite(self, name='finder_field.fits'):
        import requests
        import re
        ra=self.ra
        dec=self.dec
        filt='r'
        pos = """%(ra)s %(dec)s""" % locals()
        size=self.size
        try:
            response = requests.get(
                url="http://plpsipp1v.stsci.edu/cgi-bin/ps1cutouts",
                params={
                    "pos": pos,
                    "filter": filt,
                    "filetypes": "stack",
                    "size": int(size)*60*4,
                    "output_size": int(size)*60*4,
                    "verbose": "0",
                    "autoscale": "99.500000",
                    "catlist": "",
                },
            )
        except requests.exceptions.RequestException:
            print('HTTP Request failed')
        reFitscutouts = re.compile(
            r"""<th>(?P<imagetype>\w+)\s+(?P<skycellid>\d+.\d+)\s+(?P<ffilter>[\w\\]+)(\s+(?P<mjd>\d+\.\d+))?<br.*?href="(http:)?//plpsipp1v.*?Display</a>.*?Fits cutout" href="(?P<fiturl>(http:)?//plpsipp1v.*?\.fits)".*?</th>""", re.I)
        #if sys.version_info[0] < 3:
        #    thisIter = reFitscutouts.finditer(response.content)
        #else:
        #    thisIter = reFitscutouts.finditer(response.content.decode('utf-8'))
        thisIter = reFitscutouts.finditer(response.content)
        stackFitsUrls = []
        for item in thisIter:
            imagetype = item.group("imagetype")
            skycellid = item.group("skycellid")
            ffilter = item.group("ffilter")
            fiturl = 'http://plpsipp1v.stsci.edu%s'%item.group("fiturl")
            if fiturl[0:5] != "http:":
                fiturl = "http:" + fiturl
                mjd = item.group("mjd")
            stackFitsUrls.append(fiturl)
        if not len(stackFitsUrls):
            return(None)
        try:
            os.remove( name )
        except OSError:
            pass
        print 'Downloading image.'
        for s in stackFitsUrls:
            s = s.replace('plpsipp1v.stsci.edu//','')
            #if not os.path.dirname(name):
            #    name = '%.7f_%.7f_%s.PS1.fits'%(ra,dec,time.time())
            #else:
            #    name = '%s/%.7f_%.7f_%s.PS1.fits'%(os.path.dirname(name),ra,dec,time.time())
            wget.download(s,out=name)
            self.image = name
            break
        if os.path.exists(name):
            self.image = name
        else: return(None)


def get_decam_image(self, name='finder_field.fits'):
        """
        Get a decam finder chart, if not given input image. Size is x*0.27*0.0166 arcmin.
        """
        size=int((self.size)/0.0045)
        url = "http://legacysurvey.org/viewer/fits-cutout?ra=%.8f&dec=%.8f"%(self.ra, self.dec) +\
              "&size="+str(size)+"&layer=dr8&pixscale=0.3&bands=r" #%(size)
        #url = "http://legacysurvey.org/viewer/fits-cutout?ra=%.8f&dec=%.8f"%(self.ra, self.dec) +\
        #      "&size=500""&layer=decals-dr7&pixscale=0.27&bands=g" #%(size)
        #http://legacysurvey.org/viewer/jpeg-cutout?ra=343.9877625&dec=-43.4386083&size=100&layer=dr8&pixscale=0.3&bands=grz
        print 'Downloading image.'
        try:
            os.remove( name )
        except OSError:
            pass
        wget.download( url, out=name )
        self.image = name
__version__ = "0.1.0 (2017/07/24)"
__authors__ = ['Jakub Prchlik <jakub.prchlik@cfa.harvard.edu>']
__email__   = "jakub.prchlik@cfa.harvard.edu"

import matplotlib
#fixes multiprocess issue
matplotlib.use('agg')

import sys

from make_movie import create_movie


import subprocess
import glob
import os
import stat
import numpy as np
from datetime import date,datetime
from datetime import timedelta as dt
from multiprocessing import Pool
import grab_goes_xray_flux as ggxf
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from astropy.io import ascii
from astropy.table import vstack,Table,join

from SMEARpy import Scream

#import aia_mkimage
from aia_mkimage import aia_mkimage


#Class for making AIA movies
class aia_mkmovie:

    #initialize aia_mkmovie
    def __init__(self,start,end,wav,cadence='6m',w0=1900,h0=1144,dpi=300,usehv = False,
                 panel=False,color3=False,select=False,videowall=True,nproc=2,goes=False,wind=False,
                 archive="/data/SDO/AIA/synoptic/",dfmt = '%Y/%m/%d %H:%M:%S',
                 outf=True,synoptic=False,odir='working/',frate=10,time_stamp=True,cx=0.0,cy=0.0,
                 prompt=False,cutout=False,rotation=False,rot_time=None,download=False,
                 local=False,email=None):
        """ 
        Return an object with parameter to download data, create images, and create a movie

        Parameters
        ----------
        start : single string or datetime object
            The start time over which to create the movie.
            The argument may be a string or datetime object.
            If start is a string then the string must be in
            the same form as the dfmt parameter (dfmt parameter 
            default = %Y/%m/%d %H:%M:%S)
        end   : single string or datetime object
            The end time over which to create the movie.
            The argument may be a string or datetime object.
            If end is a string then the string must be in
            the same form as the dfmt parameter (dfmt parameter 
            default = %Y/%m/%d %H:%M:%S)
        wav   : float, string, or array-like
            The wavelengths to use in the movie.
            Use a single float or string when making a movie of
            a single wavelength (e.g. 171).
            Use an array with 3 elements if making a RGB movie.
            The 3 wavelength array must be sort R,G,B (e.g. 
            [304,193,171]). Also, set optional parameter color3=True. 
        cadence: string, int or float, optional 
            The sampling frequency for downloading and image creation.
            If cadence is an int or float the program assumes the 
            input cadence is in seconds. However, if one passes a 
            string then the cadence can be set in seconds (s),
            minutes (m), hours (h), or days (d). For example you may 
            set with cadence = '6h'. Default = '6m'
        w0: int or float, optional
            Width of the movie in pixels. If the height (h0) is larger than
            w0 the program will switch the two parameters on output. 
            However, it will also transpose the x and y axes, which allows 
            for rotated images and movies. Default = 1900
        h0: int or float, optional 
            Height of the movie in pixels. If h0 is larger than the
            width (w0) the program witll switch the two parameters on
            output. However, it will also transpose the x and y axies,
            which allows for rotated images and movies. Default = 1144
        dpi: int or float, optional
            Dots per inch in the output images. Default = 300.
        usehv: boolean, optional 
            Use helioviewer to download images (Currently not implemented)
        panel: boolean, optional
            Make a 4 panel plot. If panel set to True then color3 must be 
            False and the wavelength list must be 4 wavelengths long.
            The wav list has the following format [top right, top left,
            bottom right, bottom left]. Default = False
        color3 : boolean, optional
            Create a 3 color image. If color3 set to True panel must be
            False and the wavelength list must be 4 wavelengths long.
            The wav list has the following format [R, G, B]. Default =
            False.
        select: boolean, optional
            Select region of sun to focus on for the movie. If True
            the program will launch a GUI to manually select a region
            and set the color table
        videowall: boolean, optional
            Use the videowall height and width. Overrides setting h0 and w0.
        nproc:  int, optional
            Number of processors to use while downloading, creating images,
            and creating movies. If nproc = 1 then the program will loop 
            instead of pooling image creation. Due to problems with the 
            matplotlib backends parallel processing occasionally scrambles
            text in images. Default = 2.
        goes: boolean, optional
            Overplot the goes flux with a 1 minute cadence.
            goes only works for a full sun single wavelength image.
            Default = False
        wind: boolean, optional
            Overplot the 6 minute solar wind parameters from ACE. wind only
            works if goes is True and full sun for a single wavelength.
        archive: string, optional
            archive is the location of the aia fits files. If download
            is set then the files download to a default directory. Then
            archive does not need to be set. Default = "/data/SDO/AIA/synoptic/"
        dfmt : string, optional 
            The string format of start and end. Can be set to any date time
            stripping variable in python, link below: 
            (https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior).
            Default = '%Y/%m/%d %H:%M:%S'
        outf: boolean, optional
            Output movie filename without .mp4. (e.g. cool_flare). The movie outputs
            to the YYYY_MM_DD_HHMM_*/final/ directory. Default YYYMMDD_HHMM.
        synoptic: boolean, optional
            Check using synoptic parameters or not (synoptic are 1024x1024 images).
            Default = False.
        odir:  string, optional, deprecated
            output directory, which we now force to YYYY_MM_DD_HHMM_*
        frate: int or float,optional
            frame rate of output movie. Default = 10 frames per sec.
        time_stamp: boolean, optional 
            Include time stamp in images. Default = True
        
        cx  : float or int, optional
            Center of the field of view for creating images. If cx is set
            then the image is assumed to be a cutout. Selecting in prompt
            overrides cx. Default = 0.0.
        cy  : float or int, optional
            Center of the field of view for creating images. If cy is set
            then the image is assumed to be a cutout. Selecting in prompt
            overrides cy. Default = 0.0.
        prompt : boolean, optional
            Bring up a prompt to set w0,h0,cx,cy, and the color scale.
            Default = False
        cutout : boolean, optional
            Select a cutout from the full sun file. Automatically set if
            rotation or prompt is True. Default = False
        rotation : boolean, optional
            Use Sunpy rotation for rotation correction. Automatically True
            if rot_time is set. Default = False.
        rot_time : datetime object or string, optional
            The time cx and cy are measured. Can be set in prompt or manually.
            If manually set then the rot_time must be a datetime object or 
            a string with format dfmt. Default = None.
        download: boolean, optional
            Download the aia date from jsoc. Need to supply email if True.
            Default = False
        local   : boolean, optional
            The archive is a local directory. Only need if you previously 
            downloaded files. Then use archive to set the directory.
            Default = False
        email   : string, optional
            Email address register to JSOC, so you may download aia files.
            Default = None
        """


        #set to default scaling in aia_mkimage
        self.img_scale = None
     
        #No default x or y limits
        self.xlim = None
        self.ylim = None

        #list of acceptable wavelengths
        self.awavs = ['94','131','171','193','211','304','335','1600','1700']

        #check output directory
        if isinstance(odir,str):
            self.odir = odir
        else:
            sys.stdout.write('odir must be a string')
            sys.exit(1)
 
        #format and create output directory
        if self.odir[-1] != '/': self.odir=self.odir+'/'
        if not os.path.isdir(self.odir): os.mkdir(self.odir)

        #use synoptic image checking (default = True)
        self.synoptic = synoptic

        #make sure datetime formatter is string
        if isinstance(dfmt,str):
            self.dfmt = dfmt
        else:
            sys.stdout.write('datetime formatter must be string')
            sys.exit(1)

        #check inserted start time
        if isinstance(start,datetime):
            self.start = start
        elif isinstance(start,str):
            self.start = datetime.strptime(start,dfmt)
        else:
            sys.stdout.write('Start time must be datetime object or formatted string')


        #check output file
        if isinstance(outf,str):
            self.outf = outf
        elif outf: 
            self.outf = self.start.strftime('%Y%m%d_%H%M')
        else:
            sys.stdout('outf must be a string or undefined')
            sys.exit(1)
 

        #check inserted end time
        if isinstance(end,datetime):
            self.end = end
        elif isinstance(end,str):
            self.end = datetime.strptime(end,dfmt)
        else:
            sys.stdout.write('End time must be datetime object or formatted string')



        #check if cadence is a string
        #if not so convert the cadence 
        #assuming it is given in seconds
        if isinstance(cadence,str):
            self.cadence = cadence
        elif isinstance(cadence,(int,float)):
            self.cadence = str(cadence)+'s'
        else:
            sys.stdout.write('Cadence must be a string, integer, or float')
            sys.exit(1)

        #check image height
        if isinstance(h0,(int,float)):
            self.h0 = h0
        else:
            sys.stdout.write('h0 must be an integer or float')
            sys.exit(1)
         
        #check image height
        if isinstance(frate,(int,float)):
            self.frate = int(frate)
        else:
            sys.stdout.write('frate must be an integer or float')
            sys.exit(1)
         
        #check image width
        if isinstance(w0,(int,float)):
            self.w0 = w0
        else: 
            sys.stdout.write('w0 must be an integer or float')
            sys.exit(1)

 
        #check dpi
        if isinstance(dpi,(int,float)):
            self.dpi = dpi
        else: 
            sys.stdout.write('dpi must be an integer or float')
            sys.exit(1)

        #download files from helioviewer (probably switch to something else)
        self.usehv = usehv
        #create a panel movie
        self.panel = panel
        #create 3 color image (default = False)
        self.color3 = color3
        self.select = select
        #number of processors to use when creating images
        if isinstance(nproc, int):
            self.nproc = nproc
        else:
            sys.stdout.write('nproc must be an integer')
            sys.exit(1)
         
        #create a goes plot
        self.goes = goes
        if not goes: self.goesdat =  [] 
        #create solar wind data plot
        self.wind = wind
        if not wind: self.aceadat = []     
        #automatically turn on goes plot if wind plot is set
        if self.wind: self.goes = True

        #check if timestamp flag is set (Default = True)
        if isinstance(time_stamp,bool): 
            self.time_stamp = time_stamp
        else:
            sys.stdout.write('timestamp must be a boolean')
            sys.exit(1)

        #check if rotation flag is set (Default = False)
        if isinstance(rotation,bool): 
            self.rotation = rotation
        else:
            sys.stdout.write('rotation must be a boolean')
            sys.exit(1)

        #check if download flag is set (Default = False)
        if isinstance(download,bool): 
            self.download = download
        else:
            sys.stdout.write('download must be a boolean')
            sys.exit(1)

        #check if local flag is set (Default = False)
        if isinstance(local,bool): 
            self.local = local
        else:
            sys.stdout.write('local must be a boolean')
            sys.exit(1)


       
        #check inserted rot_time time
        if isinstance(rot_time,datetime):
            self.rot_time = rot_time
            self.rotation = True
        elif isinstance(rot_time,str):
            self.rot_time = datetime.strptime(rot_time,dfmt)
            self.rotation = True
        elif rot_time is None:
            self.rot_time = rot_time
        else:
            sys.stdout.write('rot_time must be datetime object or formatted string')
            sys.exit(1)


        #check inserted email
        if isinstance(email,str):
            self.email = email
        elif email is None:
            self.email = email
            if download: 
                sys.stdout.write('email must be supplied if download is True')
                sys.exit(1)
        else:
            sys.stdout.write('email must be a string')
            sys.exit(1)






        #check if prompt flag is set (Default = False)
        if isinstance(prompt,bool): 
            self.prompt = prompt
        else:
            sys.stdout.write('prompt must be a boolean')
            sys.exit(1)

        #check if cutout flag is set (Default = False)
        if isinstance(cutout,bool): 
            self.cutout = cutout
        else:
            sys.stdout.write('cutout must be a boolean')
            sys.exit(1)


        #Do not let panel and goes/wind plots work together
        if ((self.panel) & (self.goes)):
            sys.stdout.write('Panel and goes plot cannot be used together. Choose wisely...')
            sys.exit(1)



        #check image x center
        if isinstance(cx,(int,float)):
            self.cx = cx
        else:
            sys.stdout.write('cx must be an integer or float (Assuming 0)')
            self.cx = 0.0

        #check image y center
        if isinstance(cy,(int,float)):
            self.cy = cy
        else:
            sys.stdout.write('cy must be an integer or float (Assuming 0)')
            self.cy = 0.0




        #location of SDO files
        if isinstance(archive,str):
            #make sure archive string ends in /
            if archive[-1] != '/': archive=archive+'/'
            self.archive = archive 
        else:
            sys.stdout.write('archive must be a string')
            sys.exit(1)
      

        #dimensions for videowall (over rides setting w0 and h0
        if videowall:
            self.w0 = 1900
            self.h0 = 1144

        #use a single wavelength
        self.single = False

        #check input wavelength formatting
        #check formatting assuming float or int
        if isinstance(wav,(int,float)):
            self.single = True
            self.wav  = str(int(wav))
            if ((self.panel) | (self.color3)):
                sys.stdout.write('Single wavelength cannot be used with panel or 3 color outputs')
                sys.exit(1)
            #check to make sure wavelength is allowed
            if self.wav not in self.awavs:
                sys.stdout.write('{0} not an acceptable wavelength'.format(self.wav))
                sys.exit(1)
        #check formatting assuming string      
        elif isinstance(wav,str):
            self.single = True
            self.wav = wav
            if ((self.panel) | (self.color3)):
                sys.stdout.write('Single wavelength cannot be used with panel or 3 color outputs')
                sys.exit(1)
            #check to make sure wavelength is allowed
            if self.wav not in self.awavs:
                sys.stdout.write('{0} not an acceptable wavelength'.format(self.wav))
                sys.exit(1)

        #check formatting assuming array
        elif isinstance(wav,(list,np.ndarray)):
            self.wav = []
            for i in wav:
                if isinstance(i,(float,int)):
                    self.wav.append(str(int(i)))
                elif isinstance(i,str):
                    self.wav.append(i)
                else:
                    sys.stdout.write('Wavelengths must be string, floats, or integers')
                    sys.exit(1) 
                #check to make sure wavelength is allowed
                if self.wav[-1] not in self.awavs:
                    sys.stdout.write('{0} not an acceptable wavelength'.format(i))
                    sys.exit(1)

            if not ((self.panel) | (self.color3)):
                sys.stdout.write('Panel or 3 color not set for multiple wavelengths. Please specify one')
                sys.exit(1)
   

        #directory for file output
        wavapp = ''
        if self.single: wavapp = '_{0}'.format(self.wav)
        #loop over all wavelengths if there is more than 1 value in the list
        else: 
            for i,j in enumerate(self.wav):  wavapp = wavapp+'_{0}'.format(j)
        self.sdir = self.start.strftime('%Y%m%d_%H%M')+wavapp



        # if rotation is set and prompt is not set rot_time must be set
        if ((rot_time is None) & (not self.prompt) &(self.rotation )):
            sys.stdout.write('Rotation time must be set if rotation is set and prompt is not')
            sys.exit(1)

        #exit if both color3 and panel set
        if ((self.color3) & (self.panel)):
            sys.stdout.write('Both color3 and panel booleans cannot be set to True')
            sys.exit(1)
            
       
#create directories without erroring if they already exist c
    def create_dir(self,dirs):
        try:
            os.mkdir(dirs)
        except OSError:
            sys.stdout.write('{0} Already Exists'.format(dirs))

    

    def gather_files(self):
        # use helioviewer if requested 
        if self.usehv:
            from sunpy.net.helioviewer import HelioviewerClient
            import matplotlib.image as mpimg
            hv = HelioviewerClient()
            dayarray = glob.glob(self.sdir+'/raw/*jp2')
            forpool = np.arange(len(dayarray))
            #for i in forpool: format_img(i)
            pool2 = Pool(processes=nproc)
            outs = pool2.map(get_file,forpool)
            pool2.close()



        #video wall ratio
        self.rv = float(self.w0)/float(self.h0)


        #scale up the images with increasing dpi
        self.sc = self.dpi/100

        self.span = (self.end-self.start).total_seconds() #seconds to run the movie over




        #create a directory which will contain the raw png files
        #sdir = stard+eday.date().strftime('%Y%m%d')
        #creating a subdirectory to extra step is not need
        dirlist = [self.sdir,self.sdir+'/raw',self.sdir+'/working',self.sdir+'/working/symlinks',self.sdir+'/final',self.sdir+'/goes',self.sdir+'/ace']
        for i in dirlist: self.create_dir(i)

        

        #get all days in date time span
        if self.goes: 
            ggxf.look_xrays(self.start,self.end+dt(days=1),self.sdir)
            goesfil = glob.glob(self.sdir+'/goes/*txt')
            goesnames = [ 'YR', 'MO', 'DA', 'HHMM', 'JDay', 'Secs', 'Short', 'Long'] 
            self.goesdat = Table(names=goesnames)
        
        #loop over all day information and add to large array
            for m in goesfil:
                temp = ascii.read(m,guess=True,comment='#',data_start=2,names=goesnames)
                self.goesdat = vstack([self.goesdat,temp])
        
            #create datetime array
            self.goesdat['time_dt'] = [datetime(int(i['YR']),int(i['MO']),int(i['DA']))+dt(seconds=i['Secs']) for i in self.goesdat]
        
        if self.wind:
            aceb = glob.glob(sdir+'/ace/*mag*txt')
            acep = glob.glob(sdir+'/ace/*swe*txt')
        
            aceb_names = [ 'YR', 'MO', 'DA', 'HHMM', 'JDay', 'Secs', 'S', 'Bx','By','Bz','Bt','Lat','Long'] 
            acep_names = [ 'YR', 'MO', 'DA', 'HHMM', 'JDay', 'Secs', 'S', 'Den','Speed','Temp'] 
          
            acebdat = Table(names=aceb_names)
            acepdat = Table(names=acep_names)
        
        #put B field in large array
            for m in aceb:
                temp = ascii.read(m,guess=True,comment='#',data_start=2,names=aceb_names)
                acebdat = vstack([acebdat,temp])
        #put plasmag in large array
            for m in acep:
                temp = ascii.read(m,guess=True,comment='#',data_start=2,names=acep_names)
                acepdat = vstack([acepdat,temp])
        
        
            self.aceadat = join(acepdat,acebdat,keys=['YR','MO','DA','HHMM'])
            #create datetime array
            self.aceadat['time_dt'] = [datetime(int(i['YR']),int(i['MO']),int(i['DA']))+dt(seconds=i['Secs_1']) for i in self.aceadat]
        
        #download files locally
        if self.download:
            self.run_download()
        #use a local directory of files
        elif self.local:
            self.gather_local()
        else:
            self.scream() # use SDO archive to search




    #run file download
    def run_download(self):
        """
        Import module/class to download fits files from JSOC
        """
        import aia_download_files as adf
        #create an archive in the download directory
        self.archive = self.sdir+'/raw/'
        fobj = adf.download_files(self.start,self.end,self.wav,self.cadence,email=self.email,odir=self.archive,max_con=self.nproc)
        fobj.get_drms_files()
  
        #check the downloaded files
        self.gather_local()


    #make sure the input wavelength matches the searched wavelength
    def check_wavelength(self,fil,wav):
        """
        Check wavelength and cadence of images
        """
        from astropy.io import fits

        new_fil = []
        fil_dat = []
        datefmt = '%Y-%m-%dT%H%M%SZ'
        #retrieve file wavelength and observation time
        for i in fil:
            #try assuming download format of jsoc files
            try:
                date = datetime.strptime(i.strip(self.archive).split('.')[2],datefmt)
                wave = i.strip(self.archive).split('.')[3]
                if int(wave) == int(wav):
                    new_fil.append(i)
                    fil_dat.append(date)
            #else open the fits file and read information from the header
            except:
                data = fits.open(i)
                date = datetime.strptime(data[1].header['T_OBS'].split('.')[0], '%Y-%m-%dT%H:%M:%S')
                wave = data[1].header['wavelnth']
                if int(wave) == int(wav):
                    new_fil.append(i)
                    fil_dat.append(date)

        #convert fil_dat to numpy time array 
        timelist = np.array(fil_dat)
 
        final_list = [] #list of index to keep for image creation
        for p in self.real_cad: #loop over all cadence values to find best array values
            k = np.abs(timelist-p)
            rindex, = np.where(k == k.min()) #get the nonzero array index value
            final_list.append(new_fil[rindex[0]])

        return final_list


#retrieve desired cadence from file list
    def des_cad(self,start,end,delta):
        """Create an array from start to end with desired cadence"""
        curr = start
        while curr < end:
            yield curr
            curr += delta


    #Create a cadence from time list
    def create_cadence(self):
    #turn cadence into seconds
        if self.cadence[-1] == 's': self.cad = float(self.cadence[:-1])
        elif self.cadence[-1] == 'm': self.cad = 60.*float(self.cadence[:-1])
        elif self.cadence[-1] == 'h': self.cad = 3600.*float(self.cadence[:-1])
        elif self.cadence[-1] == 'd': self.cad = 24.*3600.*float(self.cadence[:-1])
      

    #desired cadence for the observations
        self.real_cad = [result for result in self.des_cad(self.start,self.end,dt(seconds=self.cad))]



    def gather_local(self):
        from glob import glob

        #create a list of desired cadences
        self.create_cadence()

        #Singular search if only one wavelength specified
        if self.single:
            fits_files = glob(self.archive+'*'+str(int(self.wav))+'*')
            if len(fits_files) < 1.:
                sys.stdout.write('No AIA Files Found')
                sys.exit(1)
            else:
                #make sure the wavelength header agrees with found value
                fits_files = self.check_wavelength(fits_files,self.wav)
                self.fits_files = fits_files

        #if multiple wavelengths loop over all wavelengths
        elif ((self.color3) | (self.panel)):
            self.fits_files_temp = []
           #loop over all wavelengths in array
            for i in self.wav:
                fits_files = glob(self.archive+'*'+str(int(i))+'*')
                if len(fits_files) < 1.:
                    sys.stdout.write('No AIA Files Found')
                    sys.exit(1)
                else:
                #make sure the wavelength header agrees with found value
                    fits_files = self.check_wavelength(fits_files,i)
                    self.fits_files_temp.append(fits_files)
            #transpose list array
            self.fits_files = map(list,zip(*self.fits_files_temp))


    def scream(self):

        #J. Prchlik 2016/10/06
        #Updated version calls local files
        verbose=False 
        debug = False 
        src = Scream(archive=self.archive,verbose=verbose,debug=debug)
        ##########################################################
        # Phase 1: get file names                                #
        ##########################################################
        #if self.span > 7200.:
        #    sendspan = "-{0:1.0f}m".format(self.span/60.+1.) # use minutes if greater than 7200 seconds
        #else:
        sendspan = "-{0:1.0f}s".format(self.span) # need to spend current span not total span
        paths = src.get_paths(date=self.end.strftime("%Y-%m-%d"), time=self.end.strftime("%H:%M:%S"),span=sendspan)


        #Singular search if only one wavelength specified
        if self.single:
            fits_files = src.get_filelist(date=self.end.strftime("%Y-%m-%d"),time=self.end.strftime("%H:%M:%S"),span=sendspan,wavelnth=self.wav)
            qfls, qtms = src.run_quality_check(synoptic=self.synoptic)
            if len(qfls) < 1.:
                sys.stdout.write('No AIA Files Found')
                sys.exit(1)
            else:
                self.fits_files = src.get_sample(files = qfls, sample = self.cadence, nfiles = 1)
        elif ((self.color3) | (self.panel)):
            self.fits_files_temp = []
           #loop over all wavelengths in array
            for i in self.wav:
                fits_files = src.get_filelist(date=self.end.strftime("%Y-%m-%d"),time=self.end.strftime("%H:%M:%S"),span=sendspan,wavelnth=i)
                qfls, qtms = src.run_quality_check(synoptic=self.synoptic)
                if len(qfls) < 1.:
                    sys.stdout.write('No AIA Files Found')
                    sys.exit(1)
                else:
                    self.fits_files_temp.append(src.get_sample(files = qfls, sample = self.cadence, nfiles = 1))
            #transpose list array
            self.fits_files = map(list,zip(*self.fits_files_temp))
             


    #prompt for selecting area
    def init_prompt(self):
        #check the python version to use one Tkinter syntax or another
        if sys.version_info[0] < 3:
            import Tkinter as Tk
        else:
            import tkinter as Tk
        import aia_select_cutout as asc


        #create list of files based roughly on time (indices 0-3 are different wavelengths at roughly the same time)
        if (self.panel): templist = reduce(lambda x,y: x+y, self.fits_files)
        else: templist = self.fits_files

        #init gui instance
        gui = asc.gui_c(Tk.Tk(),templist,color3=self.color3,w0=self.w0,h0=self.h0,cx=self.cx,cy=self.cy,img_scale=self.img_scale)
        #run the gui to select region
        gui.mainloop()

        #set parameters based on gui output
        self.cx = gui.cx
        self.cy = gui.cy
        self.w0 = gui.w0
        self.h0 = gui.h0
        self.xlim = [min(gui.xbox),max(gui.xbox)]
        self.ylim = [min(gui.ybox),max(gui.ybox)]

        #use image scaling to pass to mkimage
        self.img_scale=gui.img_scale
        
        #set cutout to be true if prompt selected 
        self.cutout = True

        #if rotation is set return the rotation 0 time
        if self.rotation: self.rot_time = gui.obs_time


    def create_images_movie(self):


        #create a list of class objects
        #rotate x, y if rotation is set
        #else do the normal thing
        image_list = [aia_mkimage(i,w0=self.w0,h0=self.h0,dpi=self.dpi,
                     sc=self.sc,goes=self.goes,goesdat=self.goesdat,sday=self.start,eday=self.end,
                     img_scale=self.img_scale,cutout=self.cutout,
                     ace=self.wind,aceadat=self.aceadat,single=self.single,panel=self.panel,
                     color3=self.color3,time_stamp=self.time_stamp,odir=self.sdir+'/working/',
                     cx=self.cx,cy=self.cy,
                     xlim=self.xlim,ylim=self.ylim,synoptic=self.synoptic,rot_time=self.rot_time) for i in self.fits_files]

        #J. Prchlik 2016/10/06
        #Switched jp2 to fits
        #loop is for testing purposes
        #for i in forpool: format_img(i)
        if self.nproc > 1:
            pool1 = Pool(processes=self.nproc)
            outs = pool1.map(format_img,image_list)
            pool1.close()
        #just loop is 1 processor specified
        else:
            for i in image_list: format_img(i)

        #flip w0 and h0 if h0 > w0
        if self.h0 > self.w0:
            th0 = self.h0
            tw0 = self.w0
            self.w0 = th0
            self.h0 = tw0

        #create movie object
        mo =create_movie(odir = self.sdir+'/final/',pdir = self.sdir+'/working/', ext = 'png', w0 = int(self.w0), h0=int(self.h0),frate=self.frate,outmov=self.outf)
        #run movie object
        mo.create_movie()

    def run_all(self):
        """
        Runs program from beginning to end based on input parameters
        """
        
        #get fits files
        self.gather_files()


        #if prompt set bring up a prompt
        if self.prompt: self.init_prompt()

        #run image creation
        self.create_images_movie()
 


def format_img(aia_img):
    aia_img.format_img()



    


        
        


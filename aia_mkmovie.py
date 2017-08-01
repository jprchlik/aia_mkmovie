__version__ = "0.0.1 (2017/07/24)"
__authors__ = ['Jakub Prchlik <jakub.prchlik@cfa.harvard.edu>']
__email__   = "jakub.prchlik@cfa.harvard.edu"

import matplotlib
#fixes multiprocess issue
matplotlib.use('Tkagg')

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
    def __init__(self,start,end,wav,cadence='6m',w0=1900,h0=1144,dpi=300,usehv = False,panel=False,color3=False,select=False,videowall=True,nproc=2,goes=False,wind=False,x0=0.0,y0=0.0,archive="/data/SDO/AIA/synoptic/",dfmt = '%Y/%m/%d %H:%M:%S',outf=True,synoptic=True,odir='working/',frate=10,time_stamp=True,cx=0.0,cy=0.0,prompt=False,cutout=False,rotation=False,rot_time=None):
        """ 
        Take 3 color input of R,G,B for wavelength

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
            self.outf = self.start.date().strftime('%Y%m%d_%H%M')
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


       
        #check inserted rot_time time
        if isinstance(rot_time,datetime):
            self.rot_time = rot_time
        elif isinstance(rot_time,str):
            self.rot_time = datetime.strptime(rot_time,dfmt)
        elif rot_time is None:
            self.rot_time = rot_time
        else:
            sys.stdout.write('rot_time must be datetime object or formatted string')
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
        self.sdir = self.start.date().strftime('%Y%m%d_%H%M')+wavapp



        # if rotation is set and prompt is not set rot_time must be set
        if ((rot_time is None) & (not self.prompt) &(self.rotation )):
            sys.stdout.write('Rotation time must be set if rotation is set and prompt is not')
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

#datasources = hv.get_data_sources()


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


        #Singular search if only one wavelength specificied
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
        if self.panel: templist = reduce(lambda x,y: x+y, self.fits_files)
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


        #create movie object
        mo =create_movie(odir = self.sdir+'/final/',pdir = self.sdir+'/working/', ext = 'png', w0 = self.w0, h0=self.h0,frate=self.frate,outmov=self.outf)
        #run movie object
        mo.create_movie()

    def run_all(self):
        
        self.gather_files()
        if self.prompt: self.init_prompt()
        self.create_images_movie()
 


def format_img(aia_img):
    aia_img.format_img()



    


        
        


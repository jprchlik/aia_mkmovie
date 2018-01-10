__version__ = "0.1.0 (2017/09/06)"
__authors__ = ['Jakub Prchlik <jakub.prchlik@cfa.harvard.edu>']
__email__   = "jakub.prchlik@cfa.harvard.edu"


#import dkrms
import astropy.units as u
import sys
from datetime import datetime
import numpy as np
import os




class download_files:
    def __init__(self,start,end,wav,cadence,series='aia.lev1_euv_12s',
                 segment='image',email=None,odir=None,
                 overwrite=True,max_con=1,dfmt='%Y/%m/%d %H:%M:%S'):
        """
        Both a stand alone and aia_mkmovie connected class for downloading AIA data from the JSOC archive.
        N.B. JSOC may fail without proper warning if total download size is too larger even though
        the number of files is within the required parameters. If that occurs break the download into
        smaller wavelength groupings for download. Then you may combined the data into one folder
        and use the local archive keywords in aia_mkmovie.py

        Parameters
        ----------
        start : single string or datetime object
            The start time over which to download.
            The argument may be a string or datetime object.
            If start is a string then the string must be in
            the same form as the dfmt parameter (dfmt parameter 
            default = %Y/%m/%d %H:%M:%S)
        end   : single string or datetime object
            The end time over which to download.
            The argument may be a string or datetime object.
            If end is a string then the string must be in
            the same form as the dfmt parameter (dfmt parameter 
            default = %Y/%m/%d %H:%M:%S)
        wav   : float, string, or array-like
            AIA wavelengths to download. Can be a single wavelength
            or a list of wavelengths.
        cadence: string, int or float
            The sampling frequency for downloading images.
            If cadence is an int or float the program assumes the 
            input cadence is in seconds. However, if one passes a 
            string then the cadence can be set in seconds (s),
            minutes (m), hours (h), or days (d). For example you may 
            set with cadence = '6h'. Default = '6m'
        series: string, optional
            Series to download data from. The value currently must be 
            'aia.lev1_uv_24s' for 1600 and 1700 observations or
            'aia_lev1_euv_12s' for 94,131,171,193,211,304, and 335
            observations. aia.lev1 contains both, but contains more 
            ancillary files. Default = 'aia.lev1_euv_12s' 
        segment: string, optional 
            Type of data to download from time range (e.g. spike or image).
            Default = image
        email  : string, semi-optional
            The email to send the downloaded data to. This field is required
            for JSOC download, which currently is the only option. However,
            future versions may include jhelioviewer and/or VSO downloads.  
        odir   : string, optional
            Output directory for downloaded files. Default = './'
        overwrite: boolean, optional
            Overwrite files of the same name in the local directory.
            Default = True
        max_con : int, optional
            Maximum number of connections to JSOC archive. Default = 1.
        dfmt : string, optional 
            The string format of start and end. Can be set to any date time
            stripping variable in python, link below: 
            (https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior).
        """

        #list of acceptable wavelengths
        self.awavs = [94,131,171,193,211,304,335,1600,1700]
        #list of acceptable segments
        self.asegs = ['image','spike','None']
        #list of acceptable series
        self.asers = ['aia.lev1_uv_24s','aia.lev1_euv_12s','aia.lev1']


        #check if overwrite flag is set (Default = True)
        if isinstance(overwrite,bool):
            self.overwrite = overwrite
        else:
            sys.stdout.write('overwrite must be a boolean')
            quit()

        #check if max_conn is int
        if isinstance(max_con,int):
            self.max_con = max_con
        else:
            sys.stdout.write('max_con must be a boolean')
            quit()

        #check segment
        if isinstance(segment,str):
            self.segment = segment
            if self.segment not in self.asegs:
                sys.stdout.write('segment not in acceptable segment list')
                quit()
        else:
            sys.stdout.write('segment must be a string')
            quit()

        
        #check series
        if isinstance(series,str):
            self.series = series
            if self.series not in self.asers:
                sys.stdout.write('series not in acceptable series list')
                quit()
        else:
            sys.stdout.write('series must be a string')
            quit()


        #check output directory
        if isinstance(odir,str):
            self.odir = odir
        elif odir is None:
            self.odir = './'
        else:
            sys.stdout.write('odir must be a string')
            quit()

        #check email is string
        if isinstance(email,str):
            self.email = email
        else:
            sys.stdout.write('email must be a string')
            quit()



        #check that start and end are datetime objects
       #make sure datetime formatter is string
        if isinstance(dfmt,str):
            self.dfmt = dfmt
        else:
            sys.stdout.write('datetime formatter must be string')
            quit()

        #check inserted start time
        if isinstance(start,datetime):
            self.start = start
        elif isinstance(start,str):
            self.start = datetime.strptime(start,dfmt)
        else:
            sys.stdout.write('Start time must be datetime object or formatted string')


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
            quit()

        #check input wavelength formatting
        #check formatting assuming float or int
        if isinstance(wav,(int,float)):
            self.wav  = [int(wav)*u.AA]
            #check to make sure wavelength is allowed
            if int(self.wav[0].value) not in self.awavs:
                sys.stdout.write('{0:3.0f} not an acceptable wavelength'.format(self.wav.value))
                quit()

        #check formatting assuming string      
        elif isinstance(wav,str):
            self.wav = [int(wav)*u.AA]
            #check to make sure wavelength is allowed
            if int(self.wav[0].value) not in self.awavs:
                sys.stdout.write('{0:3.0f} not an acceptable wavelength'.format(self.wav.value))
                quit()

        #check formatting assuming array
        elif isinstance(wav,(list,np.ndarray)):
            self.wav = []
            for i in wav:
                if isinstance(i,(float,int)):
                    self.wav.append(int(i)*u.AA)
                elif isinstance(i,str):
                    self.wav.append(int(i)*u.AA)
                #check to make sure wavelength is allowed
                if int(self.wav[-1].value) not in self.awavs:
                    sys.stdout.write('{0:3.0f} not an acceptable wavelength'.format(i.value))
                    quit()


        #format input wavelength
        if isinstance(wav,list):
            self.wav = [ int(i)*u.AA for i in wav]
        elif isinstance(wav,(str,int)):
            self.wav = [int(wav)*u.AA]

    def get_drms_files(self):
        """
        Downloads the requested data from class object using drms (i.e. JSOC).
        """
        import drms
        client = drms.Client(email=self.email,verbose=False)
        fmt = '%Y.%m.%d_%H:%M'
        self.t_qstr = self.series+'[{0}_TAI-{1}_TAI@{2}]'.format(self.start.strftime(fmt),self.end.strftime(fmt),self.cadence) 


        #create wavelength query string
        self.w_qstr = '[' 
        for i in self.wav: self.w_qstr = self.w_qstr+'{0},'.format(int(i.value))
        #remove last , and add bracket
        self.w_qstr = self.w_qstr[:-1]+']'
        
        #make the series string
        self.s_qstr = '{'+self.segment+'}'

        #the full query
        self.qstr = self.t_qstr+self.w_qstr+self.s_qstr

        #IF ERRORS WITH URL ERROR IT IS BECAUSE THE DOWNLOAD FILE SIZE IS TOO LARGE
        #export  the data file list 
        self.expt = client.export(self.qstr)
        #create an array of indexes to download
        index = np.arange(np.size(self.expt.urls.url))

        #get output file names to check if file already exists
        outf = self.expt.urls.record.astype(str).str.replace(':','').str.replace('{','.').str.replace('}','.').str.replace('[','.').str.replace(']','.').str.replace('-','').str.replace('\.\.','.')+'fits'
        #Find if file exits is not then set check file to true so it keeps index
        check_file = [os.path.isfile(self.odir+i) == False for i in outf]
   
        #removed indices of already downloaded files
        index = index[check_file]

        # get new files from JSOC
        #set directory to current if no path set
        outf = self.expt.download(self.odir,index,fname_from_rec=True)
     

    def get_sunpy_files(self):
        """
        Downloads the requested data from class object using sunpy (i.e. JSOC).
        SUNPY IMPLEMENTATION DOES NOT WORK!
        """
        from sunpy.net import jsoc
        client = jsoc.JSOCClient()   
        dfmt = '%Y-%m-%dT%H:%M:%S' #date format for JSOC

        #get a query response from JSOC
        response = client.query(jsoc.Time(self.start.strftime(dfmt),self.end.strftime(dfmt)),jsoc.Series(self.series),
                                jsoc.Wavelength(self.wav),jsoc.Segment(self.segment),jsoc.Notify(self.email))

        res = client.get(response,path=self.odir,sleep=60,overwrite=self.overwrite,max_conn=self.max_con)
        res.wait(progress=True)

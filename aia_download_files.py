__version__ = "0.0.1 (2017/08/01)"
__authors__ = ['Jakub Prchlik <jakub.prchlik@cfa.harvard.edu>']
__email__   = "jakub.prchlik@cfa.harvard.edu"


#import dkrms
import astropy.units as u
import sys
from datetime import datetime
import numpy as np




class download_files:
    def __init__(self,start,end,wav,cadence,series='aia.lev1_euv_12s',segment='image',email=None,odir=None,
                 overwrite=True,max_con=1,dfmt='%Y/%m%/d %H:%M:%S'):
        """
        Will connect to aia_mkmovie
        """

        #list of acceptable wavelengths
        self.awavs = [94,131,171,193,211,304,335,1600,1700]
        #list of acceptable segments
        self.asegs = ['image','spike','None']
        #list of acceptable series
        self.asers = ['aia.lev1_uv_24s','aia.lev1_euv_12s']


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
            self.odir = odir
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
            if int(self.wav.value) not in self.awavs:
                sys.stdout.write('{0:3.0f} not an acceptable wavelength'.format(self.wav.value))
                quit()

        #check formatting assuming string      
        elif isinstance(wav,str):
            self.wav = [int(wav)*u.AA]
            #check to make sure wavelength is allowed
            if int(self.wav.value) not in self.awavs:
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
            self.wav = [i*u.AA]

    def get_drms_files(self):
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
# get file from JSOC
        #set directory to current if no path set
        if self.odir is None: self.odir = './'
        outf = self.expt.download(self.odir,index,fname_from_rec=True)
     

    def get_sunpy_files(self):
        from sunpy.net import jsoc
        client = jsoc.JSOCClient()   
        dfmt = '%Y-%m-%dT%H:%M:%S' #date format for JSOC

        #get a query response from JSOC
        response = client.query(jsoc.Time(self.start.strftime(dfmt),self.end.strftime(dfmt)),jsoc.Series(self.series),
                                jsoc.Wavelength(self.wav),jsoc.Segment(self.segment),jsoc.Notify(self.email))

        res = client.get(response,path=self.odir,sleep=60,overwrite=self.overwrite,max_conn=self.max_con)
        res.wait(progress=True)

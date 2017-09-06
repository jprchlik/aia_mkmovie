import ftplib
import numpy as np
import os
import sys, getopt
from datetime import timedelta,datetime




#Funtion to retrieve events from given day on noao ftp server
def getfiles(day,ftp,sdir):
    """
    Downloads GOES X-ray fluxes from SWPC archive

    Parameters
    ----------
    day  : datetime object
        Day to download from the swpc archive.
    ftp  : ftplib FTP object
        ftp connection to the swpc archive.
    sdir : string
        Output directory
    """
    files = '{0:%Y%m%d}_Gs_xr_1m.txt'.format(day)
#create file to write ftp data
    fhandle = open(sdir+'/goes/'+files,'wb')
    try:
        ftp.retrbinary('RETR {0}'.format(files),fhandle.write)
    except:
        print '{0} not in archive'.format(files)

    fhandle.close()

#Funtion to retrieve events from given day on noao ftp server
def getfiles_ace(day,ftp,sdir,swepam=False,mag=False):
    """
    Downloads ACE text files from NOAA SWPC archives for a given day.
   
    Parameters
    ----------
    day  : datetime object
        Day to download from the swpc archive.
    ftp  : ftplib FTP object
        ftp connection to the swpc archive.
    sdir : string
        Output directory
    swepam: boolean
        Download ACE solar wind parameters
    mag   : boolean
        Download ACE magnetometer
    """
    if swepam: files = '{0:%Y%m%d}_ace_swepam_1m.txt'.format(day)
    if mag: files = '{0:%Y%m%d}_ace_mag_1m.txt'.format(day)
#create file to write ftp data
    fhandle = open(sdir+'/ace/'+files,'wb')
    try:
        ftp.retrbinary('RETR {0}'.format(files),fhandle.write)
    except:
        print '{0} not in archive'.format(files)

    fhandle.close()


def look_xrays(start,end,sdir,email=None):
    """
    Function for downloading GOES flux data 
    from the NOAA archives for a given date range. 

    Parameters
    ----------
    start : datetime object
        The start time for downloading data from the NOAA swpc archives
    end   : datetime object
        The end time for downloading data from the NOAA swpc archives
    sdir  : string
        String to path of output directory
    email : string optional
        String of email address for connecting to the swpc archives       
    """
#initialize variables
    ftp = ftplib.FTP('ftp.swpc.noaa.gov','anonymous',email)
#change ftp directory to events directory for a given year
    ftp.cwd('/pub/lists/xray/')
    
    days =  np.arange(start,end,timedelta(days=1)).astype(datetime)
#loop through days     
    for day in days:
#request file from ftp server to add to local directory
        getfiles(day,ftp,sdir)
#change ftp directory to events directory for a given year
    ftp.cwd('../ace/')
    for day in days:
        getfiles_ace(day,ftp,sdir,swepam=True)
        getfiles_ace(day,ftp,sdir,mag=True)
    
#nicely leave the ftp server
    ftp.quit()


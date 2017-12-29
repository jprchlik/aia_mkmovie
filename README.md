aia_mkmovie
=======================

aia_mkmovie package creates time series movies from SDO/AIA images.  
The creates movie replicates movies created by aia_mkmovie in SSWIDL,
but all image processing comes from python packages.


SETUP
=====

    1. Install the following packages at or above the given levels. All listed
       packages may be installed via pip or conda with many included by default
       in Anaconda.
        a. matplotlib >= 2.0.0
        b. multiprocessing >= 0.70a1
        c. astropy >= 2.0.1
        d. sunpy >= 0.8.1
        e. Tkinter Revision: 81008
        f. pyfits >= 3.4
        g. re >= 2.2.1
        h. numpy >= 1.11.3
        i. drms >= 0.5.2 (not required)
    2. Download the lastest aia_mkmovie movie package.
        a. If you are download the package for the first time you have two
           dowlonad options.
            a1. From a terminal type the following:
                git clone https://github.com/jprchlik/aia_mkmovie.git
            OR
            a2. Navigate to https://github.com/jprchlik/aia_mkmovie in a
                webbrowser. Click clone or download. Select download zip.
                Unzip the file locally. 
    3. Install the aia_mkmovie package
        Change directory to the clone directory or the unzipped directory 
        (aia_mkmovie-master). Then type the following command if you want
        it install in the default location, which you probably do:

        python setup.py install

ABOUT THE CODE
==============

The code produces efficient SDO/AIA movies. A full wiki showing the features
available is located at the following link:
    https://github.com/jprchlik/aia_mkmovie/wiki  


Also, I included docstrings for every base module.
Therefore, understanding what a specific module does is as simple as typing help(*ModuleName*).
For example help(am.aia_mkmovie) in python.

CONTACT
=======

Please post bug reports, patches, and other feedback to 

    https://github.com/jprchlik/aia_mkmovie/issues



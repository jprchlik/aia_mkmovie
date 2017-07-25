#!/usr/bin/env python

__version__ = "0.4.1 (2013.01.13)"
__authors__ = ["Jonathan P. Sattelberger <jsattelb@head.cfa.harvard.edu>"]
__email__ = "jsattelb@head.cfa.harvard.edu"


import datetime
import glob
try:
	import numpy
except ImportError:
	import sys
	print >> sys.stderr, "numpy not installed, use \"pip install numpy --upgrade\""
	sys.exit(1)
import os
try:
	import pyfits
except ImportError:
	import sys
	print >> sys.stderr, "pyfits not installed, use \"pip install pyfits --upgrade\""
	sys.exit(1)
import re
import sys

#HDWIII 2016_09_02 Edited out AIAGP6 condition post eclipse
#'' Edited out #if FITS_KEYWORD__AIFTSID >= 49152: post eclipse.
""" SMEARpy - (S)DO (M)etadata (E)xtr(A)cto(R) for (Py)thon
-----------------------------------------------------------------------------------------
"""

class Scream(object):

	"""SCREAM - (SC)ience (R)eady fil(E)s for Data (A)nalysis (M)odule
------------------------------------------------------------------------
>>> from SMEARpy import Scream
>>> src = Scream('/data/SDO/AIA/synoptic')

#>>> src = Scream(archive = '/data/SDO/AIA/synoptic',
#	date = '2013-10-01', time = '00:00:00', span = '-3d',
#	wavelnth = 171, filetype='fits')

# return a list of files from the given paths and instrument passband

>>> fls = src.get_filelist(date = '2013-10-01', time = '00:00:00', span = '-3d', wavelnth = 171)

# return a list of times for the given filelist

>>> tms = src.get_filetimes() 
	
# return a list of files and times that pass the quality assurance tests

>>> qfls, qtms = src.run_quality_check(synoptic=True)

# return a list of files within a desired cadence, say 1 file per 3 
minutes

>>> hfls = src.get_sample(files = qfls, sample = '3m', nfiles = 1)
------------------------------------------------------------------------
	"""

	#def __init__(self, archive = '', date = None, time = '00:00:00', span = '24h', \
	#	wavelnth = '', instr = None, filetype = None, \
	#	verbose = False, debug = False):

	def __init__(self, archive = '', verbose = False, debug = False):
	
		"""Required keyword arguments:
--------------------------------------------------------------------
archive: Absolute path to the primary data repository
date: Date in YYYY-MM-DD format
		
Optional arguments:
--------------------------------------------------------------------
time: Time of observation in HH:MM:SS format. Default: '00:00:00'
span: Duration of observation. This may be specified in hours, minutes or seconds by 
appending h, m, or s to the end of your query. Default: '24h'
wavelength: Specify from 0094, 0131, 0171, 0193, 0211, 0304, 1600, 1700, 4500, or leave 
blank for all. Default: ''
		
Additional arguments:
--------------------------------------------------------------------
instr: Observational instrument. Might be useful for querying for HMI as well. Default: 
'AIA'
filetype: Can be changed to other filetype associations, say 'png'. Sampling of AIA data 
requires 'fits' filetype. Default: 'fits'
		
Debugging booleans:
--------------------------------------------------------------------
verbose: Print some stuff. Default: False
debug: print even more stuff. Default: False
		"""
	
		# user assigned variables
		################################################################################
		
		self.archive = archive
		self.date = None # date
		self.time = None # time
		self.span = None # span
		
		self.wavelnth = None # wavelnth
		self.filetype = None # 'fits'
		
		#self.pattern = '%s*%s.%s' % (self.instr, self.wavelnth, self.filetype)
			
		# private variables
		################################################################################

		self._paths = None
		self._filelist = None
		self._times = None
		
		# debug options
		################################################################################
		
		self.verbose = verbose
		self.debug = debug
		
	# private: fetch timedelta from a given string-span designation
	def _get_spca(self, span = None):
	
		"""Convert a given string to datetime.timedelta. Returns datetime.timedelta and 
the extracted regular expression groups"""
	
		span = self.span if (span is None) else span

		prog = re.compile(r"^(?P<pos>[-+]|)?(?P<num>\d+)(?P<span>[dhms])$")
		m = re.match(prog, span)
		print span,m
		
		# Allowed datetime data-types
		if m.group('span') == 'd':
			td = datetime.timedelta( days = long( m.group('num') ) )
		elif m.group('span') == 'h':
			td = datetime.timedelta( hours = long( m.group('num') ) )
		elif m.group('span') == 'm':
			td = datetime.timedelta( minutes = long( m.group('num') ) )
		elif m.group('span') == 's':
			td = datetime.timedelta( seconds = long( m.group('num') ) )	
		else:
			if self.verbose == True or self.debug == True:
				print >> sys.stderr, "Unable to extract a valid timedelta from span - %s. Setting span to 0 seconds." % span
			return datetime.timedelta( seconds = 0 ), m
			
		return td, m
	
	# private: fetch datetime range and timedelta
	def _get_dtr(self, date = None, time = None, span = None):
	
		"""Get start-end datetime range and timedelta. Returns a 3 elements: start time (datetime), end time (datetime), and a time difference (datetime.timedelta)."""
		
		date = self.date if (date is None) else date
		time = self.time if (time is None) else time
		span = self.span if (span is None) else span

		try:
			dt = datetime.datetime.strptime(date +" "+ time, "%Y-%m-%d %H:%M:%S")
		except ValueError as e:
			print >> sys.stderr, "Unable to extract date-time -", e
			return None, None, None
		
		td, expr = self._get_spca(span)
		
		if expr.group('pos') == "-":
			return (dt - td), dt, td
	
		return dt, (dt + td), td
	
	# fetch file paths for the desired date, time and span.
	def get_paths(self, date = None, time = None, span = None):
	
		"""Fetch paths based on the date, time and span. Will only return paths which exist."""
		
		date = self.date if (date is None) else date
		time = self.time if (time is None) else time
		span = self.span if (span is None) else span
		
		self._paths = []

		# fetch the datetime range: start, end, and time difference
		dt, dttd, td = self._get_dtr(date = date, time = time, span = span)
		
		# unable to fetch a list of files with the given date, time or span
		if dt is None:
			return []
		
		# Archive path generation for the given datetime span
		start = datetime.datetime(dt.year, dt.month, dt.day)
		end = datetime.datetime(dttd.year, dttd.month, dttd.day)
			
		i = datetime.timedelta(days = 1)
		while start <= end:
			# Ensure the directory exists before appending.
			if os.path.isdir(self.archive +'/'+ start.strftime('%Y/%m/%d')) == True:
				self._paths.append(self.archive +'/'+ start.strftime('%Y/%m/%d'))
				if self.verbose == True or self.debug == True:
					print "Adding date: %s" % (self.archive +'/'+ start.strftime('%Y/%m/%d'))
			start += i 
		return self._paths
	
	# fetch filelist based on the desired filetype and passband
	def get_filelist(self, paths = None, date = None, time = None, span = None, wavelnth = None, instr = None, filetype = None):
		
		"""Fetch a list of files based on the desired passband, date, time, span, and search expression."""
		
		date = self.date if (date is None) else date
		self.date = date
		time = self.time if (time is None) else time
		self.time = time
		span = self.span if (span is None) else span
		self.span = span	
			
		paths = self._paths if (paths is None) else paths
		paths = self.get_paths(date = date, time = time, span = span) if (paths is None) else paths
		
		wavelnth = str(self.wavelnth).zfill(4) if (wavelnth is None) else str(wavelnth).zfill(4)
		
		# AIA is the default observable instrument
		instr = 'AIA' if instr is None else instr
		
		# Look for FITS files, but it may also be useful to filter by other filetypes, PNG, JPEG, TXT, etc.
		filetype = 'fits' if filetype is None else filetype
		self.filetype = filetype

		# glob expression
		pattern = '%s*%s.%s' % (instr, wavelnth, filetype)
		
		filelist = []
		
		# fetch the datetime range: start, end, and time difference
		dt, dttd, td = self._get_dtr(date = date, time = time, span = span)
	
		# fetch list of files
		for p in paths:
			dirs = os.listdir(p)
			for d in dirs:
				files = glob.glob(os.path.join(p, d, pattern))
				for f in files:
					filelist.append(f)
		
		# extract date-time for each filename in the list
		filelist = sorted(filelist) # sort filelist
		times = self.get_times(files = filelist) # fetch a list of times
		
		# fetch the indicies from the _times list falling between the datetime range
		indxs = [indx for indx, time in enumerate(times) if datetime.datetime.strptime(time, '%Y%m%d%H%M%S') >= dt \
			and datetime.datetime.strptime(time, '%Y%m%d%H%M%S') <= dttd]
			
		# filter the file list based on the valid time indicies, store
		self._filelist = [ filelist[indx] for indx in indxs ]
		
		# filter the time list based on the valid time indicies, store
		self._times = [ times[indx] for indx in indxs ]		
		
		return self._filelist # return to user
	
	# fits file time array
	def get_times(self, files = None):
		
		"""Function used to return a list of times based on the list of files from the desired passband"""
		
		files = self._filelist if (files is None) else files
		files = self.get_filelist() if (files is None) else files
		
		# allocate an array of length filelist, with the default datetime-string 1776-07-04 01:23:54
		times = len(files)*['17760704012354']
		
		# extract date-time for each filename in the list
		for i, f in enumerate(files):
			try:
				fn = f.split('/')[-1]
				ft = "%s%s%s%s%s%s" % (fn[3:7], fn[7:9], fn[9:11], fn[12:14], fn[14:16], fn[16:18])
				if self.verbose == True or self.debug == True:
						print "%s -> %s" % (f, ft)
				# successfully parsed filename for date and time 
				times[i] = ft
			except ValueError:
				print >> sys.stderr, "%s - Unable to extract date and time from filename" % f 
				continue
		
		return times
		
	# wrapper for get_times()
	def get_filetimes(self, files = None):
	
		"""Function used to return a list of times based on the list of files from the desired passband"""

		files = self._filelist if (files is None) else files	
		return self.get_times(files = files)
		
	# return a list of files and times that pass the synoptic quality check
	def run_quality_check(self, files = None, times = None, filetype = None, synoptic = False):
	
		"""Return a list of FITS files and FITS file times that have passed a basic file analysis."""
	
		# make sure only FITS files can be analyzed for data quality integrity
		filetype = self.filetype if filetype is None else filetype
		if not isinstance(filetype, str):
			print >> sys.stderr, "Filetype association was not set during initalization or it has been reset."
			return [], []
		
		if re.search('^FITS$', filetype, flags=re.IGNORECASE) is None:
			print >> sys.stderr, "FITS filetype is expected. Unable to perform quality check on the provided filetype."
			return [], []
		
		files = self._filelist if files is None else files
		times = self.get_times(files) if times is None else times	

		# Search for new files to process. Assume all files fail the quality check
		rlts = len(files)*[False]
		BitArray = numpy.uint64(2)**numpy.arange(0, 32, 1, dtype=numpy.uint64)

		for i, f in enumerate(files):
			rlts[i] = self.fits_quality_check(f, BitArray = BitArray, synoptic = synoptic)

		indxs = [indx for indx, rslt in enumerate(rlts) if rslt == True ]
		files = [ files[indx] for indx in indxs ]
		times = [ times[indx] for indx in indxs ]
		
		# return a list of files and times that have passed the quality check
		return files, times
	
	# fits file quality assumption		
	def fits_quality_check(self, file = '', BitArray = None, synoptic = False):

		"""Fits quality chekcer level - Acceptable files return True. Unacceptable files return False."""
	
		### HEADER SECTION ###
	
		fits_file = file
		try:
			hdulist = pyfits.open(fits_file)
		except IOError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - Unable to open file." % fits_file
			if self.debug == True: print >> sys.stderr, e
			return False # Reject Image
		
		# Verify the HDU, fix if required. Return false if verficiation or fix fails.
		try:
			hdulist.verify('silentfix')
		except: 
			print >> sys.stderr, '%s - Exception handling HDU verfication' % (fits_file)
			hdulist.close()
			return False
		
		try:
			prihdr = hdulist[0].header # Primary FITS header (Standard)
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Primary Header Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False # Reject Image

		try:
			if synoptic == True:
				sechdr = hdulist[0].header # secondary FITS header does not exist for synoptic data
			else:
				sechdr = hdulist[1].header # Secondary FITS header (AIA Specific)
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Secondary Header Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False # Reject Image
		
		### DATA SECTION ###
		
		# Data cube size quality check on the FITS file data using PyFITS
		try:
			if synoptic == True:
				(y, x) = hdulist[0].data.shape
				sz = hdulist[0].data.size	
			else:
				(y, x) = hdulist[1].data.shape
				sz = hdulist[1].data.size
		except:
			print >> sys.stderr, '%s - Unable to fetch data dimensions' % (fits_file)
			hdulist.close()
			return False
		
		if synoptic == True:
			if x != 1024 or y != 1024:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Array Dimensions do not equal 1024 x 1024 - (%s, %s)." % (fits_file, x, y)
				if self.debug == True: print >> sys.stderr, e
				hdulist.close()
				return False # Reject Image	
		else:
			if x != 4096 or y != 4096:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Array Dimensions do not equal 4096 x 4096 - (%s, %s)." % (fits_file, x, y)
				if self.debug == True: print >> sys.stderr, e
				hdulist.close()
				return False # Reject Image	

		if synoptic == True:
			if sz != 1048576:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Array elements does not equal 1048576 - (%s)." % (fits_file, sz)
				if self.debug == True: print >> sys.stderr, e
				hdulist.close()
				return False # Reject Image
		else:
			if sz != 16777216:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Array elements does not equal 16777216 - (%s)." % (fits_file, sz)
				if self.debug == True: print >> sys.stderr, e
				hdulist.close()
				return False # Reject Image
		
		### KEYWORD SECTION ###	
		
		# Header 'QUALITY' keyword test phase
		try:
			FITS_KEYWORD__QUALITY = sechdr['QUALITY']
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False

        #Remove flare short exposures (J. Prchlik 2016/12/05)
		try:
			FITS_KEYWORD__EXP_TIME = sechdr['EXPTIME']

			if FITS_KEYWORD__EXP_TIME <= 1.85: #seconds
				print >> sys.stderr,"%s - Image Quality Result: Failed - EXP_TIME less than 1.85s" % fits_file
				hdulist.close()
				return False
		except KeyError as e:
				print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
				if self.debug == True: print >> sys.stderr, e
				hdulist.close()
				return False
            
		
		# bitwise_and for QUALITY keyword
		if synoptic == True:
			# lets not worry about the bitwise check on synoptic data. It will likely be 
			# caught during the other FITS keyword checks.
			pass
		else:
			# non-synoptic data bitwise testing
			try:
				# BitArray can now created before the function call. When BitArray is not passed to
				# the function (i.e. type is None), will create the required BitArray reference.
				if BitArray is None:
					BitArray = numpy.uint64(2)**numpy.arange(0, 32, 1, dtype=numpy.uint64)
				
				BitSet = numpy.bitwise_and(FITS_KEYWORD__QUALITY, BitArray) # True, False
		
				ForbiddenBits = [0,1,2,3,4,12,13,14,15,16,17,18,20,21,31] # If any of these bits are set - reject the image
		
				if BitSet[ForbiddenBits].sum() > 0:
					print >> sys.stderr, "%s - Image Quality Result: Failed - Failed Bitwise Test." % fits_file
					hdulist.close()
					return False # Reject Image
			except TypeError as e:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Unable to Process BitSet (QUALITY %s = %s)." % \
					(fits_file, type(FITS_KEYWORD__QUALITY), FITS_KEYWORD__QUALITY)
				hdulist.close()
				if self.debug == True: print >> sys.stderr, e
				return False

		try:
			FITS_KEYWORD__PERCENTD = sechdr['percentd']
			if FITS_KEYWORD__PERCENTD < 99.9999:
				print >> sys.stderr, "%s - Image Quality Result: Failed - PERCENTD < 99.9999%%." % fits_file
				hdulist.close()
				return False # Reject Image
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False

		# Venus Transit / Pixel Count
		try:
			FITS_KEYWORD__DATAVALS = sechdr['datavals']
		
			if FITS_KEYWORD__DATAVALS < 16777200:
				print >> sys.stderr, "%s - Image Quality Result: Failed - DATAVAL < 16777200." % fits_file
				hdulist.close()
				return False # Reject Image
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False
		
		# AIFTSID - Currently used to denote types of observing modes. Anything greater than or equal to 49152
		# indicates a calibration image that should be rejected.
		try:
			FITS_KEYWORD__AIFTSID = sechdr['aiftsid']
		
			if FITS_KEYWORD__AIFTSID >= 49152:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Calibration Image Detected." % fits_file
				hdulist.close()
				return False # Reject Image
		
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False
		
		# Open Filters
		try:
			FITS_KEYWORD__WAVE_STR = sechdr['WAVE_STR']
		
			if re.search('OPEN', FITS_KEYWORD__WAVE_STR, flags=re.IGNORECASE) is not None:
				print >> sys.stderr, "%s - Image Quality Result: Failed - Open EUV Filter Mode Detected." % fits_file
				hdulist.close()
				return False # Reject Image
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False

                 
  
		
		# Eclipses and AIAGP6
		try:
			FITS_KEYWORD__ACS_ECLP = sechdr['ACS_ECLP']
		
			if FITS_KEYWORD__ACS_ECLP == "YES":
				try:
					FITS_KEYWORD__AIAGP6 = sechdr['AIAGP6']
				
					if FITS_KEYWORD__AIAGP6 != 0:
						print >> sys.stderr, "%s - Image Quality Result: Failed - AIAGP6 != 0" % fits_file
						hdulist.close()
						return False # Reject Image
				except KeyError as e:
					print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
					if self.debug == True: print >> sys.stderr, e
					hdulist.close()
					return False
		
		except KeyError as e:
			print >> sys.stderr, "%s - Image Quality Result: Failed - FITS Keyword Not Found." % fits_file
			if self.debug == True: print >> sys.stderr, e
			hdulist.close()
			return False
			
		### POSITIVE END RESULT ###
		if self.verbose == True or self.debug == True:
			print "%s - Image Quality Result: Passed" % fits_file
		
		hdulist.close() # remember to close the hdu
	
		return True # Accept Image
		
	# based on the IDL routine aia_sample.pro by Dr. Alisdair Davey <ard@head.cfa.harvard.edu>
	def get_sample(self, files = None, sample = '12s', nfiles = 1):
		"""
# files - list of filenames to sample
# sample - interval of sampling using ?(\d+)([h|m|s]) 
# selection - list of sampled filenames
# nfiles - no. of sample files in each interval, default = 1
"""
		
		if not isinstance(files, list):
			print >> sys.stderr, "list of files required"
			return None
			
		# fetch a list of times based on the provided file list
		times = self.get_times(files)
		
		# returns a float containing the total number of seconds in the timedelta. Not to be confused with timedelta.seconds
		hsample = (self._get_spca(sample)[0]).total_seconds()

		# an array of floats containing the total number of seconds in the timeselta per the list item
		basetime = times[0] # baseline time sample
		htimes = [(datetime.datetime.strptime(times[indx], '%Y%m%d%H%M%S') - datetime.datetime.strptime(basetime, '%Y%m%d%H%M%S')).total_seconds() \
			for indx, time in enumerate(times)]
		nphtimes = numpy.array(htimes) # make sure to convert the python-list to a numpy-array object
		
		# the numpy histrogram does not behave like the IDL routine. We need to determine 
		# the length maximum (hrange) value for the histogram, and create an array (hbins) with evenly 
		# spaced bins (hsample) from 0 to N (hrange), plus an additional bin (hsample).
		
		hrange = (htimes[-1] - htimes[0]) # determine the maximum length of the resultant htimes array
		# create a numpy range from 0 -> N + hsample with a step-size of hsample (1, 2, ..., 300, 600, N)
		hbins = numpy.arange(0, hrange + hsample, hsample) # add an additional sampling point to encapsulate any items on the fringe
		#for i in numpy.nditer(hbins):
		#	print i,
		hist, bins = numpy.histogram(nphtimes, bins = hbins) # histogram function
		
		selection = len(times)*[''] # list of acceptable files, empty strings will be filtered
		
		scounter = 0 # file selection ind. counter
		hloc = 0 # histogram bin width pos.
		sloc = 0 # skipped histogram bin width pos.
		
		# iterate over the numpy-array
		for hindx in numpy.nditer(hist):
			# to exclude or not to exclude
			if hindx < nfiles: 
				# was the size of the bin less than expected? add the bin-size to the skipped bins counter
				if self.verbose == True or self.debug == True:
					for f in files[(hloc + sloc) : (hloc + sloc) + hindx]:
						print >> sys.stderr, "%s - Skipped file. Calculated binning: %d. Required binning: %d. " % (f, hindx, nfiles)
				
				sloc += hindx # excluded histogram bin width
				continue # next!
							
			# array lookup includes x -> y[-1], the last value, y is not included unless it is the last value
			selection[ scounter : scounter + nfiles] = files[ (hloc + sloc) : (hloc + sloc) + nfiles ]

			hloc += hindx # increment the bin counter
			scounter += nfiles # increment the selection counter
			
		if sloc > 0:
			if self.verbose == True or self.debug == true:
				print >> sys.stderr, "Skipped %d files." % sloc
		
		indxs = [ indx for indx, select in enumerate(selection) if select != '' ]
		selection = [ selection[indx] for indx in indxs ] # filter list to exclude empty string items
		
		return selection # return a list of valid files

# END

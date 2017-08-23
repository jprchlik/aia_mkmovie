import os
import subprocess
import shutil
import numpy as np
import glob
from multiprocessing import Pool

class create_movie:

    def __init__(self,odir='',pdir='track_plots/',ext='png',nproc=1,w0=1024,h0=1024,frate=10,outmov='movie.mp4'):

        self.startd = os.getcwd()
        #concert pdir into absolute path
        pdir = os.path.abspath(pdir)
        if pdir[-1] != '/': pdir = pdir+'/'
        if odir == '': odir = os.getcwd()
        if odir[-1] != '/': odir = odir+'/'
        if outmov[-3:] != 'mp4': outmov = outmov+'.mp4'
        if ext[0] == '.' : ext = ext[1:]
 # initalize all variables
        self.w0 = w0
        self.h0 = h0
        self.frate = frate
        self.outmov = outmov  
        self.odir = odir
        self.pdir = pdir
        self.nproc= nproc
        self.ext  = ext
        self.sdir = pdir+'symlinks'
#remove the previous direcotry with symlinks
        shutil.rmtree(self.sdir)
#readd symlinks directory
        os.mkdir(self.sdir)
        
    

#gather files for creating symbolic links
    def gather_files(self):
        self.files = glob.glob('{0}*.{1}'.format(self.pdir,self.ext))
        try:
            self.lengs = str(int(1+np.log10(len(self.files)))) # get the number of sig figs for file format
        except OverflowError:
            print "No {0} files found in {1}".format(self.ext,self.pdir)

        self.sfmt  = '{0:'+self.lengs+'d}.'+self.ext #Create formating for file out index

        index = range(len(self.files))#array index for files

#number of processors greater than 1 do in parallel
#        if self.nproc > 1:
#            pool = Pool(processes=self.nproc)
#            pool.map(self.create_links,index)
#            pool.close()
#        if self.nproc == 1:
        for j in index: self.create_links(j)
#make sure nproc is greater than or equal to 1
#        if self.nproc < 1: 
#            print 'Need to have at least 1 processor allocated'
#            raise   
#use index to create formated symbolic links numerically increasing
    def create_links(self,index):
        ifile = self.files[index]
        os.symlink(ifile,'{0}/seq'.format(self.sdir)+self.sfmt.format(index).replace(' ','0'))


#write ffmpeg to file
    def write_ffmpeg(self):
        com = open(self.odir+'run_ffmpeg.csh','w')
        com.write('/usr/local/bin/ffmpeg -y -f image2 -r {2:2d} -i {5} -an -pix_fmt "yuv420p" -vcodec libx264 -level 41 -crf 18.0 -b 8192k -r {2:2d} -bufsize 8192k -maxrate 8192k -g 25 -coder 1 -profile main -preset faster -qdiff 4 -qcomp 0.7 -directpred 3 -flags +loop+mv4 -cmp +chroma -partitions +parti4x4+partp8x8+partb8x8 -subq 7 -me_range 16 -keyint_min 1 -sc_threshold 40 -i_qfactor 0.71 -rc_eq "blurCplx^(1-qComp)" -s "{0:4d}x{1:4d}" -b_strategy 1 -bidir_refine 1 -refs 6 -deblockalpha 0 -deblockbeta 0 -trellis 1 -x264opts keyint=25:min-keyint=1:bframes=1 -threads {4:1d} {3}\n'.format(int(self.w0),int(self.h0),int(self.frate),self.outmov,int(self.nproc),'../working/symlinks/seq%'+self.lengs+'d.'+self.ext))
        com.close()
#Actually run ffmpeg
    def run_ffmpeg(self):
        #change to directory with ffmpeg written
        os.chdir(self.odir)
       
        #change file to executable
        mod = subprocess.call(['/bin/csh','-c','chmod a+x run_ffmpeg.csh'])
        #run ffmpeg
        run = subprocess.call(['/bin/csh','-c','./run_ffmpeg.csh'])
        #rename the file to the file directory (cuts down on down time if running in a loop)

    def create_movie(self):
        self.gather_files()
        self.write_ffmpeg()
        self.run_ffmpeg()
        os.chdir(self.startd)
        
         
        

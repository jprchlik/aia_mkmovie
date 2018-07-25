import matplotlib
#Use TkAgg backend for plotting 
matplotlib.use('TkAgg',warn=False,force=True)


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#implement the deault mpl key bindings
from matplotlib.backend_bases import key_press_handler,MouseEvent
import numpy as np
import sys
import matplotlib.pyplot as plt
from datetime import datetime

try:                                                                                      
    import sunpy.map                                                                      
    from sunpy.cm import cm                                                               
except ImportError:                                                                       
    sys.stdout.write("sunpy not installed, use pip install sunpy --upgrade")   

#check the python version to use one Tkinter syntax or another
if sys.version_info[0] < 3:
    import Tkinter as Tk
    import tkMessageBox as box
    import tkFileDialog as Tkf
else:
    import tkinter as Tk
    from tkinter import messagebox as box
    from tkinter import filedialog as Tkf


#main gui class
class gui_c(Tk.Frame):

    #init gui class
    def __init__(self,parent,flist,w0=1900,h0=1144,cy=0,cx=0,color3=False,img_scale=None):
        """
        Shows AIA image in 3 color or single image for scaling and region selection. After you
        choose a region, the class will pass the parameters onto aia_mkimage from aia_mkmovie.

        Parameters
        ----------
        parent : Tk.Frame
            A Tk.Frame instance.
        flist  : list
            A list of files. If the list has 3 dimensions then the GUI will create a 3 color 
            image.
        w0: int or float, optional
            The pixel width of the output to pass to aia_mkmovie (Can be set in the GUI).
            If the height (h0) is larger than
            w0 the program will switch the two parameters on output. 
            However, it will also transpose the x and y axes, which allows 
            for rotated images and movies. Default = 1900
        h0: int or float, optional 
            The pixel height of the output to pass to aia_mkmovie (can be set in GUI).
            If h0 is larger than the
            width (w0) the program will switch the two parameters on
            output. However, it will also transpose the x and y axes,
            which allows for rotated images and movies. Default = 1144
        cx  : float or int, optional
            Center of the field of view for creating images. If cx is set
            then the image is assumed to be a cutout. Selecting in prompt
            overrides cx. Default = 0.0 (can be set in GUI).
        cy  : float or int, optional
            Center of the field of view for creating images. If cy is set
            then the image is assumed to be a cutout. Selecting in prompt
            overrides cy. Default = 0.0 (can be set in GUI).
        color3 : boolean, optional
            Create a 3 color image. If color3 set to True panel must be
            False and the wavelength list must be 4 wavelengths long.
            The wav list has the following format [R, G, B]. Default =
            False.
        img_scale: dictionary, optional
            Pass a dictionary where the key is a 4 character wavelength string with left padded 0s
            in Angstroms and the values are a list. The first element in the list is a color map. 
            By default the first element contains the color map given by sunpy for a given wavelength
            (e.g. for 131 the color map is cm.sdoaia131). The second and third element are respectively
            the minimum and maximum color map values. The minimum and maximum assume a arcsinh
            transformation and exposure normalized values. The program uses arcsinh for all image
            scaling because the arcsinh function behaves like a log transformation at large 
            values but does not error at negative values. If the user gives no image scale
            then a default image scale loads. The default color table works well for single
            and panel images but not for 3 color images. Will pass updated values updated 
            in GUI to aia_mkimage through aia_mkmovie.
        """
        Tk.Frame.__init__(self,parent,background='white') #create initial frame with white background

        #set the starting list to be 0
        self.order = 0
        #get maximum value in list
        self.ordermax = len(flist)-1


        #create parent variable
        self.parent = parent


        #check image height
        if isinstance(h0,(int,float)):
            self.h0 = h0
        else:
            sys.stdout.write('h0 must be an integer or float (Assuming 0)')
            self.h0 = 0.0

        #check image width
        if isinstance(w0,(int,float)):
            self.w0 = w0
        else:
            sys.stdout.write('w0 must be an integer or float (Assuming 0)')
            self.w0 = 0.0


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







        #3 color images
        if isinstance(color3,bool):
            self.color3 = color3
        else:
            sys.stdout.write('color3 must be a boolean')



        #file list
        if isinstance(flist,str):
            self.flist = [flist]
        elif isinstance(flist,list):
            self.flist = flist
            #correct for single value list  in color3
            if ((self.color3) & (len(flist) == 3)): self.flist = [self.flist]
        else:
            sys.stdout.write('flist must be a string or list')
      

        #clicked point on the figure
        self.clicked = False
        #first clicked point on the figure
        self.firstclick = False


        #Dictionary for vmax, vmin, and color
        if img_scale is None:
            self.img_scale = {'0094':[cm.sdoaia94  ,np.arcsinh(1.),np.arcsinh(150.)],
                              '0131':[cm.sdoaia131 ,np.arcsinh(1.),np.arcsinh(500.)],
                              '0171':[cm.sdoaia171 ,np.arcsinh(10.),np.arcsinh(2500.)],
                              '0193':[cm.sdoaia193 ,np.arcsinh(10.),np.arcsinh(4500.)],
                              '0211':[cm.sdoaia211 ,np.arcsinh(10.),np.arcsinh(4000.)],
                              '0304':[cm.sdoaia304 ,np.arcsinh(2.),np.arcsinh(300.)],
                              '0335':[cm.sdoaia335 ,np.arcsinh(1.),np.arcsinh(100.)],
                              '1600':[cm.sdoaia1600,np.arcsinh(20.),np.arcsinh(500.)],
                              '1700':[cm.sdoaia1700,np.arcsinh(200.),np.arcsinh(4000.)]}
        elif isinstance(img_scale,dict):
            self.img_scale = img_scale
        else:
            sys.stdout.write('img_scale must be a dictionary with color map, min value, max value')



#Start the creation of the window and GUI
        self.centerWindow()
        self.FigureWindow()
        self.initUI()
        self.aia_set()
        self.aia_plot()

#Create area and window for figure
    def FigureWindow(self):
#set the information based on screen size
        x =  self.parent.winfo_screenwidth()
        y =  self.parent.winfo_screenheight()


        aiaframe = Tk.Frame(self)

        aratio = float(x)/float(y)
#Create the figure
        self.f,self.a = plt.subplots(ncols=2,figsize=(8*aratio,8*aratio*.5))
#Separate the two plotting windows
        self.x = self.a[1]
        self.a = self.a[0]
        #turn off clicked axis for starters
        self.x.set_axis_off()
 
#Create window for the plot
        self.canvas = FigureCanvasTkAgg(self.f,master=self)
#Draw the plot
        self.canvas.draw()
#Turn on matplotlib widgets
        self.canvas.get_tk_widget().pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
#Display matplotlib widgets
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
#Connect mpl to mouse clicking
        self.f.canvas.mpl_connect('button_press_event',self.on_click_event)
#Connect mpl to mouse clicking
        #self.f.canvas.mpl_connect('key_press_event',self.on_key_event)

#create button to go up an order
        upbutton = Tk.Button(master=aiaframe,text='Increase File',command=self.increaseorder)
        upbutton.pack(side=Tk.LEFT)
#create button to go down an order
        downbutton = Tk.Button(master=aiaframe,text='Decrease File',command=self.decreaseorder)
        downbutton.pack(side=Tk.LEFT)


        aiaframe.pack(side=Tk.TOP)


#Create window in center of screen
    def centerWindow(self):
        w = 2000
        h = 1200
        sw = self.parent.winfo_screenwidth()
        sh = self.parent.winfo_screenheight()
       
        x = (sw-w)/2
        y = (sh-h)/2
        self.parent.geometry('%dx%d+%d+%d' % (w,h,x,y))


#Intialize the GUI
    def initUI(self):
#set up the title 
        self.parent.title("Select AIA Region")

#create frame for plotting
        frame = Tk.Frame(self,relief=Tk.RAISED,borderwidth=1)
        frame.pack(fill=Tk.BOTH,expand=1)
 
        self.pack(fill=Tk.BOTH,expand=1)

#set up okay and quit buttons
        quitButton = Tk.Button(self,text="Quit",command=self.onExit)
        quitButton.pack(side=Tk.RIGHT,padx=5,pady=5)

#set up center width box
        w0Text = Tk.StringVar()
        w0Text.set("Width (pixels)")
        w0Dir = Tk.Label(self,textvariable=w0Text,height=4)
        w0Dir.pack(side=Tk.LEFT)
#Add so width can be updated
        w0 = Tk.StringVar()
        w0.set('{0:5.2f}'.format(self.w0))
        self.w0val = Tk.Entry(self,textvariable=w0,width=10)
        self.w0val.bind("<Return>",self.aia_param)
        self.w0val.pack(side=Tk.LEFT,padx=5,pady=5)
       
#set up center h0 box
        h0Text = Tk.StringVar()
        h0Text.set("Height (pixels)")
        h0Dir = Tk.Label(self,textvariable=h0Text,height=4)
        h0Dir.pack(side=Tk.LEFT)
#Add so center height can be updated
        h0 = Tk.StringVar()
        h0.set('{0:5.2f}'.format(self.h0))
        self.h0val = Tk.Entry(self,textvariable=h0,width=10)
        self.h0val.bind("<Return>",self.aia_param)
        self.h0val.pack(side=Tk.LEFT,padx=5,pady=5)
       
   
#set up center x0 box
        cxText = Tk.StringVar()
        cxText.set("X0 (arcsec)")
        cxDir = Tk.Label(self,textvariable=cxText,height=4)
        cxDir.pack(side=Tk.LEFT)
#Add so center x can be updated
        self.scx = Tk.StringVar()
        self.scx.set('{0:5.2f}'.format(self.cx))
        self.cxval = Tk.Entry(self,textvariable=self.scx,width=10)
        self.cxval.bind("<Return>",self.aia_param)
        self.cxval.pack(side=Tk.LEFT,padx=5,pady=5)
       
#set up center y box
        cyText = Tk.StringVar()
        cyText.set("Y0 (arcsec)")
        cyDir = Tk.Label(self,textvariable=cyText,height=4)
        cyDir.pack(side=Tk.LEFT)
#Add so center y can be updated
        self.scy = Tk.StringVar()
        self.scy.set('{0:5.2f}'.format(self.cy))
        self.cyval = Tk.Entry(self,textvariable=self.scy,width=10)
        self.cyval.bind("<Return>",self.aia_param)
        self.cyval.pack(side=Tk.LEFT,padx=5,pady=5)

#set up order number
        orderText = Tk.StringVar()
        orderText.set("Order")
        orderDir = Tk.Label(self,textvariable=orderText,height=4)
        orderDir.pack(side=Tk.LEFT)
#Add so order number can be updated
        self.sorder = Tk.StringVar()
        self.sorder.set(str(int(self.order)))
        self.orderval = Tk.Entry(self,textvariable=self.sorder,width=5)
        self.orderval.bind("<Return>",self.on_order_box)
        self.orderval.pack(side=Tk.LEFT,padx=5,pady=5)
       


    #boxes to create if 3 color image
        if self.color3:  
###############################################
#      BLUE COLOR BOXES                       #
###############################################
    #Add so Color Min can be updated
            self.bcmin = Tk.StringVar()
            self.bcmin.set('{0:5.2f}'.format(0))
            self.bcminval = Tk.Entry(self,textvariable=self.bcmin,width=10)
            self.bcminval.bind("<Return>",self.aia_param)
            self.bcminval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Min
            bcminText = Tk.StringVar()
            bcminText.set("B Color Min.")
            bcminDir = Tk.Label(self,textvariable=bcminText,height=4)
            bcminDir.pack(side=Tk.RIGHT)
       
    #Add so Color Max can be updated
            self.bcmax = Tk.StringVar()
            self.bcmax.set('{0:5.2f}'.format(0))
            self.bcmaxval = Tk.Entry(self,textvariable=self.bcmax,width=10)
            self.bcmaxval.bind("<Return>",self.aia_param)
            self.bcmaxval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Max 
            bcmaxText = Tk.StringVar()
            bcmaxText.set("B Color Max.")
            bcmaxDir = Tk.Label(self,textvariable=bcmaxText,height=4)
            bcmaxDir.pack(side=Tk.RIGHT)

###############################################
#     GREEN COLOR BOXES                       #
###############################################
    #Add so Color Min can be updated
            self.gcmin = Tk.StringVar()
            self.gcmin.set('{0:5.2f}'.format(0))
            self.gcminval = Tk.Entry(self,textvariable=self.gcmin,width=10)
            self.gcminval.bind("<Return>",self.aia_param)
            self.gcminval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Min
            gcminText = Tk.StringVar()
            gcminText.set("G Color Min.")
            gcminDir = Tk.Label(self,textvariable=gcminText,height=4)
            gcminDir.pack(side=Tk.RIGHT)
       
    #Add so Color Max can be updated
            self.gcmax = Tk.StringVar()
            self.gcmax.set('{0:5.2f}'.format(0))
            self.gcmaxval = Tk.Entry(self,textvariable=self.gcmax,width=10)
            self.gcmaxval.bind("<Return>",self.aia_param)
            self.gcmaxval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Max 
            gcmaxText = Tk.StringVar()
            gcmaxText.set("G Color Max.")
            gcmaxDir = Tk.Label(self,textvariable=gcmaxText,height=4)
            gcmaxDir.pack(side=Tk.RIGHT)
###############################################
#     RED COLOR BOXES                         #
###############################################

    #Add so Color Min can be updated
            self.rcmin = Tk.StringVar()
            self.rcmin.set('{0:5.2f}'.format(0))
            self.rcminval = Tk.Entry(self,textvariable=self.rcmin,width=10)
            self.rcminval.bind("<Return>",self.aia_param)
            self.rcminval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Min
            rcminText = Tk.StringVar()
            rcminText.set("R Color Min.")
            rcminDir = Tk.Label(self,textvariable=rcminText,height=4)
            rcminDir.pack(side=Tk.RIGHT)
       
    #Add so Color Max can be updated
            self.rcmax = Tk.StringVar()
            self.rcmax.set('{0:5.2f}'.format(0))
            self.rcmaxval = Tk.Entry(self,textvariable=self.rcmax,width=10)
            self.rcmaxval.bind("<Return>",self.aia_param)
            self.rcmaxval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Max 
            rcmaxText = Tk.StringVar()
            rcmaxText.set("R Color Max.")
            rcmaxDir = Tk.Label(self,textvariable=rcmaxText,height=4)
            rcmaxDir.pack(side=Tk.RIGHT)

    #boxes to create if single wavelength image
        else:
    #Add so Color Min can be updated
            self.cmin = Tk.StringVar()
            self.cmin.set('{0:5.2f}'.format(0))
            self.cminval = Tk.Entry(self,textvariable=self.cmin,width=10)
            self.cminval.bind("<Return>",self.aia_param)
            self.cminval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Min
            cminText = Tk.StringVar()
            cminText.set("Color Min.")
            cminDir = Tk.Label(self,textvariable=cminText,height=4)
            cminDir.pack(side=Tk.RIGHT)
       
    #Add so Color Max can be updated
            self.cmax = Tk.StringVar()
            self.cmax.set('{0:5.2f}'.format(0))
            self.cmaxval = Tk.Entry(self,textvariable=self.cmax,width=10)
            self.cmaxval.bind("<Return>",self.aia_param)
            self.cmaxval.pack(side=Tk.RIGHT,padx=5,pady=5)
    #set up Color Max 
            cmaxText = Tk.StringVar()
            cmaxText.set("Color Max.")
            cmaxDir = Tk.Label(self,textvariable=cmaxText,height=4)
            cmaxDir.pack(side=Tk.RIGHT)
       
#set up Submenu
        menubar = Tk.Menu(self.parent)
        self.parent.config(menu=menubar)

        fileMenu = Tk.Menu(menubar)
        subMenu = Tk.Menu(fileMenu)
#create another item in menu
        fileMenu.add_separator()

        fileMenu.add_command(label='Exit',underline=0,command=self.onExit)


    #set AIA parameters
    def aia_param(self,event):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        try:
            self.h0 = float(self.h0val.get())
            self.w0 = float(self.w0val.get())
            self.cx = float(self.cxval.get())
            self.cy = float(self.cyval.get())

            #update the color parameters
            #color table update if 3 color image
            if self.color3:

                #Update R color scale
                self.rmin = float(self.rcminval.get())
                self.rmax = float(self.rcmaxval.get())
                self.img_scale[self.wav[0]][1] = self.rmin
                self.img_scale[self.wav[0]][2] = self.rmax

                #Update G color scale
                self.gmin = float(self.gcminval.get())
                self.gmax = float(self.gcmaxval.get())
                self.img_scale[self.wav[1]][1] = self.gmin
                self.img_scale[self.wav[1]][2] = self.gmax

                #Update B color scale
                self.bmin = float(self.bcminval.get())
                self.bmax = float(self.bcmaxval.get())
                self.img_scale[self.wav[2]][1] = self.bmin
                self.img_scale[self.wav[2]][2] = self.bmax
        
                #recreate 3 color image
                self.create_3color()

            #color table if single image
            else:
                self.ivmin = float(self.cminval.get())
                self.ivmax = float(self.cmaxval.get())
                self.img_scale[self.wav][1] = self.ivmin
                self.img_scale[self.wav][2] = self.ivmax
        #now replot
            self.clicked = True
            self.sub_window()
            self.aia_plot()
            self.clicked = False

#error if not floats
        except ValueError:
            self.error = 10
            self.onError()



#Exits the program
    def onExit(self):
       plt.clf()
       plt.close()
       self.quit()
       self.parent.destroy()


#Command to increase the order to plot new aia image
    def increaseorder(self):
        self.order = self.order+1
        if self.order > self.ordermax:
            self.order = 0 
        self.sorder.set(str(int(self.order)))
        self.clicked = True
        self.a.clear()
        self.aia_set()
        self.aia_plot()
        self.clicked = False
#Command to decrease order to plot new aia image
    def decreaseorder(self):
        self.order = self.order-1
        if self.order < 0:
            self.order = self.ordermax
        self.sorder.set(str(int(self.order)))
        self.clicked = True
        self.a.clear()
        self.aia_set()
        self.aia_plot()
        self.clicked = False


    #create 3 color image
    def create_3color(self):
        self.wav = []
        for j,i in enumerate(self.img):
            self.wav.append('{0:4.0f}'.format(i.wavelength.value).replace(' ','0'))
            #set normalized scaling for every observation
            self.ivmin = self.img_scale[self.wav[j]][1]
            self.ivmax = self.img_scale[self.wav[j]][2]
            prelim = (np.arcsinh(i.data/i.exposure_time.value)-self.ivmin)/self.ivmax

            #replace out of bounds points
            prelim[prelim < 0.] = 0.
            prelim[prelim > 1.] = 1.
            self.img3d[:,:,j] = prelim

            #set the string value in the plot window
            if j == 0:
                self.rcmin.set('{0:9.3}'.format(self.ivmin))
                self.rcmax.set('{0:9.3}'.format(self.ivmax))
            if j == 1:
                self.gcmin.set('{0:9.3}'.format(self.ivmin))
                self.gcmax.set('{0:9.3}'.format(self.ivmax))
            if j == 2:
                self.bcmin.set('{0:9.3}'.format(self.ivmin))
                self.bcmax.set('{0:9.3}'.format(self.ivmax))

            #Retrieved and set the time value based on R image
            self.obs_time = self.img[0].date

# get the image properties
    def img_prop(self):

        #different plotting properties if color3 set
        if self.color3:
            self.create_3color()
        else:         
            self.wav ='{0:4.0f}'.format(self.img.wavelength.value).replace(' ','0')
            #use default color tables
            self.icmap = self.img_scale[self.wav][0]
            self.ivmin = self.img_scale[self.wav][1]
            self.ivmax = self.img_scale[self.wav][2]
 
            #set the string value in the plot window
            self.cmin.set('{0:9.3}'.format(self.ivmin))
            self.cmax.set('{0:9.3}'.format(self.ivmax))
 
            #Retrieved and set the time value
            self.obs_time = self.img.date
 
    def text_loc(self):
#set text location
        #if self.w0 > self.h0:
        #    self.txtx = -(self.w0-self.h0)
        #    self.txty = (self.maxy-self.miny)*0.01
        #elif self.w0 < self.h0:
        #    self.txty = -(self.h0-self.w0)
        #    self.txtx = (self.maxx-self.minx)*0.01
        #if self.w0 == self.h0:
        self.txtx = (self.maxx-self.minx)*0.01
        self.txty = (self.maxy-self.miny)*0.01



#plot the current AIA image
    def aia_plot(self):
       #clear the current image
       self.a.clear()
       #find where to put the plotting information
       self.text_loc()


       #Make 3 color plot
       if self.color3:
           self.a.imshow(self.img3d,interpolation='none',origin='lower',extent=[self.minx,self.maxx,self.miny,self.maxy])
           #set the observation time which will be use for rotation
           self.a.text(self.minx+self.txtx,self.miny+self.txty,'AIA {0}/{1}/{2}'.format(*self.wav)+'- {0}Z'.format(self.obs_time.strftime('%Y/%m/%d - %H:%M:%S')),color='white',fontsize=14,zorder=50,fontweight='bold')
       else:
           self.a.imshow(self.data0,interpolation='none',cmap=self.icmap,origin='lower',vmin=self.ivmin,vmax=self.ivmax,extent=[self.minx,self.maxx,self.miny,self.maxy])
           #show current date
           self.a.text(self.minx+self.txtx,self.miny+self.txty,'{0}Z'.format(self.obs_time.strftime('%Y/%m/%d - %H:%M:%S')),color='white',fontsize=14,zorder=50,fontweight='bold')

       self.a.set_xlabel('Arcseconds')
       self.a.set_ylabel('Arcseconds')
       if ((self.clicked) & (self.firstclick)):
           #make sure axis is on if not turn it on
           if not self.x.get_frame_on(): self.x.set_axis_on()
           #Show the clicked region in a separate plot
           self.x.clear()
           #3 color image
           if self.color3:
               self.x.imshow(self.img3d,interpolation='none',origin='lower',extent=[self.minx,self.maxx,self.miny,self.maxy])
           else:
               self.x.imshow(self.data0,interpolation='none',cmap=self.icmap,origin='lower',vmin=self.ivmin,vmax=self.ivmax,extent=[self.minx,self.maxx,self.miny,self.maxy])

           self.x.set_xlim([min(self.xbox),max(self.xbox)])           
           self.x.set_ylim([min(self.ybox),max(self.ybox)])           
           self.x.scatter(self.cx,self.cy,marker='x',color='red',s=35,zorder=499)
           self.x.set_xlabel('Arcseconds')
           self.x.set_ylabel('Arcseconds')

           #show the selected region on the big plot
           self.a.scatter(self.cx,self.cy,marker='x',color='red',s=35,zorder=499)
           self.a.plot(self.xbox,self.ybox,color='black',linewidth=5,zorder=500)
           self.a.plot(self.xbox,self.ybox,'--',color='white',linewidth=3,zorder=501)
       self.canvas.draw()
   
#set variables spectrum of a given order
    def aia_set(self):

       #set current index depending on 3color image
       if self.color3:
           self.infile = self.flist[self.order]
       #else single file
       else:
           self.infile = self.flist[self.order]

       #put file into sunpy map
       self.img = sunpy.map.Map(self.infile)
       self.maxx,self.minx,self.maxy,self.miny = self.img_extent() #get extent of image for coverting pixel into physical
       #3 color image
       if self.color3:
           self.img3d = np.zeros((self.img[0].data.shape[0],self.img[0].data.shape[1],3))
           self.scale = [self.img[0].scale[0].value,self.img[0].scale[1].value] # get x, y image scale
       #single color image
       else:
           self.data0 = np.arcsinh(self.img.data/self.img.exposure_time.value) #reference the data plot seperately
           self.scale = [self.img.scale[0].value,self.img.scale[1].value] # get x, y image scale
       #set aia plotting preferences
       self.img_prop()

    def img_extent(self):
    #get only physical values for first image if color 3
        if self.color3:
        # get the image coordinates in pixels
            px0 = self.img[0].meta['crpix1']
            py0 = self.img[0].meta['crpix2']
        # get the image coordinates in arcsec 
            ax0 = self.img[0].meta['crval1']
            ay0 = self.img[0].meta['crval2']
        # get the image scale in arcsec 
            axd = self.img[0].meta['cdelt1']
            ayd = self.img[0].meta['cdelt2']
        #get the number of pixels
            tx,ty = self.img[0].data.shape
        else:
        # get the image coordinates in pixels
            px0 = self.img.meta['crpix1']
            py0 = self.img.meta['crpix2']
        # get the image coordinates in arcsec 
            ax0 = self.img.meta['crval1']
            ay0 = self.img.meta['crval2']
        # get the image scale in arcsec 
            axd = self.img.meta['cdelt1']
            ayd = self.img.meta['cdelt2']
        #get the number of pixels
            tx,ty = self.img.data.shape
        #get the max and min x and y values
        pminx,pmaxx = 0.-px0,tx-px0
        pminy,pmaxy = 0.-py0,ty-py0
    #convert to arcsec
        maxx,minx = ax0+pmaxx*axd,ax0+pminx*axd
        maxy,miny = ay0+pmaxy*ayd,ay0+pminy*ayd

        return maxx,minx,maxy,miny
#Basic click event
    def on_click_event(self,click):
#Click envents for continuum selection
#Make sure you click inside the plot
        try:
            #test to make sure the data are in the plot
            test = click.xdata-0.
            test = click.ydata-0.

            #store the physical value of clicked points
            self.cx = click.xdata
            self.cy = click.ydata


            #tell if the plot has beeen clicked at least once
            self.firstclick = True
            #tell the plot its been clicked
            self.clicked = True

            #create subwindow of selected region
            self.sub_window()

            #update the x and y parameters in bottom box
            self.scx.set('{0:5.2f}'.format(self.cx))
            self.scy.set('{0:5.2f}'.format(self.cy))

            #plot new cutout box
            self.aia_plot()


            #update the plot parameters

            #reset to no click
            self.clicked = False
            
#Throw error if clicked outside the plot
        except TypeError:
            self.error = 20
            self.onError()

   #create window for plotting
    def sub_window(self):
        self.xbox = [self.cx-(self.scale[0]*self.w0/2.),self.cx-(self.scale[0]*self.w0/2.),self.cx+(self.scale[0]*self.w0/2.),self.cx+(self.scale[0]*self.w0/2.),self.cx-(self.scale[0]*self.w0/2.)]
        self.ybox = [self.cy-(self.scale[1]*self.h0/2.),self.cy+(self.scale[1]*self.h0/2.),self.cy+(self.scale[1]*self.h0/2.),self.cy-(self.scale[1]*self.h0/2.),self.cy-(self.scale[1]*self.h0/2.)]

#Function for retrieving order from popup
    def on_order(self):
        m = 0
        while m == 0:
            try:
                inputO = order_popup(root)
                root.wait_window(inputO.top)
                order = int(inputO.order)
                if ((order > 0) & (order <= self.ordermax)):
                    m = 1
                else:
#Error order is out of range
                    self.error = 3
                    self.onError()
            except ValueError:
#Error order is not an integer
                self.error = 4
                self.onError()

        return order

#Function for retrieving order from text box
    def on_order_box(self,event):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        m = 0
        while m == 0:
            try:
                order = self.orderval.get()
                order = int(order)
                if ((order > 0) & (order <= self.ordermax)):
                    m = 1
                    self.order = order
                    self.aia_set()
                    self.aia_plot()
                else:
#Error order is out of range
                    self.error = 3
                    self.onError()
            except ValueError:
#Error order is not an integer
                self.error = 4
                self.onError()
            




#Tells Why Order information is incorrect
    def onError(self):
        if self.error == 1:
            box.showerror("Error","File Not Found")
        if self.error == 4:
            box.showerror("Error","Value Must be an Integer")
        if self.error == 6:
            box.showerror("Error","File is not in Fits Format")
        if self.error == 10:
            box.showerror("Error","Value Must be Float")
        if self.error == 20:
            box.showerror("Error","Must Select Inside Plot Bounds")
               
#main loop
def main():
    global root
    root = Tk.Tk()
    app = gui_c(root)
    root.mainloop()


if __name__=="__main__":
#create root frame
   main()

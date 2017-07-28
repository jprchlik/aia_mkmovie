import matplotlib
matplotlib.use('TkAgg')


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#implement the deault mpl key bindings
from matplotlib.backend_bases import key_press_handler,MouseEvent
import tkMessageBox as box
import tkFileDialog as Tkf
import numpy as np
import sys
import matplotlib.pyplot as plt

try:                                                                                      
    import sunpy.map                                                                      
    from sunpy.cm import cm                                                               
except ImportError:                                                                       
    sys.stdout.write(sys.stderr,"sunpy not installed, use pip install sunpy --upgrade")   

#check the python version to use one Tkinter syntax or another
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk


#main gui class
class gui_c(Tk.Frame):

    #init gui class
    def __init__(self,parent,flist,sday=False,eday=False,w0=1900.,h0=1200.,cy=0,cx=0):
        Tk.Frame.__init__(self,parent,background='white') #create initial frame with white background

        #file list
        if isinstance(flist,str):
            self.flist = [flist]
        elif isinstance(flist,list):
            self.flist = flist
        else:
            sys.stdout.write('flist must be a string or list')

        #set the starting list to be 0
        self.order = 0
        #get maximum value in list
        self.ordermax = len(flist)-1

        self.parent = parent
        self.w0 = w0
        self.h0 = h0
        self.cy = cy
        self.cx = cx

        #clicked point on the figure
        self.clicked = False

        #Dictionary for vmax, vmin, and color
        self.img_scale = {'0094':[cm.sdoaia94  ,np.arcsinh(1.),np.arcsinh(150.)],
                          '0131':[cm.sdoaia131 ,np.arcsinh(1.),np.arcsinh(500.)],
                          '0171':[cm.sdoaia171 ,np.arcsinh(10.),np.arcsinh(2500.)],
                          '0193':[cm.sdoaia193 ,np.arcsinh(10.),np.arcsinh(4500.)],
                          '0211':[cm.sdoaia211 ,np.arcsinh(10.),np.arcsinh(4000.)],
                          '0304':[cm.sdoaia304 ,np.arcsinh(2.),np.arcsinh(300.)],
                          '0335':[cm.sdoaia335 ,np.arcsinh(1.),np.arcsinh(100.)],
                          '1600':[cm.sdoaia1600,np.arcsinh(20.),np.arcsinh(500.)],
                          '1700':[cm.sdoaia1700,np.arcsinh(200.),np.arcsinh(4000.)]}



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

#set up center x box
        w0Text = Tk.StringVar()
        w0Text.set("Width (pixels)")
        w0Dir = Tk.Label(self,textvariable=w0Text,height=4)
        w0Dir.pack(side=Tk.LEFT)
#Add so center x can be updated
        w0 = Tk.StringVar()
        w0.set('{0:5.2f}'.format(self.w0))
        self.w0val = Tk.Entry(self,textvariable=w0,width=10)
        self.w0val.bind("<Return>",self.aia_param)
        self.w0val.pack(side=Tk.LEFT,padx=5,pady=5)
       
#set up center w0 box
        h0Text = Tk.StringVar()
        h0Text.set("Height (pixels)")
        h0Dir = Tk.Label(self,textvariable=h0Text,height=4)
        h0Dir.pack(side=Tk.LEFT)
#Add so center h0 can be updated
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
        self.orderval = Tk.Entry(self,textvariable=self.sorder,width=10)
        self.orderval.bind("<Return>",self.on_order_box)
        self.orderval.pack(side=Tk.LEFT,padx=5,pady=5)
       
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


# get the image properties
    def img_prop(self):
       self.wav ='{0:4.0f}'.format(self.img.wavelength.value).replace(' ','0')
       #use default color tables
       self.icmap = self.img_scale[self.wav][0]
       self.ivmin = self.img_scale[self.wav][1]
       self.ivmax = self.img_scale[self.wav][2]


#plot the current AIA image
    def aia_plot(self):
       self.a.clear()
       self.a.imshow(self.data0,interpolation='none',cmap=self.icmap,origin='lower',vmin=self.ivmin,vmax=self.ivmax,extent=[self.minx,self.maxx,self.miny,self.maxy])
       self.a.set_xlabel('Arcseconds')
       self.a.set_ylabel('Arcseconds')
       if self.clicked:
           #Show the clicked region in a separate plot
           self.x.clear()
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

       #set current index
       self.infile = self.flist[self.order]
       #put file into sunpy map
       self.img = sunpy.map.Map(self.infile)
       self.maxx,self.minx,self.maxy,self.miny = self.img_extent() #get extent of image for coverting pixel into physical
       self.data0 = np.arcsinh(self.img.data) #reference the data plot seperately
       self.scale = [1./self.img.scale[0].value,1./self.img.scale[1].value] # get x, y image scale
       #set aia plotting preferences
       self.img_prop()

    def img_extent(self):
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
        minx,maxx = px0-tx,tx-px0
        miny,maxy = py0-ty,ty-py0
    #convert to arcsec
        maxx,minx = maxx*axd,minx*axd
        maxy,miny = maxy*ayd,miny*ayd

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
                    self.spec_set()
                    self.spec_plot()
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

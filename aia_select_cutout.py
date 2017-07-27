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
        self.f,self.a = plt.subplots(figsize=(8*aratio,8*aratio*.5))
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
        w0Text.set("Center X Position")
        w0Dir = Tk.Label(self,textvariable=w0Text,height=4)
        w0Dir.pack(side=Tk.LEFT)
#Add so center x can be updated
        w0 = Tk.StringVar()
        w0.set('{0:5.2}'.format(self.w0))
        self.vxval = Tk.Entry(self,textvariable=w0,width=10)
        self.vxval.bind("<Return>",self.stellar_param)
        self.vxval.pack(side=Tk.LEFT,padx=5,pady=5)
       
#set up center x box
        h0Text = Tk.StringVar()
        h0Text.set("Center X Position")
        h0Dir = Tk.Label(self,textvariable=h0Text,height=4)
        h0Dir.pack(side=Tk.LEFT)
#Add so center x can be updated
        h0 = Tk.StringVar()
        h0.set('{0:5.2}'.format(self.h0))
        self.vxval = Tk.Entry(self,textvariable=h0,width=10)
        self.vxval.bind("<Return>",self.stellar_param)
        self.vxval.pack(side=Tk.LEFT,padx=5,pady=5)
       
   
#set up center x box
        cxText = Tk.StringVar()
        cxText.set("Center X Position")
        cxDir = Tk.Label(self,textvariable=cxText,height=4)
        cxDir.pack(side=Tk.LEFT)
#Add so center x can be updated
        cx = Tk.StringVar()
        cx.set('{0:5.2}'.format(fself.cx))
        self.cxval = Tk.Entry(self,textvariable=cx,width=10)
        self.cxval.bind("<Return>",self.stellar_param)
        self.cxval.pack(side=Tk.LEFT,padx=5,pady=5)
       
#set up center y boy
        cyTeyt = Tk.StringVar()
        cyTeyt.set("Center Y Position")
        cyDir = Tk.Label(self,teytvariable=cyText,height=4)
        cyDir.pack(side=Tk.LEFT)
#Add so center y can be updated
        cy = Tk.StringVar()
        cy.set('{0:5.2}'.format(fself.cy))
        self.cyval = Tk.Entry(self,textvariable=cy,width=10)
        self.cyval.bind("<Return>",self.stellar_param)
        self.cyval.pack(side=Tk.LEFT,padx=5,pady=5)
       
#set up Submenu
        menubar = Tk.Menu(self.parent)
        self.parent.config(menu=menubar)

        fileMenu = Tk.Menu(menubar)
        subMenu = Tk.Menu(fileMenu)
#create another item in menu
        fileMenu.add_separator()

        fileMenu.add_command(label='Exit',underline=0,command=self.onExit)

#Exits the program
    def onExit(self):
       self.quit()


#Command to increase the order to plot new aia image
    def increaseorder(self):
        self.order = self.order+1
        if self.order > self.ordermax:
            self.order = 1 
        self.sorder.set(str(int(self.order)))
        self.a.clear()
        self.aia_set()
        self.aia_plot()

#Command to decrease order to plot new aia image
    def decreaseorder(self):
        self.order = self.order-1
        if self.order < 1:
            self.order = self.ordermax
        self.sorder.set(str(int(self.order)))
        self.a.clear()
        self.aia_set()
        self.aia_plot()


# get the image properties
    def img_prop(self):
       self.wav ='{0:4.0f}'.format( img.wavelength.value).replace(' ','0')
       #use default color tables
       sel.ficmap = self.img_scale[self.wav][0]
       sel.fivmin = self.img_scale[self.wav][1]
       sel.fivmax = self.img_scale[self.wav][2]


#plot the current AIA image
    def aia_plot(self):
       self.a.clear()
       self.a.imshow(self.data0,interpolation='none',cmap=self.icmap,origin='lower',vmin=self.ivmin,vmax=self.ivmax,extent=[minx,maxx,miny,maxy])
       self.a.scatter(click.xdata,click.ydata,marker='x',color='red',s=35,zorder=499)
       self.a.plot(self.xbox,self.ybox,color='black',linewidth=4,zorder=500)
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
            test = click.ydata-0.
            test = click.xdata-0.
            #self.contx.append(click.xdata)
            #self.conty.append(click.ydata)
            #self.pcontx.append(click.x)
            #self.pconty.append(click.y)
            self.cx = click.xdata
            self.cy = click.ydata


            self.xbox = [click.xdata-(self.scale*self.w0/2.),click.xdata-(self.scale*self.w0/2.),click.xdata+(self.scale*self.w0/2.),click.xdata+(self.scale*self.w0/2.),click.xdata-(self.scale*self.w0/2.)]
            self.ybox = [click.ydata-(self.scale*self.h0/2.),click.ydata-(self.scale*self.h0/2.),click.ydata+(self.scale*self.h0/2.),click.ydata+(self.scale*self.h0/2.),click.ydata-(self.scale*self.h0/2.)]



            
#Throw error if clicked outside the plot
        except TypeError:
            self.error = 20
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

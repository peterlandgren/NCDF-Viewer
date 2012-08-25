#!/usr/bin/env python
#Runs just one window, structure changed

from scipy import *
#from pylab import *
from numpy import *
import numpy as N

import sys, os

from numpy import linspace, meshgrid, pi, sin, array

# Enthought library imports
from enable.api import Component, ComponentEditor, BaseTool, Window, LineStyle
from enable.component_editor import ComponentEditor

from traits.api import HasTraits, Instance, DelegatesTo, on_trait_change, Trait, Event, Array, Any, Callable, Int, Tuple, Enum, Range, Float, Bool, Button, File, Str
from traitsui.api import Item, Group, View, HGroup, VGroup, Handler, RangeEditor
from traits.util.resource import find_resource
from traitsui.file_dialog import open_file
#from traitsui.file_dialog import open_file


# Chaco imports
from chaco.api import ArrayPlotData, Plot, Blues, HPlotContainer, VPlotContainer, create_line_plot, GridDataSource, ImageData, Plot, PlotAxis, LinearMapper, GridMapper, DataRange2D, DataRange1D, CMapImagePlot, ArrayDataSource, DataView, ContourPolyPlot, CMapImagePlot, ColorBar, PlotLabel, PlotAxis, ImagePlot, OverlayPlotContainer, add_default_axes, LabelAxis, PlotGraphicsContext
from chaco.tools.api import PanTool, ZoomTool, LineInspector
from scipy.interpolate import interp1d
from scipy.integrate.quadrature import cumtrapz
from chaco.default_colormaps import *
from chaco.default_colormaps import gist_gray as jet
from chaco.overlays import *

# RH code import functions
from myrh15d import *

#Interpolate functions
from scipy import interpolate

#My functions
from ImageloadSourceFunctions import *

"""
Master program to display NCDF files and source function plots.

Can be compiled as a py2app binary program with modifications to FileDialogDemo as specified in comments.

Make sure to clean columns in file before use.

Can easily be modified to accept different file types for the main image viewer.
    To accept alterate types edit where directed in MainPlot.
    It should only be necessary to edit one spot in this file. Data should be entered into generic objects as directed.

Program written by Peter Landgren under direction of Tiago Pereira.  Contact Peter at plandgren13@my.whitworth.edu or 971-404-4320 (U.S.)
"""

class FileDialogDemo(HasTraits):
    """ Class object used to importing file names and changing certain settings up front.
    Must be modified to compile and use as a binary in py2app.
    
    """

    #The name of the selected file:
    file_name_Ray = File('/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom_IRIS_nolines/output_ray.ncdf')
    file_name_Indata = File('/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom_IRIS_nolines/output_indata.ncdf')
    #When compiling the py2app binary file use string versions below, as py2app fails when attempting to include file icon image.  Also comment out from 'traitsui.file_dialog import open_file' earlier
    #file_name_Ray = '/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom_IRIS_nolines/output_ray.ncdf'
    #file_name_Indata = '/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom_IRIS_nolines/output_indata.ncdf'
    welcomestring = Str("Welcome to the NCDF File Viewer.  Created by Peter Landgren in August 2012 at Lockheed Martin Space Systems. \n Please contact Peter Landgren at plandgren13@my.whitworth.edu with questions.")
    Wave_Ref = Float()
    Wave_Ref = 279.55
    VMax_Min = Float()
    VMax_Min = 59
    ZRangeMax = Float()
    ZRangeMin = Float()
    ZRangeMax = 3.
    ZRangeMin = 0.2




    #-- Traits View Definitions ------------------------------------------------

    view = View(VGroup(
	    Item('welcomestring', style = 'readonly', show_label = False),
	
            Item( 'file_name_Ray', springy = True ),
     
            Item( 'file_name_Indata', springy = True ),
            
            HGroup(Item('Wave_Ref', springy = True), Item('VMax_Min', springy = True), Item('ZRangeMax', springy = True), Item('ZRangeMin', springy = True)),
	    
	
        ),
	buttons = ['OK'],
	title = "Please Select Files",
        width = 1000
    )

#Created and parses Filenames and initial parameters
FileLoadWindow = FileDialogDemo()
FileLoadWindow.configure_traits()
FileLoadWindow.file_name_Ray = FileLoadWindow.file_name_Ray.strip()
FileLoadWindow.file_name_Indata = FileLoadWindow.file_name_Indata.strip()


class MainPlot(HasTraits):
    """Main class for the primary solar image and point-by-point spectra.
    
    Class creates image inspector tool object.
    
    Controls file saving.
    
    Main window that passes information on to source function plots
    """
    
    #Place to extend program for extra filetypes
    Name, fileExtension = os.path.splitext(FileLoadWindow.file_name_Ray)
    if fileExtension == '.ncdf':
	#Open .ncdf File to read
	RayFile = Rh15dout()
	RayFile.read_ray(FileLoadWindow.file_name_Ray)
	WavelengthData = RayFile.ray.wavelength
	IntensityData = RayFile.ray.intensity
    elif fileExtension == 'Alternate file type':
	pass
	#Add code here to read in new file type
    else:
	print "Ray File Extension Not Recognized"
    

    wave_ref = float(FileLoadWindow.Wave_Ref)
    intensityindex =  N.argmin(N.abs(WavelengthData[:]-wave_ref))
    #self.RayFile.read_ray('/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom/output_ray.ncdf')
    WavelengthIndex = Range(value = int(intensityindex), low = 0, high = WavelengthData.shape[0] - 1)
    WavelengthFine = Range(value = wave_ref, low = wave_ref - 0.1, high = wave_ref + 0.1)
    WavelengthView = Float(wave_ref)
    Wavelength = Range(value = wave_ref, low = float(N.min(WavelengthData[:])), high = float(N.max(WavelengthData[:])))
    save = Button('save')
    SaveFlag = Bool()
    VelocityView = Float()
    ViewPosition = Array()


    
    
    colormapminabs = float(N.min(IntensityData[:,:,intensityindex]))
    colormapmaxabs = float(N.max(IntensityData[:,:,intensityindex]))
    colormapmin = Range(low=colormapminabs, high=colormapmaxabs, value=colormapminabs)
    colormapmax = Range(low=colormapminabs, high=colormapmaxabs, value=colormapmaxabs)

    
    PrimaryPlotC = Instance(HPlotContainer)
    windowwidth = 1500
    ImageLock = Bool()
    colormap_list = Enum('jet', 'gist_gray', 'gist_heat', 'Blues', 'hot')
    Spectraplotzoom = Enum('Close', 'Far')
    Markercolor = Enum('white', 'green', 'yellow', 'red')

    
    traits_view = View(VGroup(Item('PrimaryPlotC', editor = ComponentEditor(), show_label = False, springy = True)),
		       Item('colormapmin', label = 'Minimum Intensity'), Item('colormapmax', label = 'Maximum Intensity'),
		       HGroup(Item('WavelengthIndex', label = 'Index', style = 'simple', editor = RangeEditor(mode = 'slider', low=0, high=WavelengthData.shape[0] - 1)), Item('WavelengthView', label = 'Wavelength')),
		       HGroup(Item('colormap_list', style = 'simple', label = 'Colormap'), Item('Spectraplotzoom', style = 'custom', label = 'Spectral Zoom Range'), Item('Markercolor', style = 'simple', label = 'Marker Color'), Item('VelocityView', label = 'Velocity (km/s)', style = 'readonly'), Item('ViewPosition', label = 'View Location', style = 'readonly'), Item('save', show_label = False)),
		       width = windowwidth, height = windowwidth/2.0 * 1.3, resizable = True, title = "Main Image Plot")

    #Array used for passing information to speactra line plto and also Source Function Plots.  Will be synced
    InspectorPosition = Array()

    def __init__(self):
	super(MainPlot, self).__init__()
	self.ImageLock = False
	self.SaveFlag = False
	
	self.create_PrimaryPlotC()
	
    def create_PrimaryPlotC(self):
	#Extracts the data for the main plot
        self.MainPlotData = self.IntensityData[:,:,self.intensityindex]
	self.markerplotdatax = []
	self.markerplotdatay = []
	self.SpectraMultiple = 1.0e+8
	WavelengthXPoints = [self.Wavelength, self.Wavelength]
	WavelengthYPoints = [N.min(self.IntensityData[:,:])*self.SpectraMultiple, N.max(self.IntensityData[:,:])*self.SpectraMultiple]
	WavelengthXPointsStatic = WavelengthXPoints
	WavelengthYPointsStatic = WavelengthYPoints

	AvgIntensity = N.mean(N.mean(self.IntensityData, axis = 0), axis = 0)
	print "mean computed-------------------"

        #Create main Plot (intensity plot)
        self.Mainplotdata = ArrayPlotData(Mainimagedata = self.MainPlotData, markerplotdatax = self.markerplotdatax, markerplotdatay = self.markerplotdatay)
        self.Mainplot = Plot(self.Mainplotdata)
        self.Main_img_plot = self.Mainplot.img_plot("Mainimagedata", colormap = jet)[0]
	
	#Create marker wneh x is pressed
	self.Mainplot.MarkerPlot = self.Mainplot.plot(("markerplotdatax","markerplotdatay"),type = 'line', color = 'white')

        #Create overlaid crosshairs for Main plot
        LineInspector1 = LineInspector(component = self.Main_img_plot,axis = 'index_x', write_metadata=True,is_listener = False,inspect_mode="indexed")            
        self.Main_img_plot.overlays.append(LineInspector1)
        LineInspector2 = LineInspector(component = self.Main_img_plot,axis = 'index_y', write_metadata=True,is_listener = False,inspect_mode="indexed")
        self.Main_img_plot.overlays.append(LineInspector2)

        #Create overlay tools and add them to main plot
        Main_imgtool = ImageInspectorTool(self.Main_img_plot)
        self.Main_img_plot.tools.append(Main_imgtool)
        Main_overlay = ImageInspectorOverlay(component = self.Main_img_plot, image_inspector = Main_imgtool, bgcolor = "white", border_visible = True)
        self.Main_img_plot.overlays.append(Main_overlay)

        #Sync up inspector position so it can be passed to the spectra plot
        Main_overlay.sync_trait('InspectorPosition', self, 'InspectorPosition', mutual = True)
	Main_imgtool.sync_trait('ImageLock', self, 'ImageLock')
	Main_imgtool.sync_trait('SaveFlag', self, 'SaveFlag')
	
	#Sync up max and min colormap value
	self.Main_img_plot.value_range.sync_trait('low', self, 'colormapmin', mutual = True)
	self.Main_img_plot.value_range.sync_trait('high', self, 'colormapmax', mutual = True)
	
        #Create spectra plot for a single column
	self.Spectraplotdata = ArrayPlotData(x = self.WavelengthData[:], y = (self.IntensityData[2,2] * self.SpectraMultiple), WavelengthXPointsStatic = WavelengthXPointsStatic, WavelengthYPointsStatic = WavelengthYPointsStatic, WavelengthXPoints = WavelengthXPoints, WavelengthYPoints = WavelengthYPoints, AvgIntensity = AvgIntensity * self.SpectraMultiple, is_listener = True)
        self.Spectraplot1 = Plot(self.Spectraplotdata)
        self.Spectraplot1.plot(("x","y"), type = "line", color = "blue")
        self.Spectraplot1.plot( ("WavelengthXPointsStatic","WavelengthYPointsStatic"), type = 'line', line_style = 'dash', color = 'green')
	self.Spectraplot1.plot( ("WavelengthXPoints","WavelengthYPoints"), type = 'line', line_style = 'dash', color = 'red')
	self.Spectraplot1.plot( ("x","AvgIntensity"), type = 'line', color = 'black')
        
	#Change Plot characteristics
	#Sets width around waveref to examine
	self.xlowspread = 1.
        self.xhighspread = 1.
        self.Spectraplot1.range2d.x_range.set_bounds(self.wave_ref - self.xlowspread, self.wave_ref + self.xhighspread)
        self.Spectraplot1.x_axis.title = "Wavelength (A)"
        self.Spectraplot1.y_axis.title = "Intensity (J s^-1 m^-2 Hz^-1 sr^-1)"
        self.rangearray = self.IntensityData[:,:]
        self.rangearray = self.rangearray[:,:,N.argmin(N.abs(self.WavelengthData[:]-(self.wave_ref-self.xlowspread))):N.argmin(N.abs(self.WavelengthData[:]-(self.wave_ref+self.xhighspread)))]
        self.Spectraplot1.range2d.y_range.set_bounds(N.min(self.rangearray[:,:])*self.SpectraMultiple, N.max(self.rangearray[:,:])   * self.SpectraMultiple)
	#self.Spectraplot1.range2d.y_range.set_bounds(N.min(self.rangearray[:,:])*self.SpectraMultiple, Scaling * (N.max(self.rangearray[:,:]) - N.min(self.rangearray[:,:])) + N.min(self.rangearray[:,:])  * self.SpectraMultiple)


        #add some standard tools. Note, I'm assigning the PanTool to the 
        #right mouse-button to avoid conflicting with the cursors
        self.Mainplot.tools.append(PanTool(self.Main_img_plot, drag_button="right"))
	self.Spectraplot1.tools.append(PanTool(self.Spectraplot1, drag_button="right"))
        self.Mainplot.overlays.append(ZoomTool(self.Main_img_plot))
        self.Spectraplot1.overlays.append(ZoomTool(self.Spectraplot1))

        #Changing interactive options
        zoom = ZoomTool(component=self.Mainplot, tool_mode="box", always_on=False)

        #Create Container for main plot                                  
        MainContainer = HPlotContainer(self.Mainplot, self.Spectraplot1, background = "lightgray", use_back_buffer = True)
        MainContainer.spacing = 25
        self.PrimaryPlotC = MainContainer
	

	
    def _InspectorPosition_changed(self):
	if self.ImageLock == True:
	    return
	
	if len(self.InspectorPosition) > 0:
	    xi = self.InspectorPosition[0]
	    yi = self.InspectorPosition[1]
	else:
	    xi = 0
	    yi = 0
	    return
	    
	self.ViewPosition = self.InspectorPosition    
	self.Spectraplotdata.set_data("y", self.IntensityData[self.InspectorPosition[0],self.InspectorPosition[1]] * self.SpectraMultiple)
	
    def _WavelengthIndex_changed(self):
	self.intensityindex = self.WavelengthIndex
	self.Wavelength = self.WavelengthData[self.intensityindex]
	self.WavelengthView = self.Wavelength
	self.WavelengthXPoints = [self.Wavelength, self.Wavelength]
	self.Spectraplotdata.set_data("WavelengthXPoints", self.WavelengthXPoints)
	self.intensityindex =  N.argmin(N.abs(self.WavelengthData[:]-float(self.Wavelength)))
	self.VelocityView = 299792.45 * (self.Wavelength - self.wave_ref)/self.wave_ref
	self.Mainplotdata.set_data("Mainimagedata", self.IntensityData[:,:,self.intensityindex])
	self.Mainplot.request_redraw()
	
    def _WavelengthView_changed(self):
	self.WavelengthIndex = int(N.argmin(N.abs(self.WavelengthData[:]-self.WavelengthView)))
	
	
    def _WavelengthFine_changed(self):
	self.Wavelength = self.WavelengthFine

    def _colormapmin_changed(self):
	self.Mainplot.request_redraw()
	
    def _colormapmax_changed(self):
	self.Spectraplot1.range2d.y_range.set_bounds(N.min(self.rangearray[:,:])*self.SpectraMultiple, self.colormapmax  * self.SpectraMultiple)
	self.Mainplot.request_redraw()
	
    def _colormap_list_changed(self):
	# get range of current colormap
	clr_range = self.Main_img_plot.color_mapper.range
	color_mapper = eval(self.colormap_list)(clr_range)
	self.Main_img_plot.color_mapper = color_mapper 
	self.Main_img_plot.request_redraw()
	    
    def _Spectraplotzoom_changed(self):
	if (self.Spectraplotzoom == 'Close'):
	    self.Spectraplot1.range2d.x_range.set_bounds(self.wave_ref - self.xlowspread, self.wave_ref + self.xhighspread)
	    self.Spectraplot1.range2d.y_range.set_bounds(N.min(self.rangearray[:,:])*self.SpectraMultiple, self.colormapmax  * self.SpectraMultiple)
	elif (self.Spectraplotzoom == 'Far'):
	    self.Spectraplot1.range2d.x_range.set_bounds(float(N.min(self.WavelengthData[:])), float(N.max(self.WavelengthData[:])))
	    self.Spectraplot1.range2d.y_range.set_bounds(N.min(self.IntensityData[:,:,self.intensityindex]) * self.SpectraMultiple, N.max(self.IntensityData[:,:,self.intensityindex]) * self.SpectraMultiple)
	    
	self.Spectraplot1.request_redraw()
	
    def _ImageLock_changed(self):
	if self.ImageLock == True:
	    xsize = 3
	    xi = self.InspectorPosition[0]
	    yi = self.InspectorPosition[1]
	    self.markerplotdatax = [xi-xsize, xi+xsize, xi, xi-xsize, xi+xsize]
	    self.markerplotdatay = [yi-xsize, yi+xsize, yi, yi+xsize, yi-xsize]
	elif self.ImageLock == False:
	    self.markerplotdatax = []
	    self.markerplotdatay = []
	    
	self.Mainplotdata.set_data("markerplotdatax", self.markerplotdatax)
	self.Mainplotdata.set_data("markerplotdatay", self.markerplotdatay)
	
    def _Markercolor_changed(self):
	self.Mainplot.MarkerPlot[0].color = self.Markercolor
	
    def _save_changed(self):	
    	if self.SaveFlag == True:
	    self.SaveFlag = False
	elif self.SaveFlag == False:
	    self.SaveFlag = True
	
    def _SaveFlag_changed(self):	
	DPI = 72
	print '--------------Save Initiated------------------'
	
	#Main plot save code
	size = (self.IntensityData[:,:,self.intensityindex].shape[0]*4, self.IntensityData[:,:,self.intensityindex].shape[1]*4)
	path = os.getenv('PWD')
	filenamelist = [path, '/', 'MainPlot_WaveLen', str(self.Wavelength).replace('.','_'), '.png']
	filename = ''.join(filenamelist)
	container = self.Main_img_plot
	temp = container.outer_bounds
	container.outer_bounds = list(size)
	container.do_layout(force=True)
	gc = PlotGraphicsContext(size, dpi=DPI)
	gc.render_component(container)
	gc.save(filename)
	container.outer_bounds = temp
	print "SAVED: ", filename
	
	#Spectra plot save code
	size = (1000,500)
	path = os.getenv('PWD')
	filenamelist = [path, '/', 'SpectraPlot_X', str(self.InspectorPosition[0]), '_Y', str(self.InspectorPosition[1]), '.png']
	filename = ''.join(filenamelist)
	container = self.Spectraplot1
	temp = container.outer_bounds
	container.outer_bounds = list(size)
	container.do_layout(force=True)
	gc = PlotGraphicsContext(size, dpi=DPI)
	gc.render_component(container)
	gc.save(filename)
	container.outer_bounds = temp
	print "SAVED: ", filename
	return

    
class SourcePlot(HasTraits):
    """Main class for the Source function plots.
    
    Uses inspector position from MainPlot to get position.
    
    Abstracting to alternate data types is not trivial.
	To expand to alternate type rewrite CalcSourcePlotVals, this method extracts all necessary data from the file and returns generic objects.
	
    Plot can be neglected simply by not entereding a filename inthe indata field in the FileLoadWindow upon startup.
    
    """
    

    PrimaryPlotC = Instance(HPlotContainer)
    windowwidth = 1000
    Wavelength = Float()
    ImageLock = Bool()
    SaveFlag = Bool()
    colormap_list = Enum('jet', 'gist_gray', 'gist_heat', 'Blues', 'hot')
    ColotList = ['green', 'red', 'yellow']
    TauOnelinecolor = Enum('red', 'green', 'yellow')
    Velocitylinecolor = Enum('green', 'red', 'white', 'yellow')
    VelocityZerolinecolor = Enum('yellow', 'red', 'white', 'green')
    Intensitylinecolor = Enum('white', 'red', 'yellow', 'green')
    wave_ref = float(FileLoadWindow.Wave_Ref)
    
    
    
    traits_view = View(VGroup(Group(Item('PrimaryPlotC', editor = ComponentEditor(), show_label = False, springy = True)),
			      HGroup(Item('colormap_list', style = 'simple', label = 'Colormap'), Item('TauOnelinecolor', label = 'Tau One Color', style = 'simple'), Item('Velocitylinecolor', label = 'Velocity Color', style = 'simple'), Item('VelocityZerolinecolor', label = 'Zero Velocity Color', style = 'simple'), Item('Intensitylinecolor', label = 'Intensity Color', style = 'simple'))),
			      width = windowwidth, height = windowwidth, resizable = True, title = "Source Function Plot")
    


    #Array used for passing information to speactra line plto and also Source Function Plots.  Will be synced
    InspectorPosition = Array()
    
    def __init__(self):
	super(SourcePlot, self).__init__()
	self.Wavelength = self.wave_ref
	self.RayFile = Rh15dout()
	self.RayFile.read_ray(FileLoadWindow.file_name_Ray)
	self.RayFile.read_indata(FileLoadWindow.file_name_Indata)
	#self.RayFile.read_ray('/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom/output_ray.ncdf')
	#self.RayFile.read_indata('/sanhome/tiago/rhout/cb24bih/MgII/PRD/385_PRD_newatom/output_indata.ncdf')
	WavelengthData = self.RayFile.ray.wavelength
	IntensityData = self.RayFile.ray.intensity
        

        
	self.create_SourcePlot()
	
    def CalcSourcePlotVals(self, sel, xi, yi, wave_ref = None, cscale=0.2, vmax=59,zrange=[-0.1,2],tle='',newfig=True, colour=False):
	# wave indices
	widx = sel.ray.wavelength_indices[:]
	# depth indices for temp, vz # this ONLY WORKS with non memmap arrays!!!!!!
	didx2 = sel.atmos.temperature[xi,yi] < 1.e30
	
	# protect against completely masked columns
	if not (float(N.sum(didx2)) > 0):
	    vz   = N.zeros(sel.atmos.temperature.shape[-1])
	    temp = vz.copy()
	    height = vz.copy()
	else:
	    vz = sel.atmos.velocity_z[xi,yi,didx2]
	    temp = sel.atmos.temperature[xi,yi,didx2]
	    height = sel.atmos.height[xi,yi,didx2]
	    #height = sel.atmos.height[490,26,didx2]
	
	if 'mask' in dir(sel.ray.intensity[xi,yi,0]):
	    if sel.ray.intensity[xi,yi,0].mask:
                return
		print('Intensity at (%i,%i) is masked! Making empty plot.' % (xi, yi))
		plot_form_diag_empty(height, vz, temp, vrange=[-vmax,vmax], cscale=cscale,zrange=zrange, tle=tle, newfig=newfig)
		return
	
	intensity = sel.ray.intensity[xi,yi,widx]
	wave = sel.ray.wavelength[widx]
	
	# have chi/source_function or full stuff?
	if hasattr(sel.ray, 'chi') and hasattr(sel.ray, 'source_function'):
	    # depth indices for Chi, S
	    didx1 = sel.ray.chi[xi,yi,:,0] > 0.
	    chi = sel.ray.chi[xi,yi,didx1,:]
	    S = sel.ray.source_function[xi,yi,didx1,:]
	elif hasattr(sel.ray, 'eta_line') and hasattr(sel.ray, 'chi_line'):
	    # depth indices for Chi, S
	    didx1 = sel.ray.chi_line[xi,yi,:,0] > 0.        
	    chi = sel.ray.chi_line[xi,yi,didx1,:]
	    S = sel.ray.eta_line[xi,yi,didx1,:]/sel.ray.chi_line[xi,yi,didx1,:]
	else:
	    raise ValueError('(EEE) form_diag_rr: structure has no suitable source function/chi.')
	
	
	if wave_ref is None:
	    wave_ref = N.mean(sel.ray.wavelength[widx])    
    
	#Break in old funcs
	
	vrange=[-vmax,vmax]
	sint = intensity
	# Calculate tau
	tau = tau_integ(chi, height)
	tau_one = get_tau_one(tau,height)
	
	# Convert to more convenient units
	height = height/1.e6 # Mm
	vz     = -vz/1e3     # km/s, opposite positive convention 
	
	vaxis = 299792.45 * (wave - wave_ref)/wave_ref
       
	if wave_ref > N.max(wave) or wave_ref < N.min(wave):
	    raise ValueError("Reference wavelength not contained in input wave!")
    
	wi = max(N.where(vaxis >= vrange[0])[0][0]  - 1, 0)
	wf = min(N.where(vaxis <= vrange[1])[0][-1] + 1, vaxis.shape[0]   - 1)
     
	zi = max(N.where(height <= zrange[1])[0][0]  - 1, 0)
	zf = min(N.where(height >= zrange[0])[0][-1] + 1, height.shape[0] - 1)
	
	# slice the arrays to improve the scaling
	tau     = tau[wi:wf+1,zi:zf+1]
	chi     = N.transpose(chi)[wi:wf+1,zi:zf+1]
	S       = N.transpose(S)[wi:wf+1,zi:zf+1]
	height  = height[zi:zf+1]
	vz      = vz[zi:zf+1]
	T       = temp[zi:zf+1]
	wave    = wave[wi:wf+1]
	bint    = int2bt(sint[wi:wf+1], wave)/1e3 # intensity in brightness temp [kK]
	vaxis   = vaxis[wi:wf+1]
	tau_one = tau_one[wi:wf+1]/1e6
	
	#find wavelength point closest to maximum in tau --NOPE. Now using v=0
	#ind   = N.argmax(tau_one)
	ind   = N.argmin(N.abs(vaxis))
	Sline = int2bt(S[ind],wave[ind]) # source function at ind, in brightness temp.
	
	#Flip arrays for plotting
	height = height[::-1]
	chiCT = chi[:, ::-1]
	tauCT = tau[:, ::-1]
	vz = vz[::-1]
	T = T[::-1]
	Sline = Sline[::-1]
	
    
	# colours:
	if colour:
	    # tau=1, vel, planck, S
	    ct = ['c-', 'r-', 'w:', 'y--' ]
	else:
	    ct = ['0.5', 'w-', 'w:', 'w--']
    
	color_source = Instance(ImageData)
	
	return {'height':height, 'chiCT':chiCT, 'tauCT':tauCT, 'tau':tau, 'chi':chi, 'vz':vz, 'T':T, 'Sline':Sline, 'vaxis':vaxis, 'tau_one':tau_one, 'bint':bint, 'wave':wave, 'vrange':vrange, 'zrange':zrange, 'S':S}
	
    def create_SourcePlot(self):

	plotvals = self.CalcSourcePlotVals(sel = self.RayFile, xi = 0, yi = 0, wave_ref = self.wave_ref, vmax = float(FileLoadWindow.VMax_Min), zrange=[float(FileLoadWindow.ZRangeMin), float(FileLoadWindow.ZRangeMax)])
        #plotvals = self.CalcSourcePlotVals(sel = self.RayFile, xi = 0, yi = 0, wave_ref = self.wave_ref,  zrange=[0.2,3.0])
	
	#unpacks values from CalcSourcePlotVals
	vaxis = plotvals['vaxis']
	height = plotvals['height']
	chiCT = plotvals['chiCT']
	tauCT = plotvals['tauCT']
	tau = plotvals['tau']
	chi = plotvals['chi']
	vz = plotvals['vz']
	T = plotvals['T']
	Sline = plotvals['Sline']
	tau_one = plotvals['tau_one']
	bint = plotvals['bint']
	wave = plotvals['wave']
	vrange = plotvals['vrange']
	zrange = plotvals['zrange']
	S = plotvals['S']

	
	velocityzero = 299792.45 * (self.Wavelength - self.wave_ref)/self.wave_ref
	
	
	vaxispoints = 100
	heightaxispoints = 400
        self.new_vaxis = N.linspace(vrange[0], vrange[1], vaxispoints)
        self.new_height = N.linspace(zrange[0], zrange[1], heightaxispoints)
	VelocityZeroXPoints = [velocityzero, velocityzero]
	VelocityZeroYPoints = [self.new_height[0], self.new_height[-1]]
	
	#----------------------Create Chi/Tau Plot-------------------------------
	new_chiCT = lineformInterpolate2D(vaxis, height, chiCT, self.new_vaxis , self.new_height)
        new_tauCT = lineformInterpolate2D(vaxis, height, tauCT, self.new_vaxis , self.new_height)


	self._image_indexCT = GridDataSource(xdata = self.new_vaxis, ydata = self.new_height[1:]) 
        index_mapperCT = GridMapper(range = DataRange2D(self._image_indexCT))
        self._image_valueCT = ImageData(data =  (new_chiCT[:,1:]/new_tauCT[:,1:]), value_depth = 1)
        color_mapperCT = jet(DataRange1D(self._image_valueCT))

	self.ChiTauPlotData = ArrayPlotData(imagedata = (new_chiCT[:,1:]/new_tauCT[:,1:]), index = vaxis, TauOne= tau_one, Velocity = vz, Height = height, VelocityZeroXPoints = VelocityZeroXPoints, VelocityZeroYPoints = VelocityZeroYPoints)
        self.ChiTauPlot = Plot(self.ChiTauPlotData)
        self.ChiTauPlot.img_plot1 = self.ChiTauPlot.img_plot("imagedata", colormap = color_mapperCT, xbounds = (self.new_vaxis[0], self.new_vaxis[-1]), ybounds = (self.new_height[0], self.new_height[-1]))[0]
	self.ChiTauPlot.overlays.append(ZoomTool(self.ChiTauPlot.img_plot1))
	self.ChiTauPlot.tools.append(PanTool(self.ChiTauPlot.img_plot1, drag_button="right"))

        self.ChiTauPlot.TauOnePlot = self.ChiTauPlot.plot(("index","TauOne"), type = "line", color = "red", resizeable = True)
        self.ChiTauPlot.VelocityPlot = self.ChiTauPlot.plot(("Velocity","Height"), type = "line", color = "green", resizeable = True)
	self.ChiTauPlot.VelocityZeroPlot = self.ChiTauPlot.plot(("VelocityZeroXPoints","VelocityZeroYPoints"), type = "line", color = "yellow", resizeable = True)

        self.ChiTauPlot.title = "Chi/Tau Plot"
	self.ChiTauPlot.x_axis.title = "Velocity (km/s)"
	self.ChiTauPlot.y_axis.title = "Height (Mm)"
        self.ChiTauPlot.range2d.x_range.set_bounds(self.new_vaxis[0], self.new_vaxis[-1])
        self.ChiTauPlot.range2d.y_range.set_bounds(self.new_height[0], self.new_height[-1])

        self.ChiTauPlotC = OverlayPlotContainer()
        self.ChiTauPlotC.add(self.ChiTauPlot)
	
	
	#------------------Create source function plot-------------------------
	sourceSF = S[:, ::-1]
        
        new_sourceSF = lineformInterpolate2D(vaxis, height, sourceSF, self.new_vaxis , self.new_height)

        self._image_indexSF = GridDataSource(xdata = self.new_vaxis, ydata = self.new_height)     
        index_mapperSF = GridMapper(range = DataRange2D(self._image_indexSF))

        self._image_valueSF = ImageData(data  = new_sourceSF, value_depth = 1)
        color_mapperSF = jet(DataRange1D(self._image_valueSF))

        self.SourceFuncPlotData = ArrayPlotData(imagedata = new_sourceSF, colormap = color_mapperSF, index = vaxis, TauOne = tau_one, Velocity = vz, Height = height, T = T/1e3, VelocityZeroXPoints = VelocityZeroXPoints, VelocityZeroYPoints = VelocityZeroYPoints)
        self.SourceFuncPlot = Plot(self.SourceFuncPlotData)
        self.SourceFuncPlot.img_plot1 = self.SourceFuncPlot.img_plot("imagedata", colormap = color_mapperSF, xbounds = (self.new_vaxis[0], self.new_vaxis[-1]), ybounds = (self.new_height[0], self.new_height[-1]))[0]
	self.SourceFuncPlot.overlays.append(ZoomTool(self.SourceFuncPlot.img_plot1))
        self.SourceFuncPlot.tools.append(PanTool(self.SourceFuncPlot.img_plot1, drag_button="right"))
	
        self.SourceFuncPlot.TauOnePlot = self.SourceFuncPlot.plot(("index","TauOne"), type = "line", color = "red")
        self.SourceFuncPlot.VelocityPlot = self.SourceFuncPlot.plot(("Velocity","Height"), type = "line", color = "green", resizeable = True)
	self.SourceFuncPlot.VelocityZeroPlot = self.SourceFuncPlot.plot(("VelocityZeroXPoints","VelocityZeroYPoints"), type = "line", color = "yellow", resizeable = True)    
    

        self.TemperaturePlotData = ArrayPlotData(T = T/1e3, Sline = Sline/1e3, Height = height)
        self.TemperaturePlot = Plot(self.TemperaturePlotData)
        self.TemperaturePlot.TPlot = self.TemperaturePlot.plot(("T","Height"), type = "line", color = "white", resizeable = True, line_style = 'dot', origin = 'bottom right')
        self.TemperaturePlot.Slinelot = self.TemperaturePlot.plot(("Sline","Height"), type = "line", color = "white", resizeable = True, line_style = 'dash', origin = 'bottom right')
        self.TemperaturePlot.x_grid.visible = False
        self.TemperaturePlot.y_grid.visible = False
        self.TemperaturePlot.x_axis.orientation = 'top'

        self.TemperaturePlot.y_axis.visible = False
        self.TemperaturePlot.range2d.x_range.set_bounds(0,30)
        self.TemperaturePlot.x_axis.title = "T (kK)"
        self.TemperaturePlot.default_origin = 'bottom right'
	self.TemperaturePlot.y_mapper = self.SourceFuncPlot.y_mapper
        self.SourceFuncPlotC = OverlayPlotContainer()
        self.SourceFuncPlotC.add(self.SourceFuncPlot)
        self.SourceFuncPlotC.add(self.TemperaturePlot)


        self.SourceFuncPlot.title = "Source Function Plot"
	self.SourceFuncPlot.x_axis.title = "Velocity (km/s)"
	self.SourceFuncPlot.y_axis.title = "Height (Mm)"
        self.SourceFuncPlot.range2d.x_range.set_bounds(self.new_vaxis[0], self.new_vaxis[-1])
        self.SourceFuncPlot.range2d.y_range.set_bounds(self.new_height[0], self.new_height[-1])
	
	
	#---------------------Create Tau * Exp(-tau) plot---------------------
	tauCalcTE = (tau*N.exp(-tau))
        tauCalcTE = tauCalcTE[:, ::-1]
    
        new_tauCalcTE = lineformInterpolate2D(vaxis, height, tauCalcTE, self.new_vaxis, self.new_height)
    
    
        self._image_indexTE = GridDataSource(xdata = self.new_vaxis, ydata = self.new_height)     
        index_mapperTE = GridMapper(range = DataRange2D(self._image_indexTE))
        self._image_valueTE = ImageData(data  = new_tauCalcTE, value_depth = 1)
        color_mapperTE = jet(DataRange1D(self._image_valueTE))
    
        self.TauExpPlotData = ArrayPlotData(imagedata = new_tauCalcTE, colormap = color_mapperTE, index = vaxis, TauOne = tau_one, Velocity = vz, Height = height, VelocityZeroXPoints = VelocityZeroXPoints, VelocityZeroYPoints = VelocityZeroYPoints)
        self.TauExpPlot = Plot(self.TauExpPlotData)
        self.TauExpPlot.img_plot1 = self.TauExpPlot.img_plot("imagedata", colormap = color_mapperTE, xbounds = (self.new_vaxis[0], self.new_vaxis[-1]), ybounds = (self.new_height[0], self.new_height[-1]))[0]
        self.TauExpPlot.overlays.append(ZoomTool(self.TauExpPlot.img_plot1))
	self.TauExpPlot.tools.append(PanTool(self.TauExpPlot.img_plot1, drag_button="right"))

	
        self.TauExpPlot.TauOnePlot = self.TauExpPlot.plot(("index","TauOne"), type = "line", color = "red", resizeable = True)
        self.TauExpPlot.VelocityPlot = self.TauExpPlot.plot(("Velocity","Height"), type = "line", color = "green", resizeable = True)
        self.TauExpPlot.VelocityZeroPlot = self.TauExpPlot.plot(("VelocityZeroXPoints","VelocityZeroYPoints"), type = "line", color = "yellow", resizeable = True)
	
        self.TauExpPlot.title = "Tau * Exp Plot"
	self.TauExpPlot.x_axis.title = "Velocity (km/s)"
	self.TauExpPlot.y_axis.title = "Height (Mm)"
        self.TauExpPlot.range2d.x_range.set_bounds(self.new_vaxis[0], self.new_vaxis[-1])
        self.TauExpPlot.range2d.y_range.set_bounds(self.new_height[0], self.new_height[-1])
    
        self.TauExpPlotC = OverlayPlotContainer()
        self.TauExpPlotC.add(self.TauExpPlot)

	#--------------------Create Contribution Function Plot------------------------
	new_tau_oneCF = lineformInterpolate1D(vaxis, tau_one, self.new_vaxis)
        new_vzCF = lineformInterpolate1D(height, vz, self.new_height)
        new_bint = lineformInterpolate1D(vaxis, bint, self.new_vaxis)
        
        ContributionFuncCF = (chi*N.exp(-tau)*S)[:,::-1]
        ContributionFuncCF /= N.max(ContributionFuncCF, axis=0)
        new_ContributionPlotCF = lineformInterpolate2D(vaxis, height, ContributionFuncCF, self.new_vaxis , self.new_height)
	
	
        self._image_indexCF = GridDataSource(xdata = self.new_vaxis, ydata = self.new_height)     
        index_mapperCF = GridMapper(range = DataRange2D(self._image_indexCF))
        self._image_valueCF = ImageData(data  = new_ContributionPlotCF, value_depth = 1)
        color_mapperCF = jet(DataRange1D(self._image_valueCF))
        
        self.ContributionFuncPlotData = ArrayPlotData(imagedata = new_ContributionPlotCF, colormap = color_mapperCF, index = vaxis, TauOne = tau_one, Velocity = vz, Height = height, VelocityZeroXPoints = VelocityZeroXPoints, VelocityZeroYPoints = VelocityZeroYPoints)
        self.ContributionFuncPlot = Plot(self.ContributionFuncPlotData)
        self.ContributionFuncPlot.img_plot1 = self.ContributionFuncPlot.img_plot("imagedata", colormap = color_mapperCF, xbounds = (self.new_vaxis[0], self.new_vaxis[-1]), ybounds = (self.new_height[0], self.new_height[-1]))[0]
	self.ContributionFuncPlot.overlays.append(ZoomTool(self.ContributionFuncPlot.img_plot1))
	self.ContributionFuncPlot.tools.append(PanTool(self.ContributionFuncPlot.img_plot1, drag_button="right"))

        self.ContributionFuncPlot.TauOnePlot = self.ContributionFuncPlot.plot(("index","TauOne"), type = "line", color = "red", resizeable = True)
        self.ContributionFuncPlot.VelocityPlot = self.ContributionFuncPlot.plot(("Velocity","Height"), type = "line", color = "green", resizeable = True)
        self.ContributionFuncPlot.VelocityZeroPlot = self.ContributionFuncPlot.plot(("VelocityZeroXPoints","VelocityZeroYPoints"), type = "line", color = "yellow", resizeable = True)
    
        self.IntensityPlotData = ArrayPlotData(index = self.new_vaxis, bint = new_bint)
        self.IntensityPlot = Plot(self.IntensityPlotData)
        self.IntensityPlot.intensity = self.IntensityPlot.plot(("index","bint"), color = 'white', line_width = 2.0)
        self.IntensityPlot.y_axis.orientation = 'right'
        self.IntensityPlot.y_axis.title = "I_v (kK)"
        self.IntensityPlot.default_origin = 'bottom right'
        self.IntensityPlot.x_grid.visible = False
        self.IntensityPlot.y_grid.visible = False
        self.IntensityPlot.x_axis.visible = False
        self.IntensityPlot.x_axis = self.ContributionFuncPlot.x_axis
        self.IntensityPlot.range2d.y_range.set_bounds(3,7)
        #right_axis = LabelAxis(IntensityPlot, orientation = 'right')
        #IntensityPlot.underlays.append(right_axis)
	
	
        
        self.ContributionFuncPlot.title = "Contribution Function Plot"
	self.ContributionFuncPlot.x_axis.title = "Velocity (km/s)"
	self.ContributionFuncPlot.y_axis.title = "Height (Mm)"
        self.ContributionFuncPlot.range2d.x_range.set_bounds(self.new_vaxis[0], self.new_vaxis[-1])
        self.ContributionFuncPlot.range2d.y_range.set_bounds(self.new_height[0], self.new_height[-1])
	
        self.ContributionFuncPlotC = OverlayPlotContainer()
        self.ContributionFuncPlotC.add(self.ContributionFuncPlot)
        self.ContributionFuncPlotC.add(self.IntensityPlot)
	
	
	#-------------------------------Final Container Construction-------------------------------
	self.TauExpPlot.range2d = self.ChiTauPlot.range2d
	self.ContributionFuncPlot.range2d = self.ChiTauPlot.range2d
	self.SourceFuncPlot.range2d = self.ChiTauPlot.range2d
        self.LeftPlots = VPlotContainer(self.TauExpPlotC, self.ChiTauPlotC, background = "lightgray", use_back_buffer = True)
        self.LeftPlots.spacing = 0
        self.RightPlots = VPlotContainer(self.ContributionFuncPlotC, self.SourceFuncPlotC, background = "lightgray", use_back_buffer = True)
        self.RightPlots.spacing = 0
	
        MainContainer = HPlotContainer(self.LeftPlots, self.RightPlots, background = "lightgray", use_back_buffer = True)
        MainContainer.spacing = 0
        self.PrimaryPlotC = MainContainer 
	
    def _InspectorPosition_changed(self):
	if self.ImageLock == True:
	    return
	
	if len(self.InspectorPosition) > 0:
            xi = self.InspectorPosition[0]
            yi = self.InspectorPosition[1]
        else:
            xi = 0
            yi = 0
	
	plotvals = self.CalcSourcePlotVals(sel = self.RayFile, xi=xi, yi=yi, wave_ref = self.wave_ref, vmax = float(FileLoadWindow.VMax_Min), zrange=[float(FileLoadWindow.ZRangeMin), float(FileLoadWindow.ZRangeMax)])
	#plotvals = self.CalcSourcePlotVals(sel = self.RayFile, xi = xi, yi = yi, wave_ref = self.wave_ref, vmax = 59, zrange=[0.2,3.0])

	#unpacks values from CalcSourcePlotVals
	vaxis = plotvals['vaxis']
	height = plotvals['height']
	chiCT = plotvals['chiCT']
	tauCT = plotvals['tauCT']
	tau = plotvals['tau']
	chi = plotvals['chi']
	vz = plotvals['vz']
	T = plotvals['T']
	Sline = plotvals['Sline']
	tau_one = plotvals['tau_one']
	bint = plotvals['bint']
	wave = plotvals['wave']
	vrange = plotvals['vrange']
	zrange = plotvals['zrange']
	S = plotvals['S']
	
	#ChiTau Plot updates
	new_chiCT = lineformInterpolate2D(vaxis, height, chiCT, self.new_vaxis, self.new_height)
        new_tauCT = lineformInterpolate2D(vaxis, height, tauCT, self.new_vaxis , self.new_height)
        self.ChiTauPlotData.set_data("imagedata",(new_chiCT[:,1:]/new_tauCT[:,1:]))
        self.ChiTauPlotData.set_data("index", vaxis)
        self.ChiTauPlotData.set_data("TauOne", tau_one)
        self.ChiTauPlotData.set_data("Velocity",vz)
        self.ChiTauPlotData.set_data("Height", height)
	
	#Source Plot Updates
	sourceSF = S
	sourceSF = S[:, ::-1]
	new_sourceSF = lineformInterpolate2D(vaxis, height, sourceSF, self.new_vaxis, self.new_height)
	self.SourceFuncPlotData.set_data("imagedata", new_sourceSF)
	self.SourceFuncPlotData.set_data("index",vaxis)
	self.SourceFuncPlotData.set_data("TauOne",tau_one)
	self.SourceFuncPlotData.set_data("Velocity",vz)
	self.SourceFuncPlotData.set_data("Height",height)
	self.TemperaturePlotData.set_data("T", T/1.0e3)
	self.TemperaturePlotData.set_data("Sline",Sline/1e3)
	self.TemperaturePlotData.set_data("Height",height)
	 
	#Tau Plot updates
	tauCalcTE = (tau*N.exp(-tau))
	tauCalcTE = tauCalcTE[:, ::-1]
	new_tauCalcTE = lineformInterpolate2D(vaxis, height, tauCalcTE, self.new_vaxis, self.new_height)     
	self.TauExpPlotData.set_data("imagedata", new_tauCalcTE)
	self.TauExpPlotData.set_data("index",vaxis)
	self.TauExpPlotData.set_data("TauOne",tau_one)
	self.TauExpPlotData.set_data("Velocity",vz)
	self.TauExpPlotData.set_data("Height",height)
	 
	#Source Function Plot Updates
	new_tau_oneCF = lineformInterpolate1D(vaxis, tau_one, self.new_vaxis)
	new_vzCF = lineformInterpolate1D(height, vz, self.new_height)
	new_bint = lineformInterpolate1D(vaxis, bint, self.new_vaxis)
	ContributionFuncCF = (chi*N.exp(-tau)*S)[:,::-1]
	ContributionFuncCF /= N.max(ContributionFuncCF, axis=0)
	new_ContributionPlotCF = lineformInterpolate2D(vaxis, height, ContributionFuncCF, self.new_vaxis, self.new_height)
	self.new_vaxis = N.linspace(vrange[0], vrange[1], 100)
	self.new_height = N.linspace(zrange[0], zrange[1], 400)
	new_tau_oneCF = lineformInterpolate1D(vaxis, tau_one, self.new_vaxis)
	new_vzCF = lineformInterpolate1D(height, vz, self.new_height)
	new_bint = lineformInterpolate1D(vaxis, bint, self.new_vaxis)
	self.ContributionFuncPlotData.set_data("imagedata", new_ContributionPlotCF)
	self.ContributionFuncPlotData.set_data("index",vaxis)
	self.ContributionFuncPlotData.set_data("TauOne",tau_one)
	self.ContributionFuncPlotData.set_data("Velocity",vz)
	self.ContributionFuncPlotData.set_data("Height",height)
	self.IntensityPlotData.set_data("index",self.new_vaxis)
	self.IntensityPlotData.set_data("bint",new_bint)
	
    def _Wavelength_changed(self):
	velocityzero = 299792.45 * (self.Wavelength - self.wave_ref)/self.wave_ref
	VelocityZeroXPoints = [velocityzero, velocityzero]
	try:
	    self.ChiTauPlotData.set_data("VelocityZeroXPoints", VelocityZeroXPoints)
	    self.SourceFuncPlotData.set_data("VelocityZeroXPoints", VelocityZeroXPoints)
	    self.TauExpPlotData.set_data("VelocityZeroXPoints", VelocityZeroXPoints)
	    self.ContributionFuncPlotData.set_data("VelocityZeroXPoints", VelocityZeroXPoints)
	except:
	    pass

    def _colormap_list_changed(self):
	plotlist = [self.ChiTauPlot.img_plot1, self.SourceFuncPlot.img_plot1, self.TauExpPlot.img_plot1, self.ContributionFuncPlot.img_plot1]

	for x in plotlist:
	    clr_range = x.color_mapper.range
	    color_mapper = eval(self.colormap_list)(clr_range)
	    x.color_mapper = color_mapper

	self.ChiTauPlot.request_redraw()
    
    def _TauOnelinecolor_changed(self):
	self.ChiTauPlot.TauOnePlot[0].color = self.TauOnelinecolor
	self.SourceFuncPlot.TauOnePlot[0].color = self.TauOnelinecolor
	self.TauExpPlot.TauOnePlot[0].color = self.TauOnelinecolor
	self.ContributionFuncPlot.TauOnePlot[0].color = self.TauOnelinecolor
	
    def _Velocitylinecolor_changed(self):
	self.ChiTauPlot.VelocityPlot[0].color = self.Velocitylinecolor
	self.SourceFuncPlot.VelocityPlot[0].color = self.Velocitylinecolor
	self.TauExpPlot.VelocityPlot[0].color = self.Velocitylinecolor
	self.ContributionFuncPlot.VelocityPlot[0].color = self.Velocitylinecolor
	
    def _VelocityZerolinecolor_changed(self):
	self.ChiTauPlot.VelocityZeroPlot[0].color = self.VelocityZerolinecolor
	self.SourceFuncPlot.VelocityZeroPlot[0].color = self.VelocityZerolinecolor
	self.TauExpPlot.VelocityZeroPlot[0].color = self.VelocityZerolinecolor
	self.ContributionFuncPlot.VelocityZeroPlot[0].color = self.VelocityZerolinecolor

    def _Intensitylinecolor_changed(self):
	self.IntensityPlot.intensity[0].color = self.Intensitylinecolor
	
    def _SaveFlag_changed(self):
	DPI = 72
	#Main plot save code
	size = (1000,1000)
	path = os.getenv('PWD')
	filenamelist = [path, '/', 'SourcePlots_X', str(self.InspectorPosition[0]), '_Y', str(self.InspectorPosition[1]), '.png']
	filename = ''.join(filenamelist)
	container = self.PrimaryPlotC
	temp = container.outer_bounds
	container.outer_bounds = list(size)
	container.do_layout(force=True)
	gc = PlotGraphicsContext(size, dpi=DPI)
	gc.render_component(container)
	gc.save(filename)
	container.outer_bounds = temp
	print "SAVED: ", filename

if __name__ == "__main__":
    MainPlot = MainPlot()
    
    if FileLoadWindow.file_name_Indata == '':
	MainPlot.configure_traits()
    
    if FileLoadWindow.file_name_Ray != '' and FileLoadWindow.file_name_Indata != '':
	SourcePlot = SourcePlot()
	MainPlot.edit_traits()
	
	MainPlot.sync_trait('InspectorPosition', SourcePlot, 'InspectorPosition')
	MainPlot.sync_trait('ImageLock',SourcePlot, 'ImageLock')
	MainPlot.sync_trait('Wavelength', SourcePlot, 'Wavelength') 
	MainPlot.sync_trait('Wavelength',MainPlot,'WavelengthFine')
	MainPlot.sync_trait('SaveFlag', SourcePlot, 'SaveFlag')
	

	SourcePlot.configure_traits()
    
    
    	


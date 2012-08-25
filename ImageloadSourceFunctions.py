#!/usr/bin/env python
#Attempt to model on 2d function inspector example

from scipy import *
#from pylab import *
from numpy import *
import numpy as N

import sys

from numpy import linspace, meshgrid, pi, sin, array

# Enthought library imports
from enable.api import Component, ComponentEditor, BaseTool, Window, LineStyle
from enable.component_editor import ComponentEditor

from traits.api import HasTraits, Instance, DelegatesTo, on_trait_change, Trait, Event, Array, Any, Callable, Int, Tuple, Enum, Bool
from traitsui.api import Item, Group, View, HGroup, VGroup, Handler
from traits.util.resource import find_resource
#from enthought.chaco.tools.cursor_tool import CursorTool, BaseCursorTool

# Chaco imports
from chaco.api import ArrayPlotData, Plot, Blues, HPlotContainer, VPlotContainer, create_line_plot, GridDataSource, ImageData, Plot, PlotAxis, LinearMapper, GridMapper, DataRange2D, DataRange1D, CMapImagePlot, ArrayDataSource, DataView, ContourPolyPlot, CMapImagePlot, ColorBar, PlotLabel, PlotAxis, ImagePlot, OverlayPlotContainer, add_default_axes, LabelAxis
from chaco.tools.api import PanTool, ZoomTool, LineInspector
from scipy.interpolate import interp1d
from scipy.integrate.quadrature import cumtrapz
from chaco.default_colormaps import *


# RH code import functions
#from tt.lines import *

#Import functions for linefom/plot_form_diag
#from tt.math.utilsmath import int2bt
#import matplotlib.lines as mpllines

#Interpolate functions
from scipy import interpolate






#-----------------------------------------------------------------------

def int2bt(inu, w):
    ''' Converts from radiation intensity (in J s^-1 m^-2  Hz^-1 sr^-1 units)
        to brightness temperature units (in K), at a given wavelength wave
        (in nm). '''
    
    c = 2.99792458e8
    h = 6.62606896e-34
    k = 1.3806504e-23
    wave = w * 1.e-9 # wavelength to m

    return h*c / (wave * k * N.log(2*h*c/(wave**3*inu) + 1) )

#-----------------------------------------------------------------------


    
    

#----------------------------------------------------------------------------


def tau_integ(chi_in, height):
    ''' Integrates the opacity to get the optical depth tau '''
    
    chi = N.transpose(chi_in)
    tau = N.zeros(chi.shape)
    
    try:
        # if there is a masked array, sum only over used depth points
        mm = N.invert(chi[0].mask)
        zcut = N.where(mm == 1)[0][0]
        tau[:,zcut+1:] = cumtrapz(chi[:,mm], x=-height[mm]) 
    except AttributeError:
        tau[:,1:] = cumtrapz(chi, x=-height) 
    
    
    return tau

#----------------------------------------------------------------------------
"""
def get_tau_one(tau, height):
    ''' Finds the tau = 0 iso curve from a 2D array of tau, and the height array.
        Height should be the last index on the tau array. '''
    
    tau_one = N.zeros(tau.shape[0])
    
    for i in range(tau.shape[0]):
        f = interp1d(tau[i],height)
        tau_one[i] = f(1.0)
    
    return tau_one
"""

def get_tau_one(tau, height):
    ''' Finds the tau = 0 iso curve from a 2D array of tau,
        and the height array. Height should be the last index
        on the tau array. '''
    
    idx = N.argmin(N.abs(tau-1.), axis=1)
    tau0 = tau[(N.arange(tau.shape[0]), idx)]
    tau1 = tau[(N.arange(tau.shape[0]), idx+1)]
    height0 = height[idx]
    height1 = height[idx+1]
    
    # manual linear inter/extrapolation
    theta = (1.-tau0) / (tau1-tau0)
    tau_one = height0 + theta * (height1 - height0)
    
    return tau_one

#----------------------------------------------------------------------------
def lineformInterpolate2D(vaxisF, heightF, valuesF,  new_vaxisF, new_heightF):

    fint = interpolate.interp1d(heightF, valuesF, kind='linear',bounds_error=False, fill_value=0)
    temp = fint(new_heightF)
    fint = interpolate.interp1d(vaxisF, temp.T, kind='linear',bounds_error=False, fill_value=0)
    new_valuesF = fint(new_vaxisF)

    return new_valuesF

def lineformInterpolate1D(vaxisF, valuesF, new_vaxisF):
    fint = interpolate.interp1d(vaxisF, valuesF, kind = 'linear',bounds_error=False, fill_value=0)
    new_valuesF = fint(new_vaxisF)
    return new_valuesF


#-----------------------------------------------------------------------------





def form_diag_rr(sel, xi, yi, wave_ref = None, cscale=0.2, vmax=59,zrange=[-0.1,2],tle='',newfig=True, colour=False):
	#Calls plot_form_diag from Rh15dout object. TO BE PUT AS METHOD IN rh15d.py


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
	
	if 'mask' in dir(sel.ray.intensity[xi,yi,0]):
		if sel.ray.intensity[xi,yi,0].mask:
			print('Intensity at (%i,%i) is masked! Making empty plot.' % (xi, yi))
			plot_form_diag_empty(height, vz, temp, vrange=[-vmax,vmax], cscale=cscale,
								 zrange=zrange, tle=tle, newfig=newfig)
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
	

	
	plot_PrimaryPlotC(chi, S, wave, height, vz, temp, intensity, vrange=[-vmax,vmax],colour=colour,
				   wave_ref=wave_ref, cscale=cscale, zrange=zrange, tle=tle, newfig=newfig)




#####################################################################
""" Defines the ImageInspectorTool, ImageInspectorOverlay, and
    ImageInspectorColorbarOverlay classes.
    """
# Enthought library imports
from enable.api import BaseTool, KeySpec
from traits.api import Any, Bool, Enum, Event, Tuple

# Chaco imports
from chaco.api import AbstractOverlay, ImagePlot, TextBoxOverlay


class ImageInspectorTool(BaseTool):
    """ A tool that captures the color and underlying values of an image plot.
        """
    ImageLock = Bool()
    #ImageLock = False
    SaveFlag = Bool()
    
    # This event fires whenever the mouse moves over a new image point.
    # Its value is a dict with a key "color_value", and possibly a key
    # "data_value" if the plot is a color-mapped image plot.
    new_value = Event
    
    # Indicates whether overlays listening to this tool should be visible.
    visible = Bool(True)
    
    # Stores the last mouse position.  This can be used by overlays to
    # position themselves around the mouse.
    last_mouse_position = Tuple
    
    # This key will show and hide any ImageInspectorOverlays associated
    # with this tool.
    inspector_key = KeySpec('p')
    ImageLock_key = KeySpec('x')
    Save_key = KeySpec('s')
    
    # Stores the value of self.visible when the mouse leaves the tool,
    # so that it can be restored when the mouse enters again.
    _old_visible = Enum(None, True, False) #Trait(None, Bool(True))
    
    def normal_key_pressed(self, event):
        if self.inspector_key.match(event):
            self.visible = not self.visible
            event.handled = True
        
	if self.ImageLock_key.match(event):
	    if self.ImageLock == True:
		self.ImageLock = False
		event.handled = True
	    elif self.ImageLock == False:
		self.ImageLock = True
		event.handled = True
	
	if self.Save_key.match(event):
	    if self.SaveFlag == True:
		self.SaveFlag = False
		event.handled = True
	    elif self.SaveFlag == False:
		self.SaveFlag = True
		event.handled = True

    
    def normal_mouse_leave(self, event):
        if self._old_visible is None:
            self._old_visible = self.visible
            self.visible = False
    
    def normal_mouse_enter(self, event):
        if self._old_visible is not None:
            self.visible = self._old_visible
            self._old_visible = None
    
    def normal_mouse_move(self, event):
        """ Handles the mouse being moved.
            
            Fires the **new_value** event with the data (if any) from the event's
            position.
            """
        plot = self.component
        if plot is not None:
            if isinstance(plot, ImagePlot):
                ndx = plot.map_index((event.x, event.y))
                if ndx == (None, None):
                    self.new_value = None
                    return
                
                x_index, y_index = ndx
                image_data = plot.value
                if hasattr(plot, "_cached_mapped_image") and \
                    plot._cached_mapped_image is not None:
                        self.new_value = \
                            dict(indices=ndx,
                                 data_value=image_data.data[y_index, x_index],
                                 color_value=plot._cached_mapped_image[y_index,
                                                                       x_index])
                
                else:
                    self.new_value = \
                        dict(indices=ndx,
                             color_value=image_data.data[y_index, x_index])
                
                self.last_mouse_position = (event.x, event.y)
        return


class ImageInspectorOverlay(TextBoxOverlay):
    """ An overlay that displays a box containing values from an
        ImageInspectorTool instance.
        """
    InspectorPosition=Array()
    xi = Int()
    yi = Int()
    
    # An instance of ImageInspectorTool; this overlay listens to the tool
    # for changes, and updates its displayed text accordingly.
    image_inspector = Any
    
    # Anchor the text to the mouse?  (If False, then the text is in one of the
    # corners.)  Use the **align** trait to determine which corner.
    tooltip_mode = Bool(False)
    
    # The default state of the overlay is invisible (overrides PlotComponent).
    visible = False
    
    # Whether the overlay should auto-hide and auto-show based on the
    # tool's location, or whether it should be forced to be hidden or visible.
    visibility = Enum("auto", True, False)
    
    def _image_inspector_changed(self, old, new):
        if old:
            old.on_trait_event(self._new_value_updated, 'new_value', remove=True)
            old.on_trait_change(self._tool_visible_changed, "visible", remove=True)
        if new:
            new.on_trait_event(self._new_value_updated, 'new_value')
            new.on_trait_change(self._tool_visible_changed, "visible")
            self._tool_visible_changed()
    
    def _new_value_updated(self, event):
        if event is None:
            self.text = ""
            if self.visibility == "auto":
                self.visible = False
            return
        elif self.visibility == "auto":
            self.visible = True
        
        if self.tooltip_mode:
            self.alternate_position = self.image_inspector.last_mouse_position
        else:
            self.alternate_position = None
        
        d = event
        newstring = ""
        if 'indices' in d:
            newstring += '(%d, %d)' % d['indices'] + '\n'
        if 'color_value' in d:
            newstring += "(%d, %d, %d)" % tuple(map(int,d['color_value'][:3])) + "\n"
        if 'data_value' in d:
            newstring += str(d['data_value'])
        if 'indices' in d:
            self.InspectorPosition = d['indices']
            self.xi, self.yi = d['indices']
        #print "indices stored"
        
        self.text = newstring
        self.component.request_redraw()
    
    
    def _visible_changed(self):
        self.component.request_redraw()
    
    def _tool_visible_changed(self):
        self.visibility = self.image_inspector.visible
        if self.visibility != "auto":
            self.visible = self.visibility


class ImageInspectorColorbarOverlay(AbstractOverlay):
    pass

####################################################################

    

   
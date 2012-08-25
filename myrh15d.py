# Set of programs and tools to read the outputs from RH, 1.5D version 

import os
import numpy as N
import netCDF4 as nc

class Rh15dout:
    def __init__(self,fdir='.',verbose=True):

        self.files = []
        self.params = {}
        self.verbose = verbose
        self.fdir = fdir

        if os.path.isfile('%s/output_aux.ncdf' % self.fdir):
            self.read_aux()

        if os.path.isfile('%s/output_indata.ncdf' % self.fdir):
            self.read_indata()

        if os.path.isfile('%s/output_spectrum.ncdf' % self.fdir):
            self.read_spectrum()       
            
        if os.path.isfile('%s/output_ray.ncdf' % self.fdir):
            self.read_ray()            
        
    #------------------------------------------------------------------------

    def read_aux(self, infile=None):
        ''' Reads Aux file. '''

        if infile is None: infile = '%s/output_aux.ncdf' % self.fdir

        self.files.append(read_ncdf(self, infile))
        if self.verbose: print('--- Read %s file.' % infile)

        return

    #------------------------------------------------------------------------

    def read_indata(self, infile=None):
        ''' Reads indata file. '''

        if infile is None: infile = '%s/output_indata.ncdf' % self.fdir
        
        self.files.append(read_ncdf(self, infile))
        if self.verbose: print('--- Read %s file.' % infile)

        return
        
    #------------------------------------------------------------------------

    def read_spectrum(self, infile=None):
        ''' Reads spectrum file. '''

        if infile is None: infile = '%s/output_spectrum.ncdf' % self.fdir
        
        self.spectrum = DataHolder()
        self.files.append(read_ncdf(self.spectrum, infile))
        if self.verbose: print('--- Read %s file.' % infile)

        return

    #------------------------------------------------------------------------

    def read_ray(self, infile=None):
        ''' Reads ray file. '''
        
        if infile is None: infile = '%s/output_ray.ncdf' % self.fdir

        self.ray = DataHolder()
        self.files.append(read_ncdf(self.ray, infile))
        if self.verbose: print('--- Read %s file.' % infile)

        return

    #------------------------------------------------------------------------

    def read_J(self, infile='scratch/J.dat.ncdf'):
        ''' Reads angle averaged intensity file '''

        self.files.append(read_ncdf(self, infile))
        if self.verbose: print('--- Read %s file.' % infile)

        return

    #------------------------------------------------------------------------

    def close(self):
        ''' Closes the open NetCDF files '''

        for f in self.files:
            f.close()
            
#----------------------------------------------------------------------------



"""
class NcdfAtmos:
    def __init__(self, infile):
        self.file = read_ncdf(self, infile)
        self.closed = False
        
    #------------------------------------------------------------------------
    
    def close(self):
        try:
            self.file.close()
            self.closed = True
        except RuntimeError:
            print('(WWW) NcdfAtmos: input file already closed.')
            
    #------------------------------------------------------------------------
        
    def read(self, infile):
        if not self.closed: self.close()
        
        self.file = read_ncdf(self, infile)
        return

    #------------------------------------------------------------------------
    
    def write_multi(self, outfile, xi, yi, nti=0, writeB=False, write_dscale=False,
                    zcut=0, depth_optimise=False):
        ''' Writes MULTI atmosphere file from a column of the 3D model, in RH 1.5D ncdf 
        format. Also writes the binary XDR file with magnetic fields, if writeB
        is true. '''
        

        from tt.lines.multi import watmos_multi
        from tt.lines.rh    import write_B       
    
        writeB = writeB and indata.params['has_B']
    
        # if only total H available, will have to use rhpy (which is sometimes risky...)
        if self.params['nhydr'] == 1:
            import rhpy
            nh = rhpy.nh_lte(self.temperature[nti,xi,yi,zcut:].astype('Float64'),
                             self.electron_density[nti,xi,yi,zcut:].astype('Float64'),
                             self.hydrogen_populations[nti,0,xi,yi,zcut:].astype('Float64'))
        elif self.params['nhydr'] == 6:
            nh = self.hydrogen_populations[nti,:,xi,yi,zcut:]
        else:
            raise ValueError('(EEE) write_multi: found %i hydrogen levels. For multi, need 6 or 1 '
                             % self.params['nhydr'])
        
        M_TO_CM3 = (100.)**3
        M_TO_KM = 0.001
        
        temp = self.temperature[nti,xi,yi,zcut:].copy()
        ne = self.electron_density[nti,xi,yi,zcut:].copy()/M_TO_CM3
        z = self.z[nti,zcut:].copy()*M_TO_KM*1.e5    # in cm
        vz = self.velocity_z[nti,xi,yi,zcut:].copy()*M_TO_KM
        nh = nh/M_TO_CM3
        
        if writeB:
            bx = self.B_x[nti,xi,yi,zcut:].copy()
            by = self.B_y[nti,xi,yi,zcut:].copy()
            bz = self.B_z[nti,xi,yi,zcut:].copy()
        else:
            bx = by = bz = None
        
        if depth_optimise:
            rho = self.hydrogen_populations[nti,0,xi,yi,zcut:] * 2.380491e-24 / M_TO_CM3
            res = depth_optim(z, temp, ne, vz, rho, nh=nh, bx=bx, by=by, bz=bz)
            z, temp, ne, vz, rho, nh = res[:6]
            if writeB:
                bx, by, bz = res[6:]
                
        watmos_multi(outfile, temp, ne, z*1e-5, vz=vz, nh=nh, write_dscale=write_dscale,
                              id = '%s txy-slice: (t,x,y) = (%i,%i,%i)' % 
                                   (self.params['description'],nti, xi, yi))
    
    
        if writeB:
            write_B('%s.B' % bx, by, bz)
            print('--- Wrote magnetic field to %s.B' % outfile)        
        
        return
    
    #------------------------------------------------------------------------
    
    def write_multi_3d(self, outfile, nti=0, sx=None, sy=None, sz=None, big_endian=False):
        ''' Writes atmosphere in multi_3d format (the same as the pre-Jorrit multi3d) '''
        
        from tt.lines import multi
        
        ul = 1e2  # m to cm
        uv = 1e-3 # m/s to km/s
        
        # slicing and unit conversion
        if sx is None: sx = [0, self.nx, 1]
        if sy is None: sy = [0, self.ny, 1]
        if sz is None: sz = [0, self.nz, 1]
        
        if self.params['nhydr'] > 1:
          nh=N.mean(self.hydrogen_populations[nti,:,sx[0]:sx[1]:sx[2],
                                                    sy[0]:sy[1]:sy[2],
                                                    sz[0]:sz[1]:sz[2]],axis=1)/(ul**3)
        else:
          nh=self.hydrogen_populations[nti,0,sx[0]:sx[1]:sx[2],
                                             sy[0]:sy[1]:sy[2],
                                             sz[0]:sz[1]:sz[2]]/(ul**3)
          
        rho  = nh*2.380491e-24 # nH to rho [g cm-3]
        
        x    = self.x[sx[0]:sx[1]:sx[2]]*ul
        y    = self.y[sy[0]:sy[1]:sy[2]]*ul
        z    = self.z[nti, sz[0]:sz[1]:sz[2]]*ul
        
        ne   = self.electron_density[nti,sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                         sz[0]:sz[1]:sz[2]] / (ul**3)
        temp =      self.temperature[nti,sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                         sz[0]:sz[1]:sz[2]]
        vz   =       self.velocity_z[nti,sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                         sz[0]:sz[1]:sz[2]] * uv
        
        # write to file
        multi.write_atmos3d(outfile,x,y,z,ne,temp,vz,rho=rho,big_endian=big_endian)
        
        return
        

#############################################################################
###   TOOLS                                                               ###
#############################################################################
"""
class DataHolder:
    def __init__(self):
        pass

#----------------------------------------------------------------------------

def read_ncdf(inclass, infile):
    ''' Reads NetCDF file into inclass, instance of any class.
        Variables are read into class attributes, dimensions and attributes
        are read into params dictionary. '''

    # internal attributes of NetCDF groups
    ncdf_internals =['__class__', '__delattr__', '__doc__', '__format__', 
                     '__getattr__', '__getattribute__', '__hash__', '__init__',
                     '__new__', '__reduce__', '__reduce_ex__', '__repr__',
                     '__setattr__', '__sizeof__', '__str__', '__subclasshook__',
                     '_enddef', '_grpid', '_redef', 'close', 'cmptypes',
                     'createCompoundType', 'createDimension', 'createGroup',
                     'createVLType', 'createVariable', 'delncattr', 'dimensions',
                     'file_format', 'getncattr', 'groups', 'maskanscale', 
                     'ncattrs', 'parent', 'path', 'renameDimension', 
                     'renameVariable', 'set_fill_off', 'set_fill_on', 'setncattr',
                     'sync', 'variables', 'vltypes']

    if not os.path.isfile(infile):
        raise IOError('read_ncdf: File %s not found' % infile)

    f = nc.Dataset(infile, mode='r')

    if 'params' not in dir(inclass):
        inclass.params = {}
    
    # add dimensions as attributes
    for d in f.dimensions.keys():
        inclass.params[d] = len(f.dimensions[d])
        
    # add attributes
    attrs = [a for a in dir(f) if a not in ncdf_internals]
    for att in attrs:
        inclass.params[att] = getattr(f, att)

    # add variables
    for v in f.variables.keys():
        vname = v.replace(' ','_')    # sanitise string for spaces
        setattr(inclass, vname, f.variables[v])

    # Now do the same for all groups
    for group in f.groups.keys():

        gname = group.replace(' ','_') # sanitise string for spaces
        setattr(inclass, gname, DataHolder())
        
        cur_group = f.groups[group]
        cur_class = getattr(inclass, gname)
        
        # add variables
        for v in cur_group.variables.keys():
            vname = v.replace(' ','_') # sanitise string for spaces
            setattr(cur_class, vname, cur_group.variables[v])

        # add dimensions as attributes
        for d in cur_group.dimensions.keys():
            inclass.params[d] = len(cur_group.dimensions[d])

        # add attributes
        attrs = [a for a in dir(cur_group) if a not in ncdf_internals]
        for att in attrs:
            inclass.params[att] = getattr(cur_group, att)

    return f

"""
#----------------------------------------------------------------------------

def idl2ncdf(infile, outfile, zcut=None, desc=None,comp=False):
    
    import idlsave
    from tt.lines import rhpy
    
    if type(infile) != type([]):
        infile = [infile]
        
    nt = len(infile)
    first = True
    #desc  = []

    
    
    for f in range(nt):
        fin = idlsave.read(infile[f], verbose=False)
    
        # put arrays in C order
        for k in fin.keys():
            if k not in ['desc']:
                fin[k] = N.transpose(fin[k])
                
        if first:
            temp = N.zeros((nt,)+ fin['temp'].shape,dtype='f')
            vz   = N.zeros((nt,)+ fin['temp'].shape,dtype='f')
            nh   = N.zeros((nt,)+ fin['nh'].shape,dtype='f')
            nne  = N.zeros((nt,)+ fin['temp'].shape,dtype='d')
            z    = N.zeros((nt,)+ fin['z'].shape,dtype='f')
            Bx   = N.zeros((nt,)+ fin['temp'].shape,dtype='f')
            By   = N.zeros((nt,)+ fin['temp'].shape,dtype='f')
            Bz   = N.zeros((nt,)+ fin['temp'].shape,dtype='f')
            first = False
    
        if temp[0].shape != fin['temp'].shape:
            raise ValueError('(EEE) idl2ncdf: wrong shape of array read from %s: expected %s, got %s' \
            %  (infile[f], str(temp[0].shape), str(fin['temp'].shape)))
            
        temp[f] = fin['temp']
        vz[f]   = fin['vz']
        nh[f]   = fin['nh']
        nne[f]  = fin['nne']
        z[f]    = fin['z']
        Bx[f]   = fin['bx']
        By[f]   = fin['by']
        Bz[f]   = fin['bz']       
        

    if desc is None:
        desc = fin['desc']

    make_ncdf_atmos(outfile, temp, vz, nne, nh, z, Bx=Bx, By=By, Bz=Bz, x=fin['x'],
                    y=fin['y'], desc=desc,comp=comp)
    
    
    
    print("Wrote %s" % outfile)
    return

#----------------------------------------------------------------------------



def make_ncdf_atmos(outfile, T, vz, ne, nH, z, x=None, y=None, Bz=None, By=None,
                    Bx=None,desc=None, comp=False, append=True, snap=None,complev=2):
    ''' Creates NetCDF input file for rh15d.
    
        IN:
           outfile:  string name of destination. If file exists it will be wiped.
           T:        Temperature. Its shape will determine the dimensionality written
           vz:       Same shape as T. In m/s.
           ne:       Same shape as T. In m-3.
           nH:       Shape [6, shape.T]. In m-3.
           z:        Same shape as last index of T. In m.
           x:        Same shape as first index of T. In m.
           y:        Same shape as second index of T. In m.
           snap:     Snapshot number(s)
           Bx, By, Bz: same shape as T. In T.
           append:   if True, will append to existing file (if any).
           comp:     if false, compress file.
           complev:  compression level.
    
    '''
    
    import os
    
    mode = ['w','a']
    if (append and not os.path.isfile(outfile)): append=False
    

    rootgrp = nc.Dataset(outfile, mode[append], format='NETCDF4')
    
    complev = 2
    nt = 1
    
    if len(T.shape) == 1:
        x = [0.]
        y = [0.]
        # for these, only single snapshot is supported
        T  =  T[N.newaxis,N.newaxis,N.newaxis,:]
        ne = ne[N.newaxis,N.newaxis,N.newaxis,:]
        nH = nH[N.newaxis,N.newaxis,N.newaxis,:]
        vz = vz[N.newaxis,N.newaxis,N.newaxis,:]
        z  =  z[N.newaxis,N.newaxis,N.newaxis,:]
        if Bz is not None:
            Bx = Bx[N.newaxis,N.newaxis,N.newaxis,:]
            By = By[N.newaxis,N.newaxis,N.newaxis,:]
            Bz = Bz[N.newaxis,N.newaxis,N.newaxis,:]
    if len(T.shape) == 3: # single snapshot
        T  =  T[N.newaxis,:]
        ne = ne[N.newaxis,:]
        nH = nH[N.newaxis,:]
        vz = vz[N.newaxis,:]
        z  =  z[N.newaxis,:]
        if Bz is not None:
            Bx = Bx[N.newaxis,:]
            By = By[N.newaxis,:]
            Bz = Bz[N.newaxis,:]
    elif len(T.shape) != 4:
        raise ValueError('Invalid shape for T')
    else:
        nt = T.shape[0]
        nx = T.shape[1]
        ny = T.shape[2]
        
    if nH.shape == T.shape:
        nhydr = 1
    else:
        nhydr = nH.shape[0]
        
        
    if snap is None: snap = N.arange(nt,dtype='i4')
        
    # for a new file, create dimensions and variables
    if not append:
        rootgrp.createDimension('nt', None) # create unlimited dimension
        rootgrp.createDimension('nx', T.shape[-3])
        rootgrp.createDimension('ny', T.shape[-2])
        rootgrp.createDimension('nz', T.shape[-1])
        rootgrp.createDimension('nhydr', nhydr)

        T_var =rootgrp.createVariable('temperature', 'f4', ('nt','nx','ny','nz'),zlib=comp,
                                        least_significant_digit=1,complevel=complev)
        vz_var=rootgrp.createVariable('velocity_z', 'f4', ('nt','nx','ny','nz'),zlib=comp,
                                        least_significant_digit=1,complevel=complev)
        ne_var=rootgrp.createVariable('electron_density', 'f8', ('nt','nx','ny','nz'),
                                        zlib=comp, complevel=complev)
        nh_var=rootgrp.createVariable('hydrogen_populations', 'f4', 
                                        ('nt','nhydr','nx','ny','nz'),
                                        zlib=comp, complevel=complev) 
    
        x_var  = rootgrp.createVariable('x', 'f4', ('nx',))
        y_var  = rootgrp.createVariable('y', 'f4', ('ny',))
        z_var  = rootgrp.createVariable('z', 'f4', ('nt','nz'))
        nt_var = rootgrp.createVariable('snapshot_number', 'i4', ('nt',))
    
        if Bz is not None:
            bx_var  = rootgrp.createVariable('B_x', 'f4', ('nt','nx','ny','nz'),zlib=comp,
                                             least_significant_digit=5, complevel=complev)
            by_var  = rootgrp.createVariable('B_y', 'f4', ('nt','nx','ny','nz'),zlib=comp,
                                             least_significant_digit=5, complevel=complev)
            bz_var  = rootgrp.createVariable('B_z', 'f4', ('nt','nx','ny','nz'),zlib=comp,
                                             least_significant_digit=5, complevel=complev)  
   
        if desc is None:
            rootgrp.description = \
              "Juan BIFROST snapshot 20 from file qsmag-by00_t020.snap, zslice=[156,-30]"
        else:
            rootgrp.description = desc
        
        if Bz is None:
            rootgrp.has_B = 0
        else:
            rootgrp.has_B = 1
            
        nt = [0, nt]
    else:
        # get variables
        T_var = rootgrp.variables['temperature']
        vz_var= rootgrp.variables['velocity_z']
        ne_var= rootgrp.variables['electron_density']
        nh_var= rootgrp.variables['hydrogen_populations']
        nt_var= rootgrp.variables['snapshot_number']
        x_var = rootgrp.variables['x']
        y_var = rootgrp.variables['y']
        z_var = rootgrp.variables['z']
        
        if Bz is not None:
            bx_var = rootgrp.variables['B_x']
            by_var = rootgrp.variables['B_y']
            bz_var = rootgrp.variables['B_z']
        
        nti = len(rootgrp.dimensions['nt'])
        nt  = [nti, nti+nt]

    
    T_var[nt[0]:nt[1]] = T
    vz_var[nt[0]:nt[1]] = vz
    ne_var[nt[0]:nt[1]] = ne
    nh_var[nt[0]:nt[1],:nhydr] = nH
        
    if Bz is not None:
        bx_var[nt[0]:nt[1]] = Bx
        by_var[nt[0]:nt[1]] = By
        bz_var[nt[0]:nt[1]] = Bz
        
    x_var[:] = x
    y_var[:] = y
    z_var[nt[0]:nt[1]] = z
    nt_var[nt[0]:nt[1]] = snap

        
        
    
    rootgrp.close()
    
    
    
    return

#----------------------------------------------------------------------------

def depth_optim(height, temp, ne, vz, rho, nh=None, bx=None, by=None, bz=None,
                tmax=5e4):
    ''' Performs depth optimisation of one single column (as per multi_3d).
    
        IN:
            height   [cm]
            temp     [K]
            ne       [cm-3]
            vz       [any]
            rho      [g cm-3]
            nh       [any] (optional)
            bx,by,bz [any] (optional)
            tmax     [K] maximum temperature of the first point                            
            
    '''
    
    from scipy.integrate import cumtrapz
    import scipy.interpolate as interp
    
    
    ndep = len(height)
    
    # calculate optical depth from H-bf only
    taumax=100 
    grph=2.26e-24
    crhmbf=2.9256e-17
    ee = 1.602189E-12
    bk = 1.380662E-16
      
    xhbf=1.03526e-16*ne*crhmbf/temp**1.5*N.exp(0.754*ee/bk/temp)*rho/grph
    tau = N.concatenate(([0.], cumtrapz(xhbf,-height)))
    
    idx = (tau < taumax) & (temp < tmax)
        
    # find maximum variance of T, rho, and tau for each depth
    tt = temp[idx]
    rr = rho[idx]
    ta = tau[idx]
    
    tdiv   = N.abs(N.log10(tt[1:])-N.log10(tt[:-1]))/N.log10(1.1)
    rdiv   = N.abs(N.log10(rr[1:])-N.log10(rr[:-1]))/N.log10(1.1)
    taudiv = N.abs(N.log10(ta[1:])-N.log10(ta[:-1]))/0.1
    taudiv[0] = 0.
    
    aind   = N.concatenate(([0.],N.cumsum(N.max(N.array([tdiv,rdiv,taudiv]),axis=0))))
    aind  *= (ndep-1)/aind[-1]
        
    # interpolate new height so it is constant in aind2
    nheight = interp.splev(N.arange(ndep),interp.splrep(aind,height[idx],k=3,s=0),der=0)
    
    # interpolate quantities for new depth scale
    ntemp = N.exp(interp.splev(nheight,interp.splrep(height[::-1],N.log(temp[::-1]),
                                                     k=3,s=0),der=0))
    nne   = N.exp(interp.splev(nheight,interp.splrep(height[::-1],N.log(ne[::-1]),
                                                     k=3,s=0),der=0))
    nrho  = N.exp(interp.splev(nheight,interp.splrep(height[::-1],N.log(rho[::-1]),
                                                     k=3,s=0),der=0))
    nvz   =       interp.splev(nheight,interp.splrep(height[::-1], vz[::-1],
                                                     k=3,s=0),der=0)
    
    result = [nheight, ntemp, nne, nvz, nrho]
    
    if nh is not None:
        for k in range(nh.shape[0]):
            nh[k] = N.exp(interp.splev(nheight,interp.splrep(height[::-1],
                                                N.log(nh[k,::-1]),k=3,s=0),der=0))
        result += [nh]
    
    if bx is not None:
        nbx = interp.splev(nheight,interp.splrep(height[::-1],bx[::-1],k=3,s=0),der=0)
        nby = interp.splev(nheight,interp.splrep(height[::-1],by[::-1],k=3,s=0),der=0)
        nbz = interp.splev(nheight,interp.splrep(height[::-1],bz[::-1],k=3,s=0),der=0)
        result += [nbx, nby, nbz]


    return result

#----------------------------------------------------------------------------

def atmosncdf2multi(infile, xi, yi, outfile, writeB=False,zcut=0, write_dscale=False):
    ''' Writes MULTI atmosphere file from a column of the 3D model, in RH 1.5D ncdf 
        format. Also writes the binary XDR file with magnetic fields, if writeB
        is true. 

        The ncdf infile MUST have nH with 6 levels.
        
        OBSOLETE, DOESN'T WORK WITH NEW VERSION WITH MANY SNAPSHOTS!!!

        --Tiago, 20110209
        '''

    from tt.lines.multi import watmos_multi
    from tt.lines.rh    import write_B
    
    indata = DataHolder()
    inbuf  = read_ncdf(indata,infile)

    M_TO_CM3 = (100.)**3
    M_TO_KM = 0.001
    
    

    watmos_multi(outfile, indata.temperature[xi,yi,zcut:],
                          indata.electron_density[xi,yi,zcut:]/M_TO_CM3,
                          z  = indata.z[zcut:]*M_TO_KM,
                          vz = indata.velocity_z[xi,yi,zcut:]*M_TO_KM,
                          nh = indata.hydrogen_populations[:,xi,yi,zcut:]/M_TO_CM3,
                          id = '%s xy-slice: (x,y) = (%i,%i)' % 
                                   (indata.params['description'],xi, yi),
                          write_dscale=write_dscale)
    
    
    if writeB and indata.params['has_B']:
        write_B('%s.B' % outfile, indata.B_x[xi,yi,zcut:], indata.B_y[xi,yi,zcut:],
                indata.B_z[xi,yi,zcut:])
        print('*** Wrote %s.B' % outfile)

    inbuf.close()

    return

#-----------------------------------------------------------------------------------------

def moviefromrh(rr, width=700, sdir='.', scale=True, cmap='gist_heat', iwave=[63,341],
                vmax=None, vmin=None, Mm=None, mspec=None, iris_filter=False, suffix='',
                quant='intensity', xline=None, spec_filter=False, sqr=False):
    ''' Saves sequence of PNG images from Rh15d instance, with a ray loaded.
        There are two modes. In one, only one image is saved from an individual wavelength
        (can be used for sequentially for several snapshots). For the other mode,
        images are saved for each wavelength.
        
        IN:
          rr    - Rh15d instance, with ray loaded
          width - width of the resulting images, in pixels.
          sdir  - directory where to save the images
          scale - if true, will print a 1 Mm scale on the top right
        
        '''
        
    import pylab, math
    from tt.math.fitting import gaussian
    
    nx, ny = rr.ray.intensity.shape[:-1]
    sn = rr.ray.params['snapshot_number']
    
    pylab.figure(figsize=(8, 8.*ny/nx))
    pylab.axis('off')
    pylab.subplots_adjust(left=0,bottom=0,right=1.,top=1.)
    dpi = width/8.
    lcc = 0.
    
    if nx > ny: lcc = 0.02
    
    if Mm is not None:
        Mm2pix = float(nx/Mm)
    
    if iris_filter:
        from tt.sim.synobs import img_conv
        
        if not (200 < N.mean(rr.ray.wavelength) < 300):
            raise ValueError('moviefromrh: rr instance does not seem to be of Mg II' +
                             'and IRIS filter is selected!')
        
        im = img_conv(rr.ray.intensity, rr.ray.wavelength, xMm=Mm, graph=False,
                      conv_type='IRIS_MGII_CORE')  
        im = im.T[N.newaxis]
        seq = [0]
        
    if spec_filter:
        # convolves with a Gaussian in the spectral direction, FWHM = spec_filter,
        # central wavelength given by iwave
        # selecting wavelengths within 4 * FWHM
        wfwhm = spec_filter
        wave  = rr.ray.wavelength[:]
        wcent = wave[iwave]
        widx  = (wave[:] > wcent - 2.*wfwhm) & (wave[:] < wcent + 2.*wfwhm)
        # filtering function, here set to Gaussian
        wfilt = gaussian([wcent, wfwhm/(2*math.sqrt(2*math.log(2))), 1., 0.], wave[widx])
        wfilt /= N.trapz(wfilt, x=wave[widx])
        
        im = N.trapz(getattr(rr.ray, quant)[:,:,widx] * wfilt, x=wave[widx], axis=-1).T[N.newaxis]
        seq = [iwave]
        
        
    elif type(iwave) == type([]):
        seq = range(iwave[0], iwave[1])
        im = getattr(rr.ray,quant)[:,:,iwave[0]:iwave[1]].T
        #im = rr.ray.intensity[:,:,iwave[0]:iwave[1]].T
        
        if mspec is not None:
            mspec = N.mean(N.mean(im,axis=-1),axis=-1)
            # scale mspec and nwave to fit in small box
            bsize = [nx*0.2, ny*0.2*0.6]
            mspec *= bsize[1]/N.max(mspec)
            mwave = rr.ray.wavelength[iwave[0]:iwave[1]]
            mwave = bsize[0]*(mwave - mwave[0])/(mwave[-1]-mwave[0])
            
            # add offset
            mspec += 10
            mwave += 10
            
    else:
        seq = [iwave]
        im = getattr(rr.ray, quant)[:,:,iwave].T
        #im = rr.ray.intensity[:,:,iwave].T
        im = im[N.newaxis]
        mspec = None
    
    if vmax is None:
        vmax = N.max(im)
    if vmin is None:
        vmin = N.min(im)
        

        
    
    for i, iw in enumerate(seq):
        #vmax = N.max(im[i])/3.
        
        if mspec is not None:
            pylab.plot(mwave, mspec, 'b-', lw=2)
            pylab.plot([mwave[i]], [mspec[i]],'ro',mec='w')
        
        if scale:
            # put 1 Mm scale
            Mm = Mm2pix / nx # this is in relative coordinates
            pylab.plot([0.92-Mm, 0.92], [0.92-lcc, 0.92-lcc], 'w-', lw=3,
                       transform = pylab.gca().transAxes)
            pylab.text((0.92*2-Mm)/2., 0.95, '1 Mm', color='w', ha='center', va='center',
                       fontsize='large', transform = pylab.gca().transAxes)

        if xline is not None:
            pylab.axvline(x=xline, color='y', lw=2)
            suffix +='_x%04i' % int(xline)
        
        if iris_filter:
            pylab.imshow(im[i], cmap=cmap, vmax=vmax, vmin=vmin, interpolation='nearest')
        else:
            if sqr:
                pylab.imshow(N.sqrt(im[i]), cmap=cmap, vmax=N.sqrt(vmax), vmin=N.sqrt(vmin))
            else:
                pylab.imshow(im[i], cmap=cmap, vmax=vmax, vmin=vmin)

        pylab.savefig('%s/%s_s%03i_w%04i.png' % (sdir, suffix, sn, iw), dpi=dpi)
        pylab.clf()
    
    return
"""
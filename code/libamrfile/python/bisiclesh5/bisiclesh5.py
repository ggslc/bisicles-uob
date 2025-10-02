import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt

# code to read bisicles hdf5 file

# you will need install h5py which may be troublesome in condo but works using pip

# import bisiclesh5 as bh5
# a = bh5.bisicles_hf5('plot.stable_shelves.500m.000000.2d.hdf5','thickness')

# read specified variable from file. This creates a bisicles_hd5 class that contains

# bisicles_hf5.fname - filename
# bisicles_hf5.time
# bisicles_hf5.variable_name
# bisicles_hf5.full_name
# bisicles_hf5.units
# bisicles_hf5.num_levels - number of levels
# bisicles_hf5.offset - grid offsets for each level [level]
# bisicles_hf5.num_boxes - list of number of boxes in each level [level]
# bisicles_hf5.dx - grid spacing for each level [level]
# bisicles_hf5.x - x, y for each level in global coordinate system
# bisicles_hf5.y - as list of 1-d arrays [level][box]
# bisicles_hf5.data - data for each level as list of 2-d arrays [level][box]

# also methods
# bisicles_hf5.plot() - produces a nice plot with AMR boxes
# bisicles_hf5.mesh() - displays grid points (will be slow and need to zoom to see details)
# bisicles_hf5.flatten() - creates flattened version as flat class

# b = a.flatten()

# can control level that you want to flatten to (default is finest) and region to flatten
# (default is entire extent of coarest level. These options have not be tested much.

#  flat.fname
#  flat.time
#  flat.variable_name
#  flat.full_name
#  flat.units
#  flat.dx
#  flat.x
#  flat.y
#  flat.data

# there is a guick plot method as well

# b.plot()

class bisicles_hf5:

  def __init__(self,fname,target):

    def finder(tarvarstr,h5file,n_components):

      # find location of selected variable in hf5 file

      count = 0

      while count < n_components and h5file.attrs['component_'+str(count)].decode('utf-8') != tarvarstr:
        count += 1

      if count == n_components:
        count = np.NaN
        print('string not found')

      return count

    def variable_table(vname):

      # contains table to attach full name and units to a variable

      table = np.array([['thickness','Ice thickness','m'], \
                       ['xVel','Horizontal velocity x-component','m/yr'], \
                       ['yVel','Horizontal velocity y-component','m/yr'], \
                       ['vVel','Vertical velocity','m/yr'], \
                       ['Z_surface','Upper surface elevation','m'], \
                       ['Z_bottom','Lower surface elevation','m'], \
                       ['Z_base','Bedrock elevation','m'], \
                       ['basal_friction','Basal friction','Pa'], \
                       ['C0','-','-'], \
                       ['xRhs','Gravitational driving stress in x','Pa'], \
                       ['yRhs','Gravitational driving stress in y','Pa'], \
                       ['dThickness/dt','Rate of thickness change','m/yr'], \
                       ['xfVel','Horizontal velocity x-component at upper surface','m/yr'], \
                       ['yfVel','Horizontal velocity y-component at upper surface','m/yr'], \
                       ['zfVel','Vertical velocity at upper surface','m/yr'], \
                       ['xbVel','Horizontal velocity x-component at lower surface','m/yr'], \
                       ['ybVel','Horizontal velocity y-component at lower surface','m/yr'], \
                       ['zbVel','Vertical velocity at lower surface','m/yr'], \
                       ['dragCoef','Basal friction coefficient','(units depend on law)'], \
                       ['viscosityCoef','viscosityCoef','-'], \
                       ['xxViscousTensor','xxViscousTensor','Pa'], \
                       ['yxViscousTensor','yxViscousTensor','Pa'], \
                       ['xyViscousTensor','xyViscousTensor','Pa'], \
                       ['yyViscousTensor','yyViscousTensor','Pa'], \
                       ['activeBasalThicknessSource','Mass balance at lower surface (active)','m/yr'], \
                       ['activeSurfaceThicknessSource','Mass balance at upper surface (active)','m/yr'], \
                       ['divergenceThicknessFlux','Divergence of horizontal flux','m/yr'], \
                       ['basalThicknessSource','Mass balance at lower surface','m/yr'], \
                       ['surfaceThicknessSource','Mass balance at upper surface','m/yr'], \
                       ['calvingFlux','Calving flux','-'], \
                       ['mask','mask','-'] ])

      i = np.where(table[:,0] == vname)

      return table[i,1].item(), table[i,2].item()

    # open the hf5 file, read number of variables that it contains and number of levels
    h5file = h5py.File(fname,'r')

    n_components = h5file.attrs['num_components']
    # print(n_components)

    tarcomponent = finder(target,h5file,n_components)

    n_level = h5file.attrs['num_levels']

    # create stubs for various quantities that will be read these will be
    # constructed using .append
    levelsdx = []
    levelsboxes = []
    levelsoff = []
    levelsx = []
    levelsy = []
    levelsdata = []

    # crrate a time variable in case it isnt picked up from level info
    time = np.NAN

    # visit each level in turn coarse to fine
    for level in range(n_level):

      h5level = h5file['/level_'+str(level)+'/']

      # h5level.attrs.keys()
      # ['dt', 'dx', 'prob_domain', 'ref_ratio', 'time']

      # read this levels grid spacing and calculate its offset
      dx = h5level.attrs['dx']

      # this is the offset in x,y required to align grids at different levels
      if level == 0:
        offset = dx / 2.0
      else:
        offset = offset - dx / 2.0

      # print(dx,offset)

      # read model time if it is there (years)
      if h5level.attrs.__contains__('time'):
        time = h5level.attrs['time']

      # read data (h5data) for this level, which is found in boxes so need to read
      # pointer into data for each box (h5offs)

      # print(h5level.keys())
      # ['Processors', 'boxes', 'data:datatype=0', 'data:offsets=0', 'data_attributes']
      h5box = h5level.get('boxes')
      h5box = np.array(h5box)
      h5data = h5level.get('data:datatype=0')
      h5data = np.array(h5data)
      h5offs = h5level.get('data:offsets=0')
      h5offs = np.array(h5offs)

      n_boxes = len(h5box)

      # create stubs for quantities (x,y,data) held in boxes on this level
      boxesx = []
      boxesy = []
      boxesdata = []

      # visit each box in this level
      for box in range(n_boxes):

        # find x, y for this box from bounding box information remembering
        # level offsets in x,y and add border of ghost cells
        # NOT clear why +2 but provides correct numbers when compared to h5offs
        y = np.arange(h5box['lo_i'][box]-1,h5box['hi_i'][box]+2) * dx + offset
        x = np.arange(h5box['lo_j'][box]-1,h5box['hi_j'][box]+2) * dx + offset

        # add x, y data to evolving level
        boxesx.append(x)
        boxesy.append(y)

        # box size
        nx = len(x)
        ny = len(y)

        # find locations of box data in h5
        en = h5offs[box]

        # if level < 2:
        #   print(level, box, h5offs[box], h5box['lo_i'][box], h5box['hi_i'][box], h5box['lo_j'][box], h5box['hi_j'][box], nx, ny)

        st = en + tarcomponent * nx * ny
        en = st + nx * ny

        boxesdata.append(h5data[st:en].reshape((nx,ny)))

      # once completed all boxes in a level add them to the overall structure
      levelsx.append(boxesx)
      levelsy.append(boxesy)
      levelsdata.append(boxesdata)

      # also add useful level info such as grid spacing, grid offset and number boxes
      levelsdx.append(dx)
      levelsoff.append(offset)
      levelsboxes.append(n_boxes)

    h5file.close

    # read variable name and use table to attach useful info such a long name and units
    vname = h5file.attrs['component_'+str(tarcomponent)].decode('utf-8')
    fullname, units = variable_table(vname)

    # pack all of the information into the class
    self.fname = fname
    self.time = time
    self.variable_name = vname
    self.full_name = fullname
    self.units = units
    self.num_levels = n_level # number of levels
    self.offset = levelsoff # grid offsets for each level
    self.num_boxes = levelsboxes # list of number of boxes in each level
    self.dx = levelsdx # grid spacing for each level
    self.x = levelsx # x, y for each level in global coordinate system
    self.y = levelsy # as list of arrays [level][box]
    self.data = levelsdata # data for each level as list of arrays [level][box]

  def flatten(self,lev=-1,xmin=np.NAN,xmax=np.NAN,ymin=np.NAN,ymax=np.NAN):

    class flat:

      def __init__(self,fname,time,variable_name,full_name,units,dx,x,y,data):

        self.fname = fname
        self.time = time
        self.variable_name = variable_name
        self.full_name = full_name
        self.units = units
        self.dx = dx
        self.x = x
        self.y = y
        self.data = data

      def plot(self):

        fig, ax = plt.subplots()

        X, Y = np.meshgrid(self.y - self.dx/2.0,self.x - self.dx/2.0)

        pcm = ax.pcolormesh(X,Y,self.data,shading='nearest',cmap='hsv')

        ax.set_aspect('equal','box')

        fig.colorbar(pcm,ax = ax,location = 'right',label = self.units)
        plt.title(self.full_name + ' in ' + self.units + ' at ' + '{:.2f}'.format(self.time) + ' years')

        fig.tight_layout()
        plt.show()

    def findend(xx,x,i1,i2):

      value = np.where((xx >= x[i1]) & (xx <= x[i2]))

      return value[0][0]

    # if bounding box not specfied then
    # find based on x,y bounds of the coarsest grid

    dx0 = self.dx[0]/2.0
    xmin = np.min(self.x[0][:]) + dx0
    xmax = np.max(self.x[0][:]) - dx0
    ymin = np.min(self.y[0][:]) + dx0
    ymax = np.max(self.y[0][:]) - dx0

    print(xmin,xmax,ymin,ymax)

    # create x,y coordinates for base grid at resolution specfied in call (default is finest)
    # xx = np.arange(xmin,xmax+self.dx[lev],self.dx[lev]) + self.offset[lev]
    # yy = np.arange(ymin,ymax+self.dx[lev],self.dx[lev]) + self.offset[lev]

    xx = np.arange(xmin+0.5*self.dx[lev],xmax-1.5*self.dx[lev],self.dx[lev])
    yy = np.arange(ymin+0.5*self.dx[lev],ymax-1.5*self.dx[lev],self.dx[lev])

    # if bounds have been given find nearest poinst on base grid
    if not np.isnan(xmin):
      xmin = xx[np.nonzero(xx>xmin)[0][[0]]]

    if not np.isnan(xmax):
      xmax = xx[np.nonzero(xx<xmax)[0][[-1]]]

    if not np.isnan(ymin):
      ymin = yy[np.nonzero(yy>ymin)[0][[0]]]

    if not np.isnan(ymax):
      ymax = yy[np.nonzero(yy<ymax)[0][[-1]]]

    print(xmin,xmax,ymin,ymax)
    print(xx[0],xx[1],xx[-2],xx[-1])
    print(yy[0],yy[1],yy[-2],yy[-1])

    # create suitably sized array to hold flattened data (base grid)
    data = np.zeros((len(xx),len(yy)))

    # visit each level (coarest downwards) and box
    for level in range(self.num_levels):
      for box in range(self.num_boxes[level]):

        # print(level,box)

        # find the i,j corrdinates of the bounds for this box
        # on the new base grid

        imin = findend(xx,self.x[level][box],0,1)
        imax = findend(xx,self.x[level][box],-2,-1) + 2

        jmin = findend(yy,self.y[level][box],0,1)
        jmax = findend(yy,self.y[level][box],-2,-1) + 2

        # create and interpolator for this box
        # including any ghost cells
        interp = RegularGridInterpolator((self.x[level][box],self.y[level][box]), \
                                         self.data[level][box], \
                                         method='nearest',bounds_error=False,fill_value=None)

        # create x,y coordinates for poiints on new base grid that fall within this box
        # and are to be filled
        xm,ym = np.meshgrid(yy[jmin:jmax],xx[imin:imax])

        # print(imin,imax,jmin,jmax)
        # print(np.shape(yy[jmin:jmax]),np.shape(xx[imin:imax]))
        # print(np.shape(xm),np.shape(ym),np.shape(data[imin:imax,jmin:jmax]))

        # do the nearest neighbor interpolation
        data[imin:imax,jmin:jmax] = interp((ym,xm))

        # create a class to hold flattened version and return it

    return flat(self.fname,self.time,self.variable_name,self.full_name,self.units,self.dx[lev],xx,yy,data)

  def plot(self,datamin=np.NAN,datamax=np.NAN,levelmax=999):

    def minmax(data):

      # find min and max values acorss all levels and all boxes

      levmax = []
      levmin = []

      for level in range(len(data)):

        boxmax = []
        boxmin = []

        for box in range(len(data[level])):

          boxmax.append(np.nanmax(data[level][box]))
          boxmin.append(np.nanmin(data[level][box]))

        levmax.append(max(boxmax))
        levmin.append(min(boxmin))

      datamax = max(levmax)
      datamin = min(levmin)

      return datamin, datamax

    # control number of levels plotted

    select_level = np.min([levelmax,self.num_levels])

    if np.isnan(datamin) | np.isnan(datamax):
      datamin, datamax = minmax(self.data)

    fig, ax = plt.subplots()

    color = np.linspace(0.0,1.0,self.num_levels)

    for level in range(select_level):
      for box in range(self.num_boxes[level]):

        offset = self.dx[level]/2.0

        # shifting everything by -dx/2 so that plotting become cell centred
        X, Y = np.meshgrid(self.y[level][box] - offset,self.x[level][box] - offset)

        pcm = ax.pcolormesh(X,Y,self.data[level][box],shading = 'nearest',cmap = 'hsv', \
                            vmin = datamin,vmax = datamax,zorder = 0)

        # draw bounding boxes omitting ghost cells
        # zorder needed to ensure boxes appear about colormesh
        rect = plt.Rectangle((self.y[level][box][0],self.x[level][box][0]) + offset, \
                              self.y[level][box][-1]-self.y[level][box][1], \
                              self.x[level][box][-1]-self.x[level][box][1], \
                              linewidth = 1, edgecolor = str(color[level]),facecolor = 'none',zorder = 10)

        ax.add_patch(rect)

    fig.colorbar(pcm,ax = ax, location = 'right', label = self.units)
    plt.title(self.full_name + ' in ' + self.units + ' at ' + '{:.2f}'.format(self.time) + ' years')

    # tidy up axis becasue of messing about shifting the grid above
    # ax.set_xlim([np.min(self.y[0])+self.dx[0]/2.0,np.max(self.y[0])-self.dx[0]/2.0])
    # ax.set_ylim([np.min(self.x[0])+self.dx[0]/2.0,np.max(self.x[0])-self.dx[0]/2.0])

    ax.set_aspect('equal','box')

    fig.tight_layout()
    plt.show()

  def mesh(self):

    fig, ax = plt.subplots()

    marker = np.arange(0,self.num_levels)
    color = np.linspace(0.0,1.0,self.num_levels+1)

    for level in range(self.num_levels):
      for box in range(self.num_boxes[level]):

        offset = self.dx[level]/2.0

        X, Y = np.meshgrid(self.y[level][box],self.x[level][box])

          # print(box, X.shape, self.x[level][box].shape)

        ax.scatter(X,Y,alpha = 0.5,marker = "$"+str(marker[level])+"$",color='k',zorder=0)

        rect = plt.Rectangle((self.y[level][box][0],self.x[level][box][0]) + offset, \
                              self.y[level][box][-1]-self.y[level][box][1], \
                              self.x[level][box][-1]-self.x[level][box][1], \
                              linewidth = 2, edgecolor = str(color[level]),facecolor = 'none',zorder = 10)

        ax.add_patch(rect)

    ax.set_aspect('equal','box')

    ax.set_xlim([np.min(self.y[0])+self.dx[0]/2.0,np.max(self.y[0])-self.dx[0]/2.0])
    ax.set_ylim([np.min(self.x[0])+self.dx[0]/2.0,np.max(self.x[0])-self.dx[0]/2.0])

    fig.tight_layout()
    plt.show()

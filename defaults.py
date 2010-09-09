from __future__ import division
from kspace import KSpace

default_kspace = KSpace(1)
default_resolution = 32
default_numbands = 2
default_initcode = ''
default_mode = 'tm'
isQuiet = False
template = '''%(initcode)s

(set! num-bands %(numbands)s)

(set! k-points %(kspace)s)

(set! geometry (list %(geometry)s))

(set! geometry-lattice %(lattice)s)

(set! resolution %(resolution)s)

(run-%(mode)s)'''

#~ epsilon_color_map = 'gray'
#~ 
#~ color_map_dir = '/usr/share/h5utils/colormaps'

fig_size = (6,6)
contour_lines = {'colors':'k','linestyles':['dashed','solid'], 
'linewidths':1.0}
contour_plain = {'linewidths':1.0}
contour_filled = {}
colorbar_style = {'extend':'both','shrink':0.8}

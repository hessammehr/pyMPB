    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the imKvecFormatterplied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
from kspace import KSpace

default_kspace = KSpace(1)
default_resolution = 32
default_numbands = 2
default_initcode = ''
default_postcode = ''
default_runcode = '(run-tm)'

default_number_of_tiles_to_output = 3
default_field_component_to_export = 'z'
#mpb_call = 'mpb'
mpb_call = 'mpirun -np %(num_procs)s mpbi-mpi'

mpbdata_call = ('mpb-data -rn%(resolution)s -x%(number_of_tiles_to_output)s '
                '-y%(number_of_tiles_to_output)s -o%(output_file)s '
                '%(h5_file)s')
epsh5topng_call_2D = 'h5topng -S3 -Zrcbluered -oepsilon.png %(h5_file)s'
epsh5topng_call_3D = 'h5topng -0z0 -S3 -Zrcbluered -oepsilon.png %(h5_file)s'
epsh5topng_call_3D_cross_sect = ('h5topng -0x0 -S3 -Zrcbluered ' 
                              '-oepsilonslab.png %(h5_file)s')
                              
fieldh5topng_call_2D = ('h5topng -S3 -Zcbluered -C%(eps_file)s '
                        '-o%(output_file)s %(h5_file)s')   
fieldh5topng_call_2D_no_ovl = ('h5topng -S3 -Zcbluered '
                        '-o%(output_file_no_ovl)s %(h5_file)s')   
                        
fieldh5topng_call_3D = ('h5topng -0z0 -S3 -Zcbluered -C%(eps_file)s '
                        '-o%(output_file)s %(h5_file)s')   
fieldh5topng_call_3D_no_ovl = ('h5topng -0z0 -S3 -Zcbluered '
                        '-o%(output_file_no_ovl)s %(h5_file)s')                           
                     
temporary_epsh5 = './temporary_eps.h5'                     
temporary_h5 = './temporary.h5'    
temporary_h5_folder = './patterns~/'                
                      
isQuiet = False
template = '''%(initcode)s

(set! geometry-lattice %(lattice)s)

(set! resolution %(resolution)s)

(set! num-bands %(numbands)s)

(set! k-points %(kspace)s)

(set! geometry (list %(geometry)s))

%(runcode)s

%(postcode)s

(display-eigensolver-stats)'''

#~ epsilon_color_map = 'gray'
#~ 
#~ color_map_dir = '/usr/share/h5utils/colormaps'

fig_size = (6,6)
contour_lines = {'colors':'k','linestyles':['dashed','solid'], 'linewidths':1.0}
contour_plain = {'linewidths':1.0}
contour_filled = {}
colorbar_style = {'extend':'both','shrink':0.8}

    #Copyright 2016 Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import sys
sys.path.append('../')
from os import path

import numpy as np
import matplotlib.pyplot as plt

from phc_simulations import TriHoles2D
import log
import defaults

def main():

    if len(sys.argv) > 1:
        mode=sys.argv[1]
    else:
        mode='sim'

    defaults.add_epsilon_as_inset = True

    sim = TriHoles2D(
        material='SiN',
        radius=0.34,
        numbands=4,#8,
        k_interpolation=5,#31, 
        resolution=16, 
        mesh_size=7,
        runmode=mode,
        num_processors=2,
        save_field_patterns=False,
        convert_field_patterns=False)
      
    if not sim:
        log.error('an error occurred during simulation. See the .out file')
        return
    
    log.info(' ##### success! #####\n\n')
        

if __name__ == '__main__':
    main()
